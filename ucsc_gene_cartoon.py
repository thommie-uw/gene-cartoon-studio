#!/usr/bin/env python3
"""
ucsc_gene_cartoon.py
====================

Publication-quality gene structure cartoons drawn from live UCSC Genome
Browser REST API data (https://genome.ucsc.edu/goldenPath/help/api.html).

Give it a gene symbol; it resolves the locus, pulls transcript models from a
UCSC annotation track, and renders a vector (SVG/PDF) or raster (PNG/TIFF)
figure suitable for a manuscript.

    python ucsc_gene_cartoon.py TP53 -o tp53.svg

Everything about the drawing -- colours, box heights, intron style, exon
numbering, fonts, isoform stacking, intron compression -- is controlled by a
JSON style file that you can dump, edit and reuse:

    python ucsc_gene_cartoon.py --dump-style my_style.json
    python ucsc_gene_cartoon.py CFTR --style my_style.json -o cftr.svg

Requires: requests, matplotlib.
"""

from __future__ import annotations

#: Bumped whenever app.py relies on something new in here.  app.py checks it
#: and refuses to run against a stale copy, because a half-updated pair of
#: files fails silently and confusingly (controls appear but do nothing).
__version__ = "1.6.0"

import argparse
import bisect
import hashlib
import io
import json
import math
import os
import re
import sys
import tempfile
import textwrap
import time
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

try:
    import requests
except ImportError:  # pragma: no cover
    sys.exit("This tool needs `requests`.  Install with:  pip install requests")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Polygon, Rectangle  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402


# --------------------------------------------------------------------------- #
#  UCSC REST API client
# --------------------------------------------------------------------------- #

API_ROOT = "https://api.genome.ucsc.edu"
BROWSER_ROOT = "https://genome.ucsc.edu/cgi-bin/hgTracks"

#: Annotation tracks that carry genePred-style transcript models, in the order
#: we prefer them.  ncbiRefSeqSelect == one MANE Select transcript per gene.
GENE_TRACKS = [
    "ncbiRefSeqSelect",   # MANE Select / RefSeq Select -- one clean model
    "ncbiRefSeqCurated",  # curated NM_/NR_ transcripts
    "ncbiRefSeq",         # everything incl. predicted XM_
    "refGene",            # legacy RefSeq
    "knownGene",          # GENCODE / UCSC genes (ENST ids)
    "wgEncodeGencodeBasicV44",
]

GENEPRED_TYPES = {"genePred", "genePredWithSomeExtraFields"}


class UCSCError(RuntimeError):
    pass


class UCSCClient:
    """Thin, polite, caching wrapper around the UCSC REST API."""

    def __init__(self, genome: str = "hg38", cache_dir: Optional[Path] = None,
                 timeout: int = 60, quiet: bool = False, use_cache: bool = True):
        self.genome = genome
        self.timeout = timeout
        self.quiet = quiet
        self.use_cache = use_cache
        self.cache_dir = self._resolve_cache_dir(cache_dir)
        self.session = requests.Session()
        self.session.headers["User-Agent"] = "ucsc_gene_cartoon/1.0"
        self._last_call = 0.0

    # -- plumbing ---------------------------------------------------------- #

    @staticmethod
    def _resolve_cache_dir(explicit: Optional[Path]) -> Optional[Path]:
        """
        Pick a writable cache directory, or None to run without a disk cache.

        Hosted environments (Streamlit Cloud, containers, locked-down lab
        machines) may have no home directory or a read-only one, so fall back
        to the temp dir and finally to no caching at all rather than crashing
        on import.
        """
        candidates = []
        if explicit:
            candidates.append(Path(explicit))
        env = os.environ.get("UCSC_CARTOON_CACHE")
        if env:
            candidates.append(Path(env))
        try:
            candidates.append(Path.home() / ".cache" / "ucsc_gene_cartoon")
        except (RuntimeError, OSError):
            pass
        candidates.append(Path(tempfile.gettempdir()) / "ucsc_gene_cartoon")

        for path in candidates:
            try:
                path.mkdir(parents=True, exist_ok=True)
                probe = path / ".write_test"
                probe.write_text("ok")
                probe.unlink()
                return path
            except (OSError, RuntimeError):
                continue
        return None      # memory only; st.cache_data still shields the API

    def _log(self, msg: str) -> None:
        if not self.quiet:
            print(f"[ucsc] {msg}", file=sys.stderr)

    def _get(self, endpoint: str, params: Dict[str, Any],
             use_cache: bool = True) -> Dict[str, Any]:
        use_cache = use_cache and self.use_cache and self.cache_dir is not None
        url = f"{API_ROOT}/{endpoint}"
        cache_file = None
        if use_cache:
            key = hashlib.sha1(
                json.dumps([url, sorted(params.items())]).encode()
            ).hexdigest()
            cache_file = self.cache_dir / f"{key}.json"
            if cache_file.exists():
                try:
                    return json.loads(cache_file.read_text())
                except (json.JSONDecodeError, OSError):
                    try:
                        cache_file.unlink()
                    except OSError:
                        pass

        # be a good citizen: UCSC asks for <1 request/sec sustained
        gap = time.monotonic() - self._last_call
        if gap < 0.34:
            time.sleep(0.34 - gap)

        self._log(f"GET {endpoint} {params}")
        try:
            r = self.session.get(url, params=params, timeout=self.timeout)
        except requests.RequestException as e:
            raise UCSCError(f"network error talking to UCSC: {e}") from e
        self._last_call = time.monotonic()

        if r.status_code != 200:
            raise UCSCError(
                f"UCSC returned HTTP {r.status_code} for {r.url}\n{r.text[:400]}"
            )
        try:
            data = r.json()
        except json.JSONDecodeError as e:
            raise UCSCError(f"UCSC returned non-JSON for {r.url}") from e
        if isinstance(data, dict) and "error" in data:
            raise UCSCError(f"UCSC API error: {data['error']}")

        if cache_file is not None:
            try:
                cache_file.write_text(json.dumps(data))
            except OSError:
                pass
        return data

    # -- endpoints --------------------------------------------------------- #

    def list_tracks(self) -> Dict[str, Any]:
        data = self._get("list/tracks", {"genome": self.genome})
        return data.get(self.genome, {})

    def gene_pred_tracks(self) -> List[str]:
        """Tracks in this assembly that hold transcript models we can draw."""
        out = []
        for name, meta in self.list_tracks().items():
            if isinstance(meta, dict) and meta.get("type", "").split()[0] in GENEPRED_TYPES:
                out.append(name)
        return sorted(out)

    def search(self, term: str) -> List[Tuple[str, str, str]]:
        """Return [(trackName, posName, 'chrN:start-end'), ...] for a query."""
        data = self._get("search", {"search": term, "genome": self.genome})
        hits: List[Tuple[str, str, str]] = []
        for block in data.get("positionMatches", []):
            track = block.get("trackName", "?")
            for m in block.get("matches", []):
                pos = m.get("position", "")
                if re.match(r"^\w+:\d+-\d+$", pos):
                    hits.append((track, m.get("posName", ""), pos))
        return hits

    def locate_symbol(self, symbol: str) -> Tuple[str, int, int]:
        """Resolve a gene symbol to a padded (chrom, start, end) search window."""
        hits = self.search(symbol)
        if not hits:
            raise UCSCError(
                f"UCSC found nothing for {symbol!r} in {self.genome}. "
                "Check the symbol and assembly, or pass --region chrN:start-end."
            )
        # Prefer authoritative, symbol-keyed tracks and an exact name match.
        priority = {"mane": 0, "hgnc": 1, "ncbiRefSeqSelect": 2,
                    "knownGene": 3, "refGene": 4}

        def rank(h: Tuple[str, str, str]) -> Tuple[int, int]:
            track, name, _ = h
            exact = 0 if symbol.upper() in name.upper() else 1
            return (exact, priority.get(track, 9))

        track, name, pos = sorted(hits, key=rank)[0]
        self._log(f"{symbol} -> {pos} (via {track}: {name})")
        chrom, span = pos.split(":")
        start, end = (int(x) for x in span.split("-"))
        pad = max(2000, int((end - start) * 0.1))
        return chrom, max(0, start - 1 - pad), end + pad

    def transcripts_in_region(self, track: str, chrom: str, start: int,
                              end: int) -> List["Transcript"]:
        data = self._get("getData/track", {
            "genome": self.genome, "track": track,
            "chrom": chrom, "start": start, "end": end,
        })
        rows = data.get(track)
        # Some tracks nest rows under the chromosome name.
        if isinstance(rows, dict):
            rows = rows.get(chrom, [])
        if not isinstance(rows, list):
            return []
        out = []
        for row in rows:
            try:
                out.append(Transcript.from_genepred(row))
            except (KeyError, TypeError, ValueError):
                continue
        return out

    def get_sequence(self, chrom: str, start: int, end: int,
                     rev_comp: bool = False) -> str:
        params = {"genome": self.genome, "chrom": chrom,
                  "start": start, "end": end}
        if rev_comp:
            params["revComp"] = 1
        return self._get("getData/sequence", params).get("dna", "").upper()

    def base_at(self, chrom: str, pos_1based: int) -> str:
        """The reference base at a 1-based genomic position."""
        return self.get_sequence(chrom, pos_1based - 1, pos_1based)

    # -- shareable links --------------------------------------------------- #

    def sequence_api_url(self, chrom: str, start: int, end: int) -> str:
        return (f"{API_ROOT}/getData/sequence?genome={self.genome};chrom={chrom};"
                f"start={start};end={end}")

    def browser_url(self, chrom: str, start: int, end: int) -> str:
        return (f"{BROWSER_ROOT}?db={self.genome}&position="
                f"{chrom}%3A{start + 1}-{end}")


# --------------------------------------------------------------------------- #
#  Transcript model
# --------------------------------------------------------------------------- #

@dataclass
class Transcript:
    name: str
    gene: str
    chrom: str
    strand: str
    tx_start: int          # 0-based, half-open (UCSC internal convention)
    tx_end: int
    cds_start: int
    cds_end: int
    exons: List[Tuple[int, int]]

    @classmethod
    def from_genepred(cls, row: Dict[str, Any]) -> "Transcript":
        starts = [int(x) for x in str(row["exonStarts"]).rstrip(",").split(",") if x != ""]
        ends = [int(x) for x in str(row["exonEnds"]).rstrip(",").split(",") if x != ""]
        if len(starts) != len(ends):
            raise ValueError("malformed exon arrays")
        tx_start = int(row.get("txStart", row.get("chromStart")))
        tx_end = int(row.get("txEnd", row.get("chromEnd")))
        cds_start = int(row.get("cdsStart", row.get("thickStart", tx_start)))
        cds_end = int(row.get("cdsEnd", row.get("thickEnd", tx_start)))
        return cls(
            name=str(row.get("name", "?")),
            gene=str(row.get("name2") or row.get("geneName") or row.get("name", "?")),
            chrom=str(row["chrom"]),
            strand=str(row.get("strand", "+")),
            tx_start=tx_start, tx_end=tx_end,
            cds_start=cds_start, cds_end=cds_end,
            exons=sorted(zip(starts, ends)),
        )

    # -- derived ----------------------------------------------------------- #

    @property
    def coding(self) -> bool:
        return self.cds_end > self.cds_start

    @property
    def length(self) -> int:
        return self.tx_end - self.tx_start

    @property
    def spliced_length(self) -> int:
        return sum(e - s for s, e in self.exons)

    @property
    def introns(self) -> List[Tuple[int, int]]:
        return [(self.exons[i][1], self.exons[i + 1][0])
                for i in range(len(self.exons) - 1)]

    def exon_number(self, index: int) -> int:
        """Biological exon number (5'->3'), 1-based."""
        return index + 1 if self.strand == "+" else len(self.exons) - index

    def exon_span(self, number: int) -> Optional[Tuple[int, int]]:
        """Genomic span of a biological exon number, or None if out of range."""
        if not 1 <= number <= len(self.exons):
            return None
        idx = number - 1 if self.strand == "+" else len(self.exons) - number
        return self.exons[idx]

    def exon_window(self, first: int, last: Optional[int] = None,
                    flank: int = 0) -> Tuple[int, int]:
        """
        Genomic window covering exon(s), plus optional flanking intron.

        Numbers are biological, so ``exon_window(14)`` means exon 14 as a
        reader would count it, on either strand.
        """
        last = first if last is None else last
        lo, hi = min(first, last), max(first, last)
        spans = [self.exon_span(n) for n in range(lo, hi + 1)]
        spans = [s for s in spans if s]
        if not spans:
            raise ValueError(
                f"{self.name} has exons 1-{len(self.exons)}; "
                f"asked for {first}" + (f"-{last}" if last != first else "")
            )
        g0 = min(s for s, _ in spans) - max(0, flank)
        g1 = max(e for _, e in spans) + max(0, flank)
        return max(self.tx_start - flank, g0), min(self.tx_end + flank, g1)

    def segments(self) -> List[Tuple[int, int, str]]:
        """Split exons into ('utr5' | 'cds' | 'utr3') pieces."""
        out: List[Tuple[int, int, str]] = []
        for s, e in self.exons:
            if not self.coding:
                out.append((s, e, "nc"))
                continue
            # left UTR piece
            if s < min(e, self.cds_start):
                out.append((s, min(e, self.cds_start),
                            "utr5" if self.strand == "+" else "utr3"))
            # CDS piece
            cs, ce = max(s, self.cds_start), min(e, self.cds_end)
            if cs < ce:
                out.append((cs, ce, "cds"))
            # right UTR piece
            if max(s, self.cds_end) < e:
                out.append((max(s, self.cds_end), e,
                            "utr3" if self.strand == "+" else "utr5"))
        return out


# --------------------------------------------------------------------------- #
#  cDNA (HGVS) <-> genomic coordinates
# --------------------------------------------------------------------------- #

class CDNAError(ValueError):
    pass


#: c.123  c.-45  c.*67  c.376-2  c.375+1  c.*12+3   (with optional prefix/suffix)
_HGVS_RE = re.compile(
    r"""
    ^\s*
    (?:[^:]*:)?                 # optional NM_000546.6:  or  TP53:
    \s*
    (?:(?P<kind>[cnrgmp])\.)?   # coordinate system prefix
    \s*
    (?P<star>\*)?               # 3' UTR marker
    (?P<minus>-)?               # 5' UTR marker
    (?P<base>\d+)               # the cDNA position itself
    (?:(?P<offsign>[+-])(?P<offset>\d+))?   # intronic offset
    (?P<rest>.*)$               # ref>alt, del, dup, ins ... ignored for placing
    """,
    re.VERBOSE,
)


@dataclass(frozen=True)
class CDNAPosition:
    """A parsed HGVS c. position."""
    base: int          # signed: negative = 5'UTR, positive = CDS
    utr3: bool = False  # True for c.*N
    offset: int = 0     # intronic offset, e.g. +1 / -2
    kind: str = "c"

    def __str__(self) -> str:
        core = f"*{self.base}" if self.utr3 else str(self.base)
        off = f"{self.offset:+d}" if self.offset else ""
        return f"{self.kind}.{core}{off}"


def parse_hgvs(text: Any) -> CDNAPosition:
    """
    Parse an HGVS-ish cDNA position.

    Accepts ``c.743G>A``, ``NM_000546.6:c.524G>A``, ``c.-124C>T``,
    ``c.*38A>G``, ``c.376-2A>G``, ``c.1521_1523del`` (uses the first position)
    and bare integers.  Protein (``p.``) and genomic (``g.``) notations are
    rejected -- they can't be placed on a cDNA axis.
    """
    if text is None or (isinstance(text, float) and text != text):
        raise CDNAError("no cDNA position given (the cell is empty)")
    s = str(text).strip()
    if not s:
        raise CDNAError("no cDNA position given (the cell is empty)")
    if re.match(r"^\s*(?:[^:]*:)?\s*p\.", s, re.I):
        raise CDNAError(
            f"{s!r} is protein (p.) notation; this tool places variants on a "
            "cDNA axis, so it needs c./n. positions such as c.743G>A"
        )
    m = _HGVS_RE.match(s)
    if not m:
        raise CDNAError(f"cannot parse {s!r} as a cDNA position")

    kind = (m.group("kind") or "c").lower()
    if kind in ("p", "g", "m"):
        raise CDNAError(
            f"{s!r} is {kind}. notation; this tool needs cDNA (c./n./r.) positions"
        )
    if kind == "r":
        kind = "c"

    base = int(m.group("base"))
    if m.group("minus"):
        base = -base
    offset = 0
    if m.group("offset"):
        offset = int(m.group("offset"))
        if m.group("offsign") == "-":
            offset = -offset
    return CDNAPosition(base=base, utr3=bool(m.group("star")),
                        offset=offset, kind=kind)


class CDNAMapper:
    """
    Convert cDNA (HGVS ``c.``) positions to genomic coordinates for one
    transcript, and back again.

    Numbering follows the HGVS recommendations:

    * ``c.1`` is the A of the initiator ATG, ``c.-1`` the base before it
    * ``c.*1`` is the base after the stop codon
    * ``c.375+1`` / ``c.376-2`` are intronic, offset from the flanking exon base
    * ``n.`` numbering (non-coding transcripts) counts from the transcript 5' end

    All genomic coordinates in and out of this class are **1-based**, matching
    what the UCSC browser and HGVS both display.
    """

    def __init__(self, tx: Transcript):
        self.tx = tx
        self.strand = tx.strand
        # exon blocks in transcript order, as (tx_lo, tx_hi, g_first, step)
        self._blocks: List[Tuple[int, int, int, int]] = []
        exons = tx.exons if tx.strand == "+" else list(reversed(tx.exons))
        pos = 1
        for s, e in exons:
            length = e - s
            if tx.strand == "+":
                g_first, step = s + 1, 1
            else:
                g_first, step = e, -1
            self._blocks.append((pos, pos + length - 1, g_first, step))
            pos += length
        self.tx_length = pos - 1

        if tx.coding:
            g_cds_start = tx.cds_start + 1 if tx.strand == "+" else tx.cds_end
            g_cds_end = tx.cds_end if tx.strand == "+" else tx.cds_start + 1
            self.n_cds_start = self.genomic_to_tx(g_cds_start)
            self.n_cds_end = self.genomic_to_tx(g_cds_end)
        else:
            # n. numbering: c.1 == the first transcribed base
            self.n_cds_start = 1
            self.n_cds_end = self.tx_length
        if self.n_cds_start is None or self.n_cds_end is None:
            raise CDNAError(f"{tx.name}: CDS boundaries fall outside its exons")

    # -- primitives -------------------------------------------------------- #

    def tx_to_genomic(self, n: int) -> Optional[int]:
        """Transcript coordinate (1-based from the 5' end) -> genomic."""
        for lo, hi, g0, step in self._blocks:
            if lo <= n <= hi:
                return g0 + step * (n - lo)
        return None

    def genomic_to_tx(self, g: int) -> Optional[int]:
        """Genomic (1-based) -> transcript coordinate, or None if intronic."""
        for lo, hi, g0, step in self._blocks:
            n = lo + step * (g - g0)
            if lo <= n <= hi:
                return n
        return None

    # -- the useful direction ---------------------------------------------- #

    def cdna_to_tx(self, p: CDNAPosition) -> int:
        if p.utr3:
            n = self.n_cds_end + p.base
        elif p.base < 0:
            n = self.n_cds_start + p.base          # base is negative
        elif p.base == 0:
            raise CDNAError("c.0 does not exist in HGVS numbering")
        else:
            n = self.n_cds_start + p.base - 1
        return n

    def to_genomic(self, position: Any) -> int:
        """
        Map a cDNA position (string, or a parsed :class:`CDNAPosition`) to a
        1-based genomic coordinate.  Raises :class:`CDNAError` if it falls
        outside the transcript.
        """
        p = position if isinstance(position, CDNAPosition) else parse_hgvs(position)
        n = self.cdna_to_tx(p)
        if n < 1:
            # e.g. c.-219 when the annotated 5'UTR is only 171 nt: a promoter
            # position, upstream of where this transcript is annotated to start
            utr5 = self.n_cds_start - 1
            raise CDNAError(
                f"{p} lies {1 - n} nt upstream of the start of "
                f"{self.tx.name}, which has only {utr5} nt of 5'UTR "
                "annotated. Promoter positions like this fall outside the "
                "transcript, so there is nowhere on the gene model to draw "
                "them. A transcript with a longer annotated 5' end would "
                "place it."
            )
        if n > self.tx_length:
            raise CDNAError(
                f"{p} is past the 3' end of {self.tx.name} "
                f"(the transcript is {self.tx_length} nt, and c.1 sits at "
                f"n.{self.n_cds_start})"
            )
        g = self.tx_to_genomic(n)
        if g is None:                                    # pragma: no cover
            raise CDNAError(f"{p} could not be placed on {self.tx.name}")
        if p.offset:
            # intronic: step off the exon edge, in transcript direction
            g += p.offset if self.strand == "+" else -p.offset
        return g

    # -- reporting --------------------------------------------------------- #

    def exon_of(self, position: Any) -> Optional[int]:
        """Biological exon number containing a cDNA position (None if intronic)."""
        p = position if isinstance(position, CDNAPosition) else parse_hgvs(position)
        if p.offset:
            return None
        g = self.to_genomic(p)
        for idx, (s, e) in enumerate(self.tx.exons):
            if s < g <= e:
                return self.tx.exon_number(idx)
        return None

    def protein_position(self, position: Any) -> Optional[int]:
        """Codon number for a coding cDNA position (1-based), else None."""
        p = position if isinstance(position, CDNAPosition) else parse_hgvs(position)
        if p.utr3 or p.base < 0 or p.offset:
            return None
        return (p.base - 1) // 3 + 1


# --------------------------------------------------------------------------- #
#  Variants
# --------------------------------------------------------------------------- #

@dataclass
class Variant:
    """One variant to place on the cartoon."""
    label: str = ""
    cdna: str = ""                     # what the user typed, e.g. "c.743G>A"
    genomic: Optional[int] = None      # 1-based; filled in by placement
    category: str = ""                 # free text; drives colour and legend
    count: int = 1                     # recurrence, scales the head
    color: Optional[str] = None        # explicit override
    exon: Optional[int] = None
    protein: Optional[int] = None
    note: str = ""
    error: str = ""                    # why it couldn't be placed

    @property
    def placed(self) -> bool:
        return self.genomic is not None and not self.error


def place_variants(variants: Sequence[Variant], transcript: Transcript,
                   ) -> Tuple[List[Variant], List[Variant]]:
    """
    Resolve each variant's cDNA position against a transcript.

    Returns ``(placed, failed)``.  Variants that already carry a genomic
    coordinate are passed through untouched, so a file can mix the two.
    """
    mapper = CDNAMapper(transcript)
    placed, failed = [], []
    for v in variants:
        if v.genomic is not None and not v.cdna:
            placed.append(v)
            continue
        try:
            v.genomic = mapper.to_genomic(v.cdna)
            v.exon = mapper.exon_of(v.cdna)
            v.protein = mapper.protein_position(v.cdna)
            v.error = ""
            placed.append(v)
        except (CDNAError, ValueError) as e:
            v.error = str(e)
            failed.append(v)
    return placed, failed


# --------------------------------------------------------------------------- #
#  Style
# --------------------------------------------------------------------------- #

@dataclass
class Style:
    """Every knob for the drawing.  Dump to JSON, edit, reuse."""

    # ---- canvas ----
    figure_width_in: float = 9.0
    row_height_in: float = 0.85          # vertical space per transcript
    margin_in: float = 0.45
    dpi: int = 400
    background: str = "none"             # "none" = transparent, or a colour
    font_family: str = "DejaVu Sans"
    font_size: float = 9.0

    # ---- glyph geometry (fraction of a row) ----
    cds_height: float = 0.42
    utr_height: float = 0.22
    noncoding_height: float = 0.26
    intron_linewidth: float = 1.2
    exon_edge_width: float = 0.8
    exon_shape: str = "box"              # box | round
    corner_radius: float = 0.30          # round: radius as a fraction of height

    # ---- colours ----
    cds_color: str = "#2C6FA6"
    cds_edge: str = "#123A5C"
    utr_color: str = "#BBD5E8"
    utr_edge: str = "#123A5C"
    noncoding_color: str = "#B7B7B7"
    noncoding_edge: str = "#5A5A5A"
    intron_color: str = "#4A4A4A"
    text_color: str = "#1A1A1A"
    highlight_color: str = "#D1495B"

    # ---- introns ----
    intron_mode: str = "compress"        # linear | compress | equal
    intron_width_frac: float = 0.06      # compressed intron width, rel. to mean exon
    intron_style: str = "line"           # line | chevron_line | angled
    angled_height: float = 0.30          # peak height for intron_style = angled
    chevron_spacing: float = 0.022       # axis fraction between strand arrows
    chevron_size: float = 0.006
    show_intron_breaks: bool = True      # "//" marks on compressed introns

    # ---- labels ----
    show_exon_numbers: bool = True
    exon_number_size: float = 6.5
    exon_number_min_width: float = 0.006  # 'inside' only: skip narrow boxes
    exon_number_position: str = "above"   # above | inside | below
    show_transcript_labels: bool = True
    transcript_label_size: float = 8.0
    transcript_label_side: str = "left"   # left | right
    show_title: bool = True
    title_size: float = 12.0
    title_weight: str = "bold"
    show_subtitle: bool = True            # locus / assembly / track line
    subtitle_size: float = 7.5
    show_strand_arrow: bool = True        # 5'->3' arrow under the gene
    show_legend: bool = True
    legend_size: float = 7.5

    # ---- axis ----
    axis_mode: str = "auto"               # auto | coords | scalebar | none
    scalebar_bp: int = 0                  # 0 = choose automatically
    axis_size: float = 7.0
    orient_five_prime_left: bool = True   # flip minus-strand genes

    # ---- annotations ----
    annotation_height: float = 0.18
    annotation_label_size: float = 7.0
    annotation_default_color: str = "#E8A33D"

    # ---- variants ----
    #: "stacked"  -- markers sit directly above their position and pile up
    #:              vertically, label alongside each marker.  No stems.
    #: "lanes"    -- one labelled row per category, so mutation types can be
    #:              compared along the gene instead of being intermixed.
    #: "lollipop" -- classic stems, with heads nudged sideways when crowded.
    variant_style: str = "stacked"
    lane_gap: float = 0.16                 # extra space between category lanes
    lane_label_size: float = 7.5
    lane_rule: bool = True                 # faint baseline under each lane
    lane_rule_color: str = "#E4E4E4"
    variant_marker: str = "o"              # o | s | D | v | ^
    variant_head_size: float = 42.0        # area in pt^2 of a single-count head
    variant_head_edge: str = "#FFFFFF"
    variant_head_edge_width: float = 0.5
    variant_base_gap: float = 0.20         # gap between gene and the first tier
    #: 0 = derive the row pitch from the marker size (recommended: shrinking
    #: the markers then genuinely shortens the figure).  Set > 0 to pin it.
    variant_stack_gap: float = 0.0
    variant_row_spacing: float = 1.30      # row pitch, in marker diameters
    variant_pack: float = 1.25             # horizontal clearance, in diameters
    variant_label_gap: float = 0.004       # marker-to-label gap (axis fraction)
    # lollipop mode only
    variant_stem_color: str = "#9A9A9A"
    variant_stem_width: float = 0.8
    variant_stem_height: float = 0.55      # rows between gene and the heads
    variant_collision: str = "spread"      # spread (nudge sideways) | stack
    variant_scale_by_count: bool = True    # bigger head for recurrent variants
    variant_max_scale: float = 3.2         # cap on head area multiplier
    show_variant_labels: bool = True
    variant_label_size: float = 6.5
    variant_label_rotation: float = 45.0   # degrees; 0 = horizontal
    variant_label_max: int = 40            # skip labels beyond this many variants
    show_variant_legend: bool = True
    variant_default_color: str = "#D1495B"
    variant_palette: List[str] = field(default_factory=lambda: [
        "#D1495B", "#2C6FA6", "#E8A33D", "#4A7C59", "#7A5195",
        "#00898A", "#BC5090", "#8C8C8C",
    ])
    variant_axis: bool = False             # draw a c. position axis under the gene

    @classmethod
    def load(cls, path: Optional[str]) -> "Style":
        s = cls()
        if path:
            data = json.loads(Path(path).read_text())
            known = {f.name for f in fields(cls)}
            unknown = set(data) - known
            if unknown:
                print(f"[style] ignoring unknown keys: {sorted(unknown)}",
                      file=sys.stderr)
            for k, v in data.items():
                if k in known:
                    setattr(s, k, v)
        return s

    def dump(self, path: str) -> None:
        Path(path).write_text(json.dumps(asdict(self), indent=2) + "\n")


# --------------------------------------------------------------------------- #
#  Genomic -> axis coordinate mapping
# --------------------------------------------------------------------------- #

class CoordinateMapper:
    """
    Piecewise-linear map from genomic bp to axis x in [0, 1].

    In 'linear' mode this is a single stretch, so the picture is to scale but
    a gene with a 100 kb intron shows exons as invisible slivers.  In
    'compress'/'equal' mode introns are squashed to a constant width -- the
    convention used in almost every published gene-structure figure.
    """

    def __init__(self, transcripts: Sequence[Transcript], style: Style,
                 flip: bool = False,
                 bounds: Optional[Tuple[int, int]] = None):
        self.style = style
        self.flip = flip
        if bounds:
            self.g_start, self.g_end = bounds
        else:
            self.g_start = min(t.tx_start for t in transcripts)
            self.g_end = max(t.tx_end for t in transcripts)
        self.bounded = bounds is not None
        self.mode = style.intron_mode
        self.blocks: List[Tuple[int, int, float, float, bool]] = []  # gs,ge,xs,xe,is_exon
        self._starts: List[int] = []
        self._build(transcripts)

    # -- construction ------------------------------------------------------ #

    def _merged_exons(self, transcripts: Sequence[Transcript]) -> List[Tuple[int, int]]:
        spans = sorted(e for t in transcripts for e in t.exons)
        merged: List[List[int]] = []
        for s, e in spans:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        # a focused view only lays out what falls inside the window
        out = []
        for s, e in merged:
            s, e = max(s, self.g_start), min(e, self.g_end)
            if s < e:
                out.append((s, e))
        return out

    def _build(self, transcripts: Sequence[Transcript]) -> None:
        total = self.g_end - self.g_start
        if self.mode == "linear" or total == 0:
            self.blocks = [(self.g_start, self.g_end, 0.0, 1.0, True)]
            self._starts = [self.g_start]
            self.bp_per_x = total if total else 1
            return

        exons = self._merged_exons(transcripts)
        if not exons:
            self.blocks = [(self.g_start, self.g_end, 0.0, 1.0, True)]
            self._starts = [self.g_start]
            self.bp_per_x = total or 1
            return

        # Interleave exon blocks and gap blocks across the full span.
        raw: List[Tuple[int, int, bool]] = []
        cursor = self.g_start
        for s, e in exons:
            if s > cursor:
                raw.append((cursor, s, False))
            raw.append((s, e, True))
            cursor = e
        if cursor < self.g_end:
            raw.append((cursor, self.g_end, False))

        exon_lens = [e - s for s, e, is_ex in raw if is_ex]
        mean_exon = sum(exon_lens) / len(exon_lens)
        gap_units = self.style.intron_width_frac * len(exon_lens) * mean_exon
        gap_units = max(gap_units, mean_exon * 0.15)

        # weight per block, in arbitrary units
        weights = []
        for s, e, is_ex in raw:
            if is_ex:
                weights.append(float(e - s) if self.mode == "compress"
                               else mean_exon)
            else:
                weights.append(gap_units)
        wsum = sum(weights) or 1.0

        x = 0.0
        for (s, e, is_ex), w in zip(raw, weights):
            xw = w / wsum
            self.blocks.append((s, e, x, x + xw, is_ex))
            x += xw
        self._starts = [b[0] for b in self.blocks]
        # bp per unit x, measured on exon blocks (the honest scale in compress mode)
        ex = [(b[1] - b[0]) / (b[3] - b[2]) for b in self.blocks
              if b[4] and b[3] > b[2]]
        self.bp_per_x = sum(ex) / len(ex) if ex else total

    # -- use --------------------------------------------------------------- #

    def x(self, pos: int) -> float:
        pos = max(self.g_start, min(self.g_end, pos))
        i = bisect.bisect_right(self._starts, pos) - 1
        i = max(0, min(i, len(self.blocks) - 1))
        gs, ge, xs, xe, _ = self.blocks[i]
        frac = 0.0 if ge == gs else (pos - gs) / (ge - gs)
        val = xs + frac * (xe - xs)
        return 1.0 - val if self.flip else val

    def span(self, start: int, end: int) -> Tuple[float, float]:
        a, b = self.x(start), self.x(end)
        return (a, b) if a <= b else (b, a)

    def visible(self, start: int, end: int) -> bool:
        """True if a feature overlaps the drawn window at all."""
        return end > self.g_start and start < self.g_end

    def clip(self, start: int, end: int) -> Optional[Tuple[int, int]]:
        """Trim a feature to the window, or None if it falls outside."""
        s, e = max(start, self.g_start), min(end, self.g_end)
        return (s, e) if s < e else None

    @property
    def gap_blocks(self) -> List[Tuple[int, int]]:
        return [(b[0], b[1]) for b in self.blocks if not b[4]]


# --------------------------------------------------------------------------- #
#  Renderer
# --------------------------------------------------------------------------- #

class GeneCartoon:
    def __init__(self, transcripts: Sequence[Transcript], style: Style,
                 gene: str, genome: str, track: str,
                 annotations: Optional[List[Dict[str, Any]]] = None,
                 links: Optional[Dict[str, str]] = None,
                 variants: Optional[Sequence[Variant]] = None,
                 focus: Optional[Tuple[int, int]] = None,
                 focus_label: str = ""):
        if not transcripts:
            raise ValueError("no transcripts to draw")
        self.tx = list(transcripts)
        self.focus = focus
        self.focus_label = focus_label
        self.st = style
        self.gene = gene
        self.genome = genome
        self.track = track
        self.annotations = annotations or []
        self.links = links or {}
        all_v = [v for v in (variants or []) if v.placed]
        if focus:
            g0, g1 = focus
            self.variants = [v for v in all_v if g0 < v.genomic <= g1]
            self.offscreen_variants = [v for v in all_v
                                       if not (g0 < v.genomic <= g1)]
        else:
            self.variants = all_v
            self.offscreen_variants = []
        self._variant_colors: Dict[str, str] = self._assign_variant_colors()
        self._lane_rows: List[Tuple[str, int, int]] = []
        self._gids: Dict[int, int] = {id(v): i
                                      for i, v in enumerate(self.variants)}
        self.strand = self.tx[0].strand
        self.flip = style.orient_five_prime_left and self.strand == "-"
        self.map = CoordinateMapper(self.tx, style, flip=self.flip,
                                    bounds=focus)

    # -- helpers ----------------------------------------------------------- #

    # x and y are in wildly different units (axis fraction vs. row height), so
    # anything that should look isotropic on paper has to be converted through
    # inches.  These two are set in render() once the axes limits are known.
    _x_per_in: float = 1.0     # axis x-units per inch
    _y_per_in: float = 1.0     # axis y-units per inch
    #: how far in from the axis edge left/right labels are drawn
    LABEL_DX: float = 0.015

    def _rect(self, ax, x0: float, x1: float, ycen: float, h: float,
              face: str, edge: str, url: Optional[str] = None, zorder: int = 3):
        w = max(x1 - x0, 1e-5)
        if self.st.exon_shape == "round" and self.st.corner_radius > 0:
            # keep the radius isotropic on paper, and never more than half of
            # the shorter side -- otherwise thin exons turn into ellipses
            w_in, h_in = w / self._x_per_in, h / self._y_per_in
            r_in = min(self.st.corner_radius * h_in, w_in / 2, h_in / 2)
            rx, ry = r_in * self._x_per_in, r_in * self._y_per_in
            p = Polygon(_rounded_rect_pts(x0, ycen - h / 2, w, h, rx, ry),
                        closed=True, linewidth=self.st.exon_edge_width,
                        facecolor=face, edgecolor=edge, zorder=zorder,
                        joinstyle="round")
        else:
            p = Rectangle((x0, ycen - h / 2), w, h,
                          linewidth=self.st.exon_edge_width,
                          facecolor=face, edgecolor=edge, zorder=zorder)
        if url:
            p.set_url(url)
        ax.add_patch(p)
        return p

    def _text_width(self, text: str, size: float) -> float:
        """Rough width of a string in axis x-units."""
        return _text_width_in(text, size) * self._x_per_in

    def _text_h(self, size_pt: float) -> float:
        """
        Height of one line of text, in row units.

        Text is sized in points but the layout works in row units, and the two
        are only related through row_height_in. Reserving fixed row units for
        text is what used to break the figure when the font size or row height
        was changed.
        """
        return (size_pt * 1.35 / 72.0) * self._y_per_in

    # -- variants ---------------------------------------------------------- #

    def _assign_variant_colors(self) -> Dict[str, str]:
        """
        One colour per category, in order of first appearance.

        An explicit ``Variant.color`` wins, so a colour chosen in the GUI (or
        set in the spreadsheet) also drives the legend -- otherwise the key
        and the markers would disagree.
        """
        pal = self.st.variant_palette or [self.st.variant_default_color]
        out: Dict[str, str] = {}
        for v in self.variants:
            cat = v.category or ""
            if cat not in out:
                out[cat] = v.color or pal[len(out) % len(pal)]
            elif v.color and out[cat] != v.color:
                out[cat] = v.color
        return out

    def _gid_of(self, v: Variant) -> int:
        return self._gids.get(id(v), 0)

    def _marker_diameter(self) -> float:
        """Widest marker, in points (scatter sizes are areas, hence sqrt)."""
        size = self.st.variant_head_size
        if self.st.variant_scale_by_count and self.variants:
            cmax = max(max(1, v.count) for v in self.variants)
            if cmax > 1:
                size *= self.st.variant_max_scale
        return math.sqrt(max(size, 1.0))

    def _stack_gap(self) -> float:
        """
        Vertical pitch between stacked rows, in row units.

        Tied to the marker diameter by default, so turning the markers down
        actually shortens the figure instead of leaving holes between them.
        """
        st = self.st
        if st.variant_stack_gap > 0:
            return st.variant_stack_gap
        d_in = self._marker_diameter() / 72.0
        gap = d_in * self._y_per_in * st.variant_row_spacing
        if st.show_variant_labels and len(self.variants) <= st.variant_label_max:
            # rows also have to clear the labels sitting beside each marker,
            # or small markers stack their labels on top of one another
            gap = max(gap, self._text_h(st.variant_label_size) * 1.30)
        return gap

    def variant_tooltips(self) -> List[Dict[str, Any]]:
        """
        One record per drawn variant, keyed by the SVG element id.

        The app uses this to attach hover text and click-to-pin behaviour to
        the very same SVG it offers for download -- no second renderer, so
        what you interact with is exactly what you publish.
        """
        chrom = self.tx[0].chrom
        out = []
        for i, v in enumerate(self.variants):
            out.append({
                "gid": f"variant-{i}",
                "label": str(v.label or v.cdna or "variant"),
                "cdna": str(v.cdna or ""),
                "category": str(v.category or ""),
                "count": int(v.count or 1),
                "position": f"{chrom}:{v.genomic:,}" if v.genomic else "",
                "exon": "" if v.exon is None else f"exon {v.exon}",
                "codon": "" if v.protein is None else f"codon {v.protein}",
                "note": str(v.note or ""),
                "color": self.variant_color(v),
            })
        return out

    def variant_color(self, v: Variant) -> str:
        if v.color:
            return v.color
        if not v.category:
            return self.st.variant_default_color
        return self._variant_colors.get(v.category, self.st.variant_default_color)

    def _variant_layout(self) -> List[Tuple[Variant, float, float, int]]:
        """
        Position each variant's head.

        Returns ``(variant, anchor_x, head_x, tier)``.  ``anchor_x`` is the
        true genomic position -- where the stem meets the gene -- while
        ``head_x`` is where the head is drawn.  In "spread" mode (the default)
        colliding heads are nudged sideways and joined to their anchor by a
        bent stem, which keeps every label legible; in "stack" mode they pile
        up vertically instead.
        """
        st = self.st
        items = sorted(((v, self.map.x(v.genomic - 1)) for v in self.variants),
                       key=lambda p: p[1])
        if not items:
            return []

        head_w = (self._marker_diameter() / 72.0) * self._x_per_in * st.variant_pack
        if st.variant_style == "lollipop" and st.show_variant_labels \
                and len(items) <= st.variant_label_max:
            # lollipop labels share one baseline, so neighbouring heads must be
            # far enough apart for the rotated text not to run into each other
            rot = math.radians(st.variant_label_rotation or 90.0)
            need = (self._text_h(st.variant_label_size) / max(self._y_per_in, 1e-9)
                    ) * self._x_per_in / max(math.sin(rot), 0.2)
            head_w = max(head_w, need * 1.30)

        show_lbl = (st.show_variant_labels
                    and len(items) <= st.variant_label_max)

        def pack(members: Sequence[Tuple[Variant, float]]
                 ) -> List[Tuple[Variant, float, int]]:
            """Greedy row packing: a row is free once the previous marker
            (and its label) has ended."""
            ends: List[float] = []
            out = []
            for v, x in members:
                right = head_w / 2
                if show_lbl and v.label:
                    right += st.variant_label_gap + self._text_width(
                        v.label, st.variant_label_size) + head_w * 0.25
                lo = x - head_w / 2
                tier = 0
                while tier < len(ends) and ends[tier] > lo:
                    tier += 1
                if tier == len(ends):
                    ends.append(x + right)
                else:
                    ends[tier] = x + right
                out.append((v, x, tier))
            return out

        if st.variant_style == "stacked":
            # Markers stay on their true position and pile upwards.
            return [(v, x, x, tier) for v, x, tier in pack(items)]

        if st.variant_style == "lanes":
            # One row per category, in the same order as the legend, with the
            # busiest lane nearest the gene.  Lanes pack internally, so a
            # crowded type gets the height it needs without pushing the others.
            order = [c for c in self._variant_colors]
            self._lane_rows = []
            out, base = [], 0
            for cat in order:
                members = [(v, x) for v, x in items if (v.category or "") == cat]
                if not members:
                    continue
                packed = pack(members)
                n_sub = max(tier for _, _, tier in packed) + 1
                for v, x, tier in packed:
                    out.append((v, x, x, base + tier))
                self._lane_rows.append((cat, base, n_sub))
                base += n_sub
            return out

        if st.variant_collision == "stack":
            out, occupied = [], []
            for v, x in items:
                tier = 0
                while tier < len(occupied) and occupied[tier] > x - head_w:
                    tier += 1
                if tier == len(occupied):
                    occupied.append(x)
                else:
                    occupied[tier] = x
                out.append((v, x, x, tier))
            return out

        # -- spread: push right, then re-centre each packed cluster -- #
        xs = [x for _, x in items]
        hx = []
        for x in xs:
            hx.append(x if not hx else max(x, hx[-1] + head_w))

        i = 0
        while i < len(hx):
            j = i
            while j + 1 < len(hx) and hx[j + 1] - hx[j] <= head_w + 1e-9:
                j += 1
            if j > i:                                   # a real cluster
                want = sum(xs[i:j + 1]) / (j - i + 1)   # centre on the true mean
                have = (hx[i] + hx[j]) / 2
                shift = want - have
                lo_gap = hx[i] - (hx[i - 1] + head_w) if i > 0 else float("inf")
                shift = max(shift, -lo_gap)          # don't collide leftwards
                for k in range(i, j + 1):
                    hx[k] += shift
            i = j + 1

        return [(v, x, h, 0) for (v, x), h in zip(items, hx)]

    def _variant_extent(self, y_base: float) -> float:
        """Top y the lollipop track will reach -- used to reserve space."""
        return self._draw_variants(None, y_base)

    def _draw_variants(self, ax, y_base: float) -> float:
        """
        Draw the lollipop track above ``y_base`` and return the top y used.
        Pass ``ax=None`` to measure without drawing.
        """
        st = self.st
        layout = self._variant_layout()
        if not layout:
            return y_base

        gap = self._stack_gap()
        cmax = max(max(1, v.count) for v, _, _, _ in layout)
        label_this = (st.show_variant_labels
                      and len(layout) <= st.variant_label_max)

        def head_area(v: Variant) -> float:
            size = st.variant_head_size
            if st.variant_scale_by_count and cmax > 1:
                frac = (max(1, v.count) - 1) / (cmax - 1)
                size *= 1.0 + frac * (st.variant_max_scale - 1.0)
            return size

        # ---- no stems: markers directly above their position ---- #
        if st.variant_style in ("stacked", "lanes"):
            lanes = st.variant_style == "lanes"
            # in lanes mode each lane is nudged up by the lanes below it
            lane_of: Dict[int, int] = {}
            if lanes:
                for i, (_, base, n_sub) in enumerate(self._lane_rows):
                    for k in range(n_sub):
                        lane_of[base + k] = i

            def y_of(tier: int) -> float:
                extra = lane_of.get(tier, 0) * st.lane_gap if lanes else 0.0
                return (y_base + st.variant_base_gap
                        + tier * gap + extra)

            top = y_base
            for i, (v, x, _, tier) in enumerate(layout):
                y = y_of(tier)
                top = max(top, y)
                if ax is None:
                    continue
                sc = ax.scatter([x], [y], s=head_area(v),
                                marker=st.variant_marker,
                                facecolor=self.variant_color(v),
                                edgecolor=st.variant_head_edge,
                                linewidth=st.variant_head_edge_width, zorder=6,
                                clip_on=False)
                # id survives into the SVG, which is what makes the figure
                # hoverable/clickable in the app without a second renderer
                sc.set_gid(f"variant-{self._gid_of(v)}")
                if label_this and v.label:
                    ax.text(x + st.variant_label_gap
                            + (math.sqrt(head_area(v)) / 72.0
                               * self._x_per_in) / 2,
                            y, str(v.label), ha="left", va="center",
                            fontsize=st.variant_label_size,
                            color=st.text_color, zorder=7)

            if lanes and ax is not None:
                for cat, base, n_sub in self._lane_rows:
                    y_lo, y_hi = y_of(base), y_of(base + n_sub - 1)
                    if st.lane_rule:
                        ax.add_line(Line2D(
                            [0.0, 1.0],
                            [y_lo - gap * 0.55] * 2,
                            color=st.lane_rule_color, linewidth=0.6, zorder=1))
                    ax.text(-self.LABEL_DX, (y_lo + y_hi) / 2, cat,
                            ha="right", va="center",
                            fontsize=st.lane_label_size,
                            color=self._variant_colors.get(cat, st.text_color),
                            fontweight="bold", zorder=7)
            return top + gap * 0.6

        # every head sits on one baseline (or its tier, in stack mode)
        y_head = y_base + st.variant_stem_height
        top = max(y_head + tier * gap
                  for _, _, _, tier in layout)

        if ax is not None:
            for v, anchor, hx, tier in layout:
                y = y_head + tier * gap
                # bend the stem where the head had to be nudged sideways
                knee = y_base + st.variant_stem_height * 0.45
                ax.add_line(Line2D([anchor, anchor, hx, hx],
                                   [y_base, knee, y - 0.04, y],
                                   color=st.variant_stem_color,
                                   linewidth=st.variant_stem_width, zorder=4,
                                   solid_capstyle="round",
                                   solid_joinstyle="round"))
                ax.scatter([hx], [y], s=head_area(v), marker=st.variant_marker,
                           facecolor=self.variant_color(v),
                           edgecolor=st.variant_head_edge,
                           linewidth=st.variant_head_edge_width, zorder=6,
                           clip_on=False)

        if not (label_this and any(v.label for v, _, _, _ in layout)):
            return top

        # labels share a baseline above the tallest head, so rotated text runs
        # parallel and never crosses a neighbour
        y_label = top + 0.13
        if ax is not None:
            rot = st.variant_label_rotation
            for v, _, hx, _ in layout:
                if v.label:
                    ax.text(hx, y_label, str(v.label),
                            ha="left" if rot else "center", va="bottom",
                            rotation=rot, rotation_mode="anchor",
                            fontsize=st.variant_label_size,
                            color=st.text_color, zorder=7)

        longest = max((len(str(v.label)) for v, _, _, _ in layout if v.label),
                      default=0)
        rise = math.sin(math.radians(st.variant_label_rotation)) \
            if st.variant_label_rotation else 1.0
        text_h = longest * st.variant_label_size * 0.62 / 72.0 * self._y_per_in
        return y_label + (text_h * rise if st.variant_label_rotation
                          else st.variant_label_size / 72.0 * self._y_per_in)

    #: vertical step between wrapped legend rows, in row units
    VLEGEND_LINE = 0.22

    def _variant_legend_layout(self) -> Tuple[List[Tuple[str, float, int]], int]:
        """Place category keys, wrapping onto extra rows when they overrun."""
        st = self.st
        cats = [c for c in self._variant_colors if c]
        if not cats:
            return [], 0
        gap = self._text_width("nn", st.legend_size)
        items: List[Tuple[str, float, int]] = []
        x, line = 0.0, 0
        for cat in cats:
            w = gap * 0.5 + self._text_width(cat, st.legend_size) + gap
            if x > 0.0 and x + w > 1.0:
                line += 1
                x = 0.0
            items.append((cat, x, line))
            x += w
        return items, line + 1

    def _draw_variant_legend(self, ax, y: float) -> None:
        st = self.st
        items, _ = self._variant_legend_layout()
        gap = self._text_width("nn", st.legend_size)
        for cat, x, line in items:
            yy = y - line * self.VLEGEND_LINE
            ax.scatter([x], [yy], s=st.variant_head_size,
                       marker=st.variant_marker,
                       facecolor=self._variant_colors[cat],
                       edgecolor=st.variant_head_edge,
                       linewidth=st.variant_head_edge_width, clip_on=False)
            ax.text(x + gap * 0.5, yy, cat, ha="left", va="center",
                    fontsize=st.legend_size, color="#444444")

    def _chevron(self, ax, x: float, ycen: float, pointing_right: bool):
        s = self.st.chevron_size
        h = self.st.cds_height * 0.30
        d = 1 if pointing_right else -1
        ax.add_line(Line2D([x - d * s, x, x - d * s],
                           [ycen + h, ycen, ycen - h],
                           color=self.st.intron_color,
                           linewidth=self.st.intron_linewidth * 0.9,
                           solid_capstyle="round", zorder=2))

    def _draw_intron(self, ax, x0: float, x1: float, y: float, right: bool):
        st = self.st
        if x1 - x0 < 1e-6:
            return
        if st.intron_style == "angled":
            xm = (x0 + x1) / 2
            ax.add_line(Line2D([x0, xm, x1],
                               [y, y + st.angled_height * st.cds_height / 0.42 * 0.5, y],
                               color=st.intron_color, linewidth=st.intron_linewidth,
                               solid_joinstyle="miter", zorder=2))
        else:
            ax.add_line(Line2D([x0, x1], [y, y], color=st.intron_color,
                               linewidth=st.intron_linewidth, zorder=2))
        if st.intron_style == "chevron_line" and x1 - x0 > st.chevron_spacing * 0.6:
            n = max(1, int((x1 - x0) / st.chevron_spacing))
            for k in range(n):
                xc = x0 + (k + 0.5) * (x1 - x0) / n
                self._chevron(ax, xc, y, right)

    def _draw_break_marks(self, ax, t: Transcript, y: float):
        """Little '//' marks showing where introns were compressed."""
        st = self.st
        if (not st.show_intron_breaks or st.intron_mode == "linear"
                or st.intron_style == "angled"):
            return
        for gs, ge in self.map.gap_blocks:
            # only mark gaps that fall inside *this* transcript
            if ge <= t.tx_start or gs >= t.tx_end or ge - gs < 200:
                continue
            x0, x1 = self.map.span(gs, ge)
            xm = (x0 + x1) / 2
            h = self.st.cds_height * 0.34
            for off in (-0.0035, 0.0035):
                ax.add_line(Line2D([xm + off - 0.0022, xm + off + 0.0022],
                                   [y - h, y + h],
                                   color=self.st.intron_color,
                                   linewidth=self.st.intron_linewidth * 0.8,
                                   zorder=4))

    # -- main draw --------------------------------------------------------- #

    def figure(self):
        """Build and return the matplotlib Figure without writing anything."""
        st = self.st
        plt.rcParams.update({
            "font.family": st.font_family,
            "font.size": st.font_size,
            "svg.fonttype": "none",     # keep text editable in Illustrator
            "pdf.fonttype": 42,
            "text.color": st.text_color,
        })

        n = len(self.tx)
        half = max(st.cds_height, st.utr_height, st.noncoding_height) / 2

        self._y_per_in = 1.0 / st.row_height_in

        # ---- left margin: must fit *every* label written down the left ---- #
        left_labels: List[Tuple[str, float]] = []
        if st.show_transcript_labels and st.transcript_label_side == "left":
            left_labels += [(t.name, st.transcript_label_size) for t in self.tx]
        if st.variant_style == "lanes" and self.variants:
            left_labels += [(c, st.lane_label_size)
                            for c in self._variant_colors if c]
        need_in = max((_text_width_in(s, pt) for s, pt in left_labels),
                      default=0.0)

        right_in = 0.06
        if st.show_transcript_labels and st.transcript_label_side == "right":
            right_in = max(_text_width_in(t.name, st.transcript_label_size)
                           for t in self.tx) + 0.06

        # x spans [0,1] for the gene plus whatever the margins need. Solved
        # iteratively because the axis span and the margin depend on each other.
        x_lo, x_hi = -0.03, 1.03
        for _ in range(5):
            span = x_hi - x_lo
            # labels start LABEL_DX in from the axis edge, so the margin has to
            # cover the text *plus* that offset
            x_lo = -((need_in / st.figure_width_in) * span
                     + (self.LABEL_DX + 0.012 if need_in else 0.03))
            x_hi = 1.0 + ((right_in / st.figure_width_in) * span
                          + (self.LABEL_DX if right_in > 0.06 else 0.0))
        self._x_per_in = (x_hi - x_lo) / st.figure_width_in

        # ---- headings: shrink to fit rather than run off the page ---- #
        t0 = self.tx[0]
        if self.focus:
            _g0, _g1 = self.focus
        else:
            _g0 = min(t.tx_start for t in self.tx)
            _g1 = max(t.tx_end for t in self.tx)
        subtitle = (f"{self.genome}  {t0.chrom}:{_g0 + 1:,}–{_g1:,}  "
                    f"({self.strand} strand, {_g1 - _g0:,} bp)  ·  {self.track}")
        if self.focus_label:
            subtitle = f"{self.focus_label}  ·  " + subtitle

        def fit(text: str, size: float) -> float:
            avail = st.figure_width_in * 0.90
            need = _text_width_in(text, size)
            return size * avail / need if need > avail else size

        title_size = fit(self.gene, st.title_size)
        subtitle_size = fit(subtitle, st.subtitle_size)

        # ---------- vertical layout, all in "row units" ---------- #
        # Rows run from y = n-1 (top) down to y = 0 (bottom).
        num_above = st.show_exon_numbers and st.exon_number_position == "above"
        num_below = st.show_exon_numbers and st.exon_number_position == "below"
        # exon numbers may stagger onto a second tier, hence two text heights
        num_band = 0.10 + 2 * self._text_h(st.exon_number_size) + 0.17
        y = (n - 1) + half + (num_band if num_above else 0.10)

        # variant lollipops sit directly above the top transcript
        y_variants = y if self.variants else None
        if self.variants:
            y = self._variant_extent(y)

        ann_y: Dict[int, float] = {}
        if self.annotations:
            rows = sorted({int(a.get("row", 0)) for a in self.annotations})
            top_tier = 0
            for i, r in enumerate(rows):
                members = [a for a in self.annotations if int(a.get("row", 0)) == r]
                tiers = self._annotation_tiers(members)
                span = 0.46 + 0.20 * (max(tiers.values()) if tiers else 0)
                ann_y[r] = y + 0.34 + sum(
                    0.46 for _ in rows[:i]) + 0.20 * top_tier
                top_tier += max(tiers.values()) if tiers else 0
                if i == len(rows) - 1:
                    y = ann_y[r] + span - 0.12

        y_sub = y + 0.16
        y_title = y_sub + (self._text_h(subtitle_size) + 0.06
                           if st.show_subtitle else 0.0)
        y_max = y_title + (self._text_h(title_size) + 0.08
                           if st.show_title else 0.05)

        b = -half - ((0.10 + self._text_h(st.exon_number_size) + 0.17)
                     if num_below else 0.06)                  # bottom of models
        if st.show_strand_arrow:
            y_arrow = b - 0.16 - (self._text_h(st.axis_size)
                                  if num_below else 0.0)
            b = y_arrow - self._text_h(st.axis_size) - 0.06
        else:
            y_arrow = None
        if st.axis_mode != "none":
            y_axis = b - 0.16
            # coords mode prints tick labels *and* the chromosome name
            lines = 2 if self._axis_kind() == "coords" else 1
            b = y_axis - 0.12 - lines * (self._text_h(st.axis_size) + 0.06)
        else:
            y_axis = None
        if st.show_legend:
            y_legend = b - 0.10 - self._text_h(st.legend_size) / 2
            b = y_legend - self._text_h(st.legend_size) / 2 - 0.08
        else:
            y_legend = None
        # in lanes mode the row labels already name every category
        if (self.variants and st.show_variant_legend
                and st.variant_style != "lanes"
                and any(self._variant_colors)):
            y_vlegend = b - 0.24
            _, vlines = self._variant_legend_layout()
            b = y_vlegend - 0.20 - (vlines - 1) * self.VLEGEND_LINE
        else:
            y_vlegend = None
        y_min = b - 0.08

        height = (y_max - y_min) * st.row_height_in
        fig = plt.figure(figsize=(st.figure_width_in, height))
        if st.background != "none":
            fig.patch.set_facecolor(st.background)
        else:
            fig.patch.set_alpha(0.0)

        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_axis_off()
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_min, y_max)

        # ---- transcripts ---- #
        for i, t in enumerate(self.tx):
            self._draw_transcript(ax, t, n - 1 - i)

        # ---- variants ---- #
        if y_variants is not None:
            self._draw_variants(ax, y_variants)

        # ---- annotations ---- #
        for r, yr in ann_y.items():
            self._draw_annotations(ax, yr,
                                   [a for a in self.annotations
                                    if int(a.get("row", 0)) == r])

        # ---- title / subtitle ---- #
        if st.show_title:
            ax.text(0.5, y_title, self.gene, ha="center", va="bottom",
                    fontsize=title_size, fontweight=st.title_weight,
                    color=st.text_color, url=self.links.get("browser"))
        if st.show_subtitle:
            ax.text(0.5, y_sub, subtitle, ha="center", va="bottom",
                    fontsize=subtitle_size, color="#666666",
                    url=self.links.get("sequence"))

        if y_arrow is not None:
            self._draw_strand_arrow(ax, y_arrow)
        if y_axis is not None:
            self._draw_axis(ax, y_axis)
        if y_legend is not None:
            self._draw_legend(ax, y_legend)
        if y_vlegend is not None:
            self._draw_variant_legend(ax, y_vlegend)

        return fig

    def _save_kw(self) -> Dict[str, Any]:
        return dict(dpi=self.st.dpi,
                    transparent=(self.st.background == "none"),
                    bbox_inches="tight", pad_inches=self.st.margin_in / 2)

    def render(self, out_paths: Sequence[str], quiet: bool = False) -> List[str]:
        """Write the cartoon to one or more files; format follows the extension."""
        fig = self.figure()
        for p in out_paths:
            Path(p).parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(p, **self._save_kw())
            if not quiet:
                print(f"wrote {p}")
        plt.close(fig)
        return list(out_paths)

    def to_svg(self) -> str:
        """The cartoon as an SVG string -- handy for inline display."""
        fig = self.figure()
        buf = io.StringIO()
        fig.savefig(buf, format="svg", **self._save_kw())
        plt.close(fig)
        return buf.getvalue()

    def to_bytes(self, fmt: str = "png") -> bytes:
        """The cartoon as raw bytes in any matplotlib format (png/pdf/svg/tif)."""
        fig = self.figure()
        buf = io.BytesIO()
        fig.savefig(buf, format=fmt, **self._save_kw())
        plt.close(fig)
        return buf.getvalue()

    def _draw_transcript(self, ax, t: Transcript, y: float) -> None:
        st = self.st
        seq_url = self.links.get("sequence")

        # introns
        pointing_right = (t.strand == "+") != self.flip
        for gs, ge in t.introns:
            piece = self.map.clip(gs, ge)
            if piece:
                self._draw_intron(ax, *self.map.span(*piece), y,
                                  pointing_right)
        if not t.introns:
            piece = self.map.clip(t.tx_start, t.tx_end)
            if piece:
                self._draw_intron(ax, *self.map.span(*piece), y,
                                  pointing_right)
        self._draw_break_marks(ax, t, y)

        # exon segments
        for gs, ge, kind in t.segments():
            piece = self.map.clip(gs, ge)
            if not piece:
                continue
            x0, x1 = self.map.span(*piece)
            if kind == "cds":
                h, fc, ec = st.cds_height, st.cds_color, st.cds_edge
            elif kind == "nc":
                h, fc, ec = st.noncoding_height, st.noncoding_color, st.noncoding_edge
            else:
                h, fc, ec = st.utr_height, st.utr_color, st.utr_edge
            self._rect(ax, x0, x1, y, h, fc, ec, url=seq_url)

        # exon numbers
        if st.show_exon_numbers:
            inside = st.exon_number_position == "inside"
            # Pack labels into as few tiers as will hold them, left to right.
            # NB: iterate in *drawn* order, not genomic order. A minus-strand
            # gene is flipped, so assuming ascending x mis-tiers every label.
            entries = []
            for idx, (gs, ge) in enumerate(t.exons):
                piece = self.map.clip(gs, ge)
                if not piece:
                    continue
                x0, x1 = self.map.span(*piece)
                entries.append(((x0 + x1) / 2, t.exon_number(idx), x1 - x0))
            entries.sort(key=lambda e: e[0])

            ends: List[float] = []          # rightmost x used, per tier
            for xm, num, box_w in entries:
                if inside:
                    # a number can only sit inside a box wide enough to hold it
                    if box_w < st.exon_number_min_width:
                        continue
                    ax.text(xm, y, str(num), ha="center", va="center",
                            fontsize=st.exon_number_size, color="white",
                            fontweight="bold", zorder=5)
                    continue

                half_w = self._text_width(str(num), st.exon_number_size) / 2
                lo, hi = xm - half_w - 0.003, xm + half_w + 0.003
                tier = 0
                while tier < len(ends) and ends[tier] > lo:
                    tier += 1
                if tier >= 2:               # two tiers is the readable limit;
                    continue                # skipping beats overprinting
                if tier == len(ends):
                    ends.append(hi)
                else:
                    ends[tier] = hi

                sign = -1 if st.exon_number_position == "below" else 1
                dy = sign * (st.cds_height / 2 + 0.09 + tier * 0.17)
                if tier:  # leader line from the label down to its exon
                    ax.add_line(Line2D(
                        [xm, xm],
                        [y + sign * (st.cds_height / 2 + 0.02),
                         y + dy - sign * 0.015],
                        color="#AAAAAA", linewidth=0.45, zorder=1))
                ax.text(xm, y + dy, str(num), ha="center",
                        va="bottom" if sign > 0 else "top",
                        fontsize=st.exon_number_size, color=st.text_color,
                        zorder=5)

        # transcript label
        if st.show_transcript_labels:
            label = t.name if len(self.tx) > 1 or t.name != self.gene else t.name
            if st.transcript_label_side == "left":
                ax.text(-self.LABEL_DX, y, label, ha="right", va="center",
                        fontsize=st.transcript_label_size, color=st.text_color)
            else:
                ax.text(1.0 + self.LABEL_DX, y, label, ha="left", va="center",
                        fontsize=st.transcript_label_size, color=st.text_color)

    def _draw_strand_arrow(self, ax, y: float) -> None:
        st = self.st
        x0, x1 = (0.0, 1.0)
        ax.annotate("", xy=(x1, y), xytext=(x0, y),
                    arrowprops=dict(arrowstyle="-|>", color="#8A8A8A",
                                    linewidth=0.9, shrinkA=0, shrinkB=0))
        ax.text(x0, y + 0.055, "5′", ha="left", va="bottom",
                fontsize=st.axis_size, color="#8A8A8A")
        ax.text(x1, y + 0.055, "3′", ha="right", va="bottom",
                fontsize=st.axis_size, color="#8A8A8A")

    def _annotation_tiers(self, items: Sequence[Dict[str, Any]]) -> Dict[int, int]:
        """Assign each labelled annotation a stacking tier so labels don't collide."""
        st = self.st
        out: Dict[int, int] = {}
        placed: List[Tuple[float, float, int]] = []
        for a in sorted(items, key=lambda d: int(d["start"])):
            label = a.get("label")
            if not label:
                continue
            start, end = int(a["start"]), int(a.get("end", a["start"]))
            if a.get("coords", "genomic") == "browser":
                start -= 1
            x0, x1 = self.map.span(start, max(end, start + 1))
            w = self._text_width(label, st.annotation_label_size) + 0.01
            xm = (x0 + x1) / 2
            lo, hi = xm - w / 2, xm + w / 2
            tier = 0
            while any(pt == tier and not (hi < p0 or lo > p1)
                      for p0, p1, pt in placed):
                tier += 1
            placed.append((lo, hi, tier))
            out[id(a)] = tier
        return out

    def _draw_annotations(self, ax, y: float,
                          items: Sequence[Dict[str, Any]]) -> None:
        st = self.st
        tiers = self._annotation_tiers(items)
        for a in sorted(items, key=lambda d: int(d["start"])):
            color = a.get("color", st.annotation_default_color)
            start, end = int(a["start"]), int(a.get("end", a["start"]))
            if a.get("coords", "genomic") == "browser":  # 1-based inclusive
                start -= 1
            end = max(end, start + 1)
            if not self.map.visible(start, end):
                continue
            x0, x1 = self.map.span(*self.map.clip(start, end))
            kind = a.get("style", "box")
            if kind == "marker":
                xm = (x0 + x1) / 2
                ax.add_patch(Polygon(
                    [(xm, y - 0.10), (xm - 0.007, y + 0.10), (xm + 0.007, y + 0.10)],
                    closed=True, facecolor=color, edgecolor="none", zorder=6))
            elif kind == "bracket":
                h = st.annotation_height / 2
                ax.add_line(Line2D([x0, x0, x1, x1], [y - h, y + h, y + h, y - h],
                                   color=color, linewidth=1.2, zorder=6))
            else:
                self._rect(ax, x0, x1, y, st.annotation_height, color,
                           a.get("edge", color), zorder=6)
            label = a.get("label")
            if label:
                xm = (x0 + x1) / 2
                tier = tiers.get(id(a), 0)
                dy = st.annotation_height / 2 + 0.06 + tier * 0.20
                if tier:
                    ax.add_line(Line2D([xm, xm],
                                       [y + st.annotation_height / 2, y + dy],
                                       color="#AAAAAA", linewidth=0.45, zorder=1))
                ax.text(xm, y + dy, label, ha="center", va="bottom",
                        fontsize=st.annotation_label_size, color=st.text_color)

    def _axis_kind(self) -> str:
        m = self.st.axis_mode
        if m == "auto":
            return "coords" if self.st.intron_mode == "linear" else "scalebar"
        return m

    def _draw_axis(self, ax, y: float) -> None:
        st = self.st
        mode = self._axis_kind()
        if mode == "none":
            return

        if mode == "coords":
            g0, g1 = self.map.g_start, self.map.g_end
            ticks = _nice_ticks(g0, g1, 5)
            ax.add_line(Line2D([0, 1], [y, y], color="#555555", linewidth=0.8))
            # In compressed mode the mapping is piecewise, so evenly spaced
            # coordinates land unevenly on the axis. Draw every tick mark, but
            # only label the ones that have room -- otherwise they overprint.
            placed: List[Tuple[float, float]] = []
            for tpos in ticks:
                xt = self.map.x(tpos)
                ax.add_line(Line2D([xt, xt], [y, y - 0.07],
                                   color="#555555", linewidth=0.8))
                label = f"{tpos:,}"
                half_w = self._text_width(label, st.axis_size) / 2 + 0.004
                if any(not (xt + half_w < lo or xt - half_w > hi)
                       for lo, hi in placed):
                    continue
                placed.append((xt - half_w, xt + half_w))
                ax.text(xt, y - 0.11, label, ha="center", va="top",
                        fontsize=st.axis_size, color="#555555")
            ax.text(0.5, y - 0.13 - self._text_h(st.axis_size),
                    self.tx[0].chrom, ha="center", va="top",
                    fontsize=st.axis_size, color="#555555")
        else:
            bp = st.scalebar_bp or _nice_scalebar(self.map.bp_per_x)
            width = bp / self.map.bp_per_x
            if width > 0.45:                      # too wide -> drop a decade
                bp = max(1, bp // 10)
                width = bp / self.map.bp_per_x
            x0 = 0.0
            ax.add_line(Line2D([x0, x0 + width], [y, y],
                               color="#333333", linewidth=1.6,
                               solid_capstyle="butt"))
            for xe in (x0, x0 + width):
                ax.add_line(Line2D([xe, xe], [y - 0.05, y + 0.05],
                                   color="#333333", linewidth=1.6))
            ax.text(x0 + width / 2, y - 0.09, _fmt_bp(bp), ha="center",
                    va="top", fontsize=st.axis_size, color="#333333")
            if st.intron_mode != "linear":
                ax.text(x0 + width + 0.02, y, "(exon scale; introns not to scale)",
                        ha="left", va="center", fontsize=st.axis_size * 0.9,
                        color="#999999", style="italic")

    def _draw_legend(self, ax, y: float) -> None:
        st = self.st
        items = [("CDS", st.cds_color, st.cds_edge, st.cds_height * 0.75),
                 ("UTR", st.utr_color, st.utr_edge, st.utr_height * 0.9)]
        if any(not t.coding for t in self.tx):
            items.append(("non-coding exon", st.noncoding_color,
                          st.noncoding_edge, st.noncoding_height * 0.9))
        items.append(("intron", None, st.intron_color, 0))

        swatch = 0.026
        gap = self._text_width("nn", st.legend_size)
        x = 0.0
        for label, fc, ec, h in items:
            if fc is None:
                ax.add_line(Line2D([x, x + swatch], [y, y], color=ec,
                                   linewidth=st.intron_linewidth))
            else:
                self._rect(ax, x, x + swatch, y, h, fc, ec)
            ax.text(x + swatch + gap * 0.45, y, label, ha="left", va="center",
                    fontsize=st.legend_size, color="#444444")
            x += (swatch + gap * 0.45
                  + self._text_width(label, st.legend_size) + gap)


# --------------------------------------------------------------------------- #
#  small utilities
# --------------------------------------------------------------------------- #

def _nice_ticks(lo: int, hi: int, target: int) -> List[int]:
    span = max(hi - lo, 1)
    raw = span / max(target, 1)
    mag = 10 ** int(f"{raw:e}".split("e")[1])
    for m in (1, 2, 2.5, 5, 10):
        step = m * mag
        if span / step <= target * 1.4:
            break
    step = int(step) or 1
    first = ((lo // step) + 1) * step
    return [int(t) for t in range(first, hi, step)]


def _nice_scalebar(bp_per_x: float) -> int:
    target = bp_per_x * 0.20            # aim for ~20% of the axis
    mag = 10 ** int(f"{max(target, 1):e}".split("e")[1])
    for m in (1, 2, 5, 10):
        if m * mag >= target:
            return int(m * mag)
    return int(10 * mag)


def _text_width_in(text: Any, size_pt: float) -> float:
    """Rough width of a string in inches at a given point size."""
    return len(str(text)) * size_pt * 0.60 / 72.0


def _rounded_rect_pts(x: float, y: float, w: float, h: float,
                      rx: float, ry: float, n: int = 6) -> List[Tuple[float, float]]:
    """Vertices of a rectangle with independently scaled corner radii."""
    import math
    rx, ry = max(rx, 0.0), max(ry, 0.0)
    if rx <= 0 or ry <= 0:
        return [(x, y), (x + w, y), (x + w, y + h), (x, y + h)]
    corners = [  # centre, start angle
        ((x + w - rx, y + ry), 270.0),
        ((x + w - rx, y + h - ry), 0.0),
        ((x + rx, y + h - ry), 90.0),
        ((x + rx, y + ry), 180.0),
    ]
    pts: List[Tuple[float, float]] = []
    for (cx, cy), a0 in corners:
        for k in range(n + 1):
            a = math.radians(a0 + 90.0 * k / n)
            pts.append((cx + rx * math.cos(a), cy + ry * math.sin(a)))
    return pts


def _fmt_bp(bp: int) -> str:
    if bp >= 1_000_000 and bp % 1_000_000 == 0:
        return f"{bp // 1_000_000} Mb"
    if bp >= 1000 and bp % 1000 == 0:
        return f"{bp // 1000} kb"
    return f"{bp:,} bp"


def select_transcripts(txs: List[Transcript], gene: str, mode: str,
                       explicit: Optional[List[str]]) -> List[Transcript]:
    """Filter the transcripts pulled from a region down to the ones we draw."""
    if explicit:
        wanted = {w.upper() for w in explicit}
        keep = [t for t in txs
                if t.name.upper() in wanted
                or t.name.split(".")[0].upper() in wanted]
        if not keep:
            raise UCSCError(
                f"none of {explicit} found here; available: "
                + ", ".join(sorted({t.name for t in txs})[:25])
            )
        return keep

    same_gene = [t for t in txs if t.gene.upper() == gene.upper()] or txs
    same_gene.sort(key=lambda t: (-t.spliced_length, -len(t.exons), t.name))

    if mode == "all":
        return same_gene
    if mode == "canonical":
        return same_gene[:1]
    # "longest-coding": prefer a coding model
    coding = [t for t in same_gene if t.coding]
    return (coding or same_gene)[:1]


# --------------------------------------------------------------------------- #
#  Notebook / library API
# --------------------------------------------------------------------------- #

class GeneFigure:
    """
    Result of :func:`draw_gene`.

    In a Jupyter cell it displays itself as a crisp inline SVG.  It also
    carries the parsed models and the UCSC links, and can save to any format::

        g = draw_gene("TP53")
        g                       # renders inline
        g.transcripts[0].exons  # the underlying model
        g.sequence_url          # UCSC /getSequence link for the locus
        g.save("tp53.svg", "tp53.tif")
    """

    def __init__(self, cartoon: GeneCartoon, transcripts: List[Transcript],
                 links: Dict[str, str], track: str, client: "UCSCClient"):
        self.cartoon = cartoon
        self.transcripts = transcripts
        self.links = links
        self.track = track
        self.gene = cartoon.gene
        self.genome = cartoon.genome
        self._client = client
        self.variants: List[Variant] = list(cartoon.variants)
        self.failed_variants: List[Variant] = []

    @property
    def mapper(self) -> CDNAMapper:
        """cDNA <-> genomic mapper for the primary transcript."""
        return CDNAMapper(self.transcripts[0])

    def variant_table(self) -> List[Dict[str, Any]]:
        """Placed and failed variants as dicts, ready for a DataFrame."""
        rows = []
        for v in list(self.variants) + list(self.failed_variants):
            rows.append({"label": v.label, "cDNA": v.cdna,
                         "genomic": v.genomic, "exon": v.exon,
                         "codon": v.protein, "category": v.category,
                         "count": v.count, "status": v.error or "ok"})
        return rows

    # -- display ----------------------------------------------------------- #

    def _repr_svg_(self) -> str:
        return self.cartoon.to_svg()

    def __repr__(self) -> str:
        return (f"<GeneFigure {self.gene} ({self.genome}, {self.track}): "
                f"{len(self.transcripts)} transcript(s)>")

    # -- convenience ------------------------------------------------------- #

    @property
    def sequence_url(self) -> str:
        return self.links["sequence"]

    @property
    def browser_url(self) -> str:
        return self.links["browser"]

    def sequence(self) -> str:
        """Download the locus sequence (5'->3' for the transcribed strand)."""
        t = self.transcripts[0]
        g0 = min(x.tx_start for x in self.transcripts)
        g1 = max(x.tx_end for x in self.transcripts)
        return self._client.get_sequence(t.chrom, g0, g1,
                                         rev_comp=(t.strand == "-"))

    def save(self, *paths: str, quiet: bool = False) -> List[str]:
        return self.cartoon.render(list(paths), quiet=quiet)

    def to_bytes(self, fmt: str = "png") -> bytes:
        return self.cartoon.to_bytes(fmt)

    def to_svg(self) -> str:
        return self.cartoon.to_svg()

    def table(self) -> List[Dict[str, Any]]:
        """Transcript summary as a list of dicts (feed straight to pandas)."""
        return [{"transcript": t.name, "gene": t.gene, "chrom": t.chrom,
                 "start": t.tx_start + 1, "end": t.tx_end, "strand": t.strand,
                 "exons": len(t.exons), "exonic_bp": t.spliced_length,
                 "cds_bp": sum(e - s for s, e, k in t.segments() if k == "cds"),
                 "coding": t.coding}
                for t in self.transcripts]


def draw_gene(gene: Optional[str] = None, genome: str = "hg38",
              track: Optional[str] = None, region: Optional[str] = None,
              transcripts: str = "longest-coding",
              transcript_ids: Optional[List[str]] = None,
              style: Optional[Any] = None,
              annotations: Optional[List[Dict[str, Any]]] = None,
              variants: Optional[Sequence[Variant]] = None,
              exon: Optional[Any] = None, flank: int = 200,
              quiet: bool = True, **style_kwargs: Any) -> GeneFigure:
    """
    One-liner for notebooks and scripts.

        draw_gene("TP53")
        draw_gene("CFTR", intron_mode="linear", figure_width_in=12)
        draw_gene("TP53", track="ncbiRefSeqCurated", transcripts="all")

    `style` may be a Style object or a path to a style JSON file; any extra
    keyword arguments override individual style fields.
    """
    client = UCSCClient(genome, quiet=quiet)

    if region:
        m = re.match(r"^(\w+):([\d,]+)-([\d,]+)$", region)
        if not m:
            raise UCSCError("region must look like chr17:7668421-7687490")
        chrom = m.group(1)
        start = int(m.group(2).replace(",", "")) - 1
        end = int(m.group(3).replace(",", ""))
    elif gene:
        chrom, start, end = client.locate_symbol(gene)
    else:
        raise UCSCError("pass either a gene symbol or region=")

    if track:
        candidates = [track]
    else:
        available = set(client.gene_pred_tracks())
        candidates = [t for t in GENE_TRACKS if t in available] or ["ncbiRefSeq"]

    txs: List[Transcript] = []
    used = candidates[0]
    for tr in candidates:
        txs = client.transcripts_in_region(tr, chrom, start, end)
        if gene:
            txs = [t for t in txs if t.gene.upper() == gene.upper()] or txs
        if txs:
            used = tr
            break
    if not txs:
        raise UCSCError(f"no transcript models in {chrom}:{start + 1}-{end} "
                        f"(tried {', '.join(candidates)})")

    label = gene or txs[0].gene
    drawn = select_transcripts(txs, label, transcripts, transcript_ids)

    st = style if isinstance(style, Style) else Style.load(style)
    known = {f.name for f in fields(Style)}
    for k, v in style_kwargs.items():
        if k not in known:
            raise TypeError(f"unknown style option {k!r}")
        setattr(st, k, v)

    g0 = min(t.tx_start for t in drawn)
    g1 = max(t.tx_end for t in drawn)
    links = {"sequence": client.sequence_api_url(chrom, g0, g1),
             "browser": client.browser_url(chrom, g0, g1)}

    placed, failed = [], []
    if variants:
        placed, failed = place_variants(list(variants), drawn[0])

    bounds, focus_label = resolve_focus(drawn[0], exon, flank)
    cartoon = GeneCartoon(drawn, st, label, genome, used,
                          annotations=annotations, links=links,
                          variants=placed, focus=bounds,
                          focus_label=focus_label)
    fig = GeneFigure(cartoon, drawn, links, used, client)
    fig.variants, fig.failed_variants = placed, failed
    return fig


def resolve_focus(tx: Transcript, exon: Optional[Any] = None,
                  flank: int = 200) -> Tuple[Optional[Tuple[int, int]], str]:
    """
    Work out a focus window from an exon spec.

    ``exon`` may be a number (``14``), a string range (``"14-16"``), a
    ``(first, last)`` pair, or None for the whole gene.  Returns
    ``(bounds, label)`` ready for :class:`GeneCartoon`.
    """
    if exon is None or exon == "" or str(exon).lower() in ("all", "whole gene"):
        return None, ""

    if isinstance(exon, (tuple, list)) and len(exon) == 2:
        first, last = int(exon[0]), int(exon[1])
    else:
        s = str(exon).strip().lower().replace("exon", "").strip()
        m = re.match(r"^(\d+)\s*[-–:]\s*(\d+)$", s)
        if m:
            first, last = int(m.group(1)), int(m.group(2))
        else:
            if not s.isdigit():
                raise ValueError(
                    f"cannot read {exon!r} as an exon; use e.g. 14 or 14-16")
            first = last = int(s)

    bounds = tx.exon_window(first, last, flank=flank)
    label = f"exon {first}" if first == last else f"exons {first}–{last}"
    return bounds, label


def load_annotations(path: Optional[str]) -> List[Dict[str, Any]]:
    if not path:
        return []
    data = json.loads(Path(path).read_text())
    if isinstance(data, dict):
        data = data.get("annotations", [])
    return list(data)


# --------------------------------------------------------------------------- #
#  CLI
# --------------------------------------------------------------------------- #

EPILOG = textwrap.dedent("""\
    examples
    --------
      # MANE Select model of TP53, compressed introns, SVG for Illustrator
      %(prog)s TP53 -o figures/tp53.svg

      # every RefSeq isoform of CFTR, true genomic scale, 600-dpi TIFF + SVG
      %(prog)s CFTR --track ncbiRefSeqCurated --transcripts all \\
               --intron-mode linear -o cftr.svg -o cftr.tif --dpi 600

      # mouse gene, explicit transcripts, custom look, domain annotations
      %(prog)s Trp53 --genome mm39 --transcript NM_011640 \\
               --style lab_style.json --annotations domains.json -o trp53.pdf

      # start from the defaults and tweak
      %(prog)s --dump-style lab_style.json

    annotation file (JSON list)
    ---------------------------
      [{"start": 7676594, "end": 7676000, "label": "DNA-binding domain",
        "color": "#E8A33D", "style": "box", "coords": "genomic"},
       {"start": 7674220, "label": "R248Q", "style": "marker",
        "color": "#D1495B"}]

      coords: "genomic" (0-based, UCSC internal) or "browser" (1-based)
      style:  box | bracket | marker
    """)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="ucsc_gene_cartoon.py",
        description="Draw a publication-quality gene structure cartoon from "
                    "UCSC Genome Browser REST API data.",
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("gene", nargs="?", help="gene symbol, e.g. TP53")

    src = p.add_argument_group("data source")
    src.add_argument("-g", "--genome", default="hg38",
                     help="UCSC assembly (default: hg38)")
    src.add_argument("-t", "--track", default=None,
                     help="annotation track (default: first available of "
                          + ", ".join(GENE_TRACKS[:4]) + ")")
    src.add_argument("--region", default=None,
                     help="draw a region instead of resolving a symbol, "
                          "chrN:start-end (1-based, as in the browser)")
    src.add_argument("--transcripts", choices=["canonical", "longest-coding", "all"],
                     default="longest-coding",
                     help="which isoforms to draw (default: longest-coding)")
    src.add_argument("--transcript", action="append", dest="transcript_ids",
                     help="draw a specific transcript by accession; repeatable")
    src.add_argument("--exon", default=None, metavar="N|A-B",
                     help="zoom in on one exon (--exon 14) or a run of them "
                          "(--exon 14-16); numbers are biological, 5'->3'")
    src.add_argument("--flank", type=int, default=200, metavar="BP",
                     help="intron shown either side of a focused exon "
                          "(default: 200 bp)")
    src.add_argument("--no-cache", action="store_true",
                     help="bypass the local API response cache")
    src.add_argument("--cache-dir", default=None)

    out = p.add_argument_group("output")
    out.add_argument("-o", "--out", action="append", default=None,
                     help="output file; extension picks the format "
                          "(.svg .pdf .eps .png .tif). Repeatable.")
    out.add_argument("--dpi", type=int, default=None,
                     help="raster resolution (default 400)")
    out.add_argument("--width", type=float, default=None,
                     help="figure width in inches (default 9)")
    out.add_argument("--save-sequence", metavar="FASTA", default=None,
                     help="also download the locus sequence to a FASTA file")
    out.add_argument("--save-json", metavar="JSON", default=None,
                     help="dump the parsed transcript models")

    sty = p.add_argument_group("style")
    sty.add_argument("--style", default=None, help="style JSON file")
    sty.add_argument("--dump-style", metavar="FILE", default=None,
                     help="write the default style JSON and exit")
    sty.add_argument("--intron-mode", choices=["linear", "compress", "equal"],
                     default=None,
                     help="linear = true scale; compress = fixed-width introns "
                          "(default); equal = equal-width exons too")
    sty.add_argument("--intron-style", choices=["line", "chevron_line", "angled"],
                     default=None)
    sty.add_argument("--exon-shape", choices=["box", "round"], default=None)
    sty.add_argument("--cds-color", default=None)
    sty.add_argument("--utr-color", default=None)
    sty.add_argument("--axis", dest="axis_mode",
                     choices=["auto", "coords", "scalebar", "none"], default=None)
    sty.add_argument("--no-exon-numbers", action="store_true")
    sty.add_argument("--no-legend", action="store_true")
    sty.add_argument("--no-title", action="store_true")

    ann = p.add_argument_group("annotations")
    ann.add_argument("--annotations", default=None,
                     help="JSON file of extra features to overlay")

    misc = p.add_argument_group("misc")
    misc.add_argument("--list-tracks", action="store_true",
                      help="list transcript tracks in the assembly and exit")
    misc.add_argument("-q", "--quiet", action="store_true")
    return p


def apply_overrides(style: Style, a: argparse.Namespace) -> Style:
    direct = ["intron_mode", "intron_style", "exon_shape", "cds_color",
              "utr_color", "axis_mode", "dpi"]
    for k in direct:
        v = getattr(a, k, None)
        if v is not None:
            setattr(style, k, v)
    if a.width is not None:
        style.figure_width_in = a.width
    if a.no_exon_numbers:
        style.show_exon_numbers = False
    if a.no_legend:
        style.show_legend = False
    if a.no_title:
        style.show_title = style.show_subtitle = False
    return style


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    if args.dump_style:
        Style().dump(args.dump_style)
        print(f"wrote default style to {args.dump_style}")
        return 0

    client = UCSCClient(args.genome, cache_dir=args.cache_dir,
                        quiet=args.quiet, use_cache=not args.no_cache)

    if args.list_tracks:
        for t in client.gene_pred_tracks():
            print(t)
        return 0

    if not args.gene and not args.region:
        build_parser().print_help()
        return 2

    try:
        # ---- locate ---- #
        if args.region:
            m = re.match(r"^(\w+):([\d,]+)-([\d,]+)$", args.region)
            if not m:
                raise UCSCError("--region must look like chr17:7668421-7687490")
            chrom = m.group(1)
            start = int(m.group(2).replace(",", "")) - 1
            end = int(m.group(3).replace(",", ""))
        else:
            chrom, start, end = client.locate_symbol(args.gene)

        # ---- track ---- #
        if args.track:
            tracks = [args.track]
        else:
            available = set(client.gene_pred_tracks())
            tracks = [t for t in GENE_TRACKS if t in available] or ["ncbiRefSeq"]

        txs: List[Transcript] = []
        used_track = tracks[0]
        for tr in tracks:
            txs = client.transcripts_in_region(tr, chrom, start, end)
            if args.gene:
                txs = [t for t in txs if t.gene.upper() == args.gene.upper()] or txs
            if txs:
                used_track = tr
                break
        if not txs:
            raise UCSCError(
                f"no transcript models found in {chrom}:{start + 1}-{end} "
                f"(tried {', '.join(tracks)}). Try --track / --list-tracks."
            )

        gene_label = args.gene or txs[0].gene
        drawn = select_transcripts(txs, gene_label, args.transcripts,
                                   args.transcript_ids)

        if not args.quiet:
            print(f"{gene_label}: {len(drawn)} transcript(s) from {used_track}",
                  file=sys.stderr)
            for t in drawn:
                print(f"  {t.name:<18} {t.chrom}:{t.tx_start + 1:,}-{t.tx_end:,} "
                      f"{t.strand}  {len(t.exons)} exons  "
                      f"{t.spliced_length:,} nt spliced"
                      f"{'' if t.coding else '  [non-coding]'}", file=sys.stderr)

        # ---- links ---- #
        g0 = min(t.tx_start for t in drawn)
        g1 = max(t.tx_end for t in drawn)
        links = {
            "sequence": client.sequence_api_url(chrom, g0, g1),
            "browser": client.browser_url(chrom, g0, g1),
        }
        if not args.quiet:
            print(f"  sequence: {links['sequence']}", file=sys.stderr)
            print(f"  browser : {links['browser']}", file=sys.stderr)

        # ---- optional extras ---- #
        if args.save_sequence:
            dna = client.get_sequence(chrom, g0, g1,
                                      rev_comp=(drawn[0].strand == "-"))
            hdr = (f">{gene_label} {args.genome} {chrom}:{g0 + 1}-{g1} "
                   f"strand={drawn[0].strand}")
            body = "\n".join(dna[i:i + 60] for i in range(0, len(dna), 60))
            Path(args.save_sequence).write_text(f"{hdr}\n{body}\n")
            print(f"wrote {args.save_sequence} ({len(dna):,} nt)")

        if args.save_json:
            Path(args.save_json).write_text(
                json.dumps({"gene": gene_label, "genome": args.genome,
                            "track": used_track, "links": links,
                            "transcripts": [asdict(t) for t in drawn]}, indent=2))
            print(f"wrote {args.save_json}")

        # ---- draw ---- #
        style = apply_overrides(Style.load(args.style), args)
        bounds, focus_label = resolve_focus(drawn[0], args.exon, args.flank)
        if bounds and not args.quiet:
            print(f"  focus  : {focus_label} "
                  f"({chrom}:{bounds[0] + 1:,}-{bounds[1]:,})", file=sys.stderr)
        outs = args.out or [f"{gene_label}_{args.genome}.svg"]
        GeneCartoon(drawn, style, gene_label, args.genome, used_track,
                    annotations=load_annotations(args.annotations),
                    links=links, focus=bounds,
                    focus_label=focus_label).render(outs)

    except (UCSCError, ValueError, OSError, json.JSONDecodeError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
