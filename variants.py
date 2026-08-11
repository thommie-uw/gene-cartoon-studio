#!/usr/bin/env python3
"""
variants.py -- read variant tables (CSV / TSV / Excel) for gene cartoons.

The point is that a user drops in whatever spreadsheet they already have and
it works: column names are detected case- and punctuation-insensitively, and
anything that can't be placed comes back with a per-row reason rather than
blowing up the whole import.

    from variants import load_variant_table
    vf = load_variant_table("my_variants.xlsx")
    vf.variants        # list[Variant], ready for draw_gene(variants=...)
    vf.columns         # which column was used for what
    vf.warnings        # anything the user should know
"""

from __future__ import annotations

import io
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union

try:
    import pandas as pd
except ImportError:  # pragma: no cover
    raise SystemExit("variants.py needs pandas:  pip install pandas openpyxl")

from ucsc_gene_cartoon import CDNAError, Variant, parse_hgvs

# --------------------------------------------------------------------------- #
#  Column detection
# --------------------------------------------------------------------------- #

#: role -> candidate header names, in priority order.  Matching is done on a
#: normalised form of the header (lowercased, punctuation stripped).
COLUMN_HINTS: Dict[str, List[str]] = {
    "cdna": [
        "cdna", "cdnaposition", "cdnachange", "chgvs", "hgvsc", "hgvscdna",
        "cdot", "cnomenclature", "nucleotide", "nucleotidechange", "cdnavariant",
        "variantcdna", "hgvs", "variant", "mutation", "position", "cpos",
    ],
    "label": [
        "label", "name", "variantname", "protein", "proteinchange", "phgvs",
        "hgvsp", "pdot", "aachange", "aminoacidchange", "aminoacid",
        "proteinvariant", "id", "variantid",
    ],
    "category": [
        "category", "class", "classification", "clinicalsignificance",
        "significance", "clinsig", "consequence", "effect", "type",
        "varianttype", "group", "impact", "acmg", "pathogenicity",
    ],
    "count": [
        "count", "n", "occurrences", "occurrence", "recurrence", "frequency",
        "freq", "cases", "samples", "patients", "observations", "tally",
    ],
    "color": ["color", "colour", "hex", "hexcolor", "hexcolour"],
    "genomic": [
        "genomic", "genomicposition", "gpos", "gdot", "hgvsg", "chrompos",
        "pos", "start", "coordinate", "grch38", "hg38",
    ],
    "note": ["note", "notes", "comment", "comments", "description", "source"],
}

#: never auto-detect these as the cDNA column even though they look plausible
_CDNA_LAST_RESORT = {"position", "variant", "mutation", "hgvs", "cpos"}

#: substrings used for a second, fuzzier pass over headers we didn't match
#: exactly -- catches things like "No. of cases" or "Variant classification"
FUZZY_HINTS: Dict[str, List[str]] = {
    "count": ["case", "sample", "patient", "count", "occur", "freq",
              "recurr", "observ", "tally", "tumour", "tumor"],
    "category": ["signif", "classif", "consequence", "pathogen", "effect",
                 "impact", "acmg", "categor"],
    "label": ["protein", "aminoacid", "aachange", "mutationaa", "variantname",
              "label", "aacode"],
    "cdna": ["cdna", "hgvsc", "nucleotide", "cchange"],
    "note": ["note", "comment", "descript"],
}

#: a value only counts as cDNA-looking if it carries real HGVS markers --
#: a bare integer is far more likely to be an ID or a count
_CDNA_SHAPE = re.compile(
    r"(^|:)\s*[cnr]\.|[+\-]\d|\*\d|\d\s*[ACGT]\s*>\s*[ACGT]|del|dup|ins",
    re.I,
)


def _norm(name: Any) -> str:
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


def detect_columns(df: "pd.DataFrame") -> Dict[str, Optional[str]]:
    """Guess which column plays which role.  Returns ``{role: column or None}``."""
    normalised = {_norm(c): c for c in df.columns}
    out: Dict[str, Optional[str]] = {r: None for r in COLUMN_HINTS}
    taken: set = set()

    for role, hints in COLUMN_HINTS.items():
        for hint in hints:
            if role == "cdna" and hint in _CDNA_LAST_RESORT:
                continue
            col = normalised.get(hint)
            if col is not None and col not in taken:
                out[role] = col
                taken.add(col)
                break

    # second pass: fuzzy substring match on anything still unfilled
    for role, needles in FUZZY_HINTS.items():
        if out.get(role) is not None:
            continue
        for col in df.columns:
            if col in taken:
                continue
            n = _norm(col)
            if not any(needle in n for needle in needles):
                continue
            # a recurrence column has to actually hold numbers, or "Sample ID"
            # gets mistaken for a count
            if role == "count" and not _mostly_numeric(df[col]):
                continue
            out[role] = col
            taken.add(col)
            break

    # cDNA fallbacks: vague headers, then any column that *looks* like c. notation
    if out["cdna"] is None:
        for hint in _CDNA_LAST_RESORT:
            col = normalised.get(hint)
            if col is not None and col not in taken:
                out["cdna"] = col
                taken.add(col)
                break
    if out["cdna"] is None:
        for col in df.columns:
            if col in taken:
                continue
            if _looks_like_cdna(df[col]):
                out["cdna"] = col
                taken.add(col)
                break
    return out


def _mostly_numeric(series: "pd.Series", threshold: float = 0.8) -> bool:
    sample = [s for s in series.dropna().head(25)]
    if not sample:
        return False
    ok = 0
    for s in sample:
        try:
            float(str(s).strip().replace(",", ""))
            ok += 1
        except ValueError:
            pass
    return ok / len(sample) >= threshold


def _looks_like_cdna(series: "pd.Series", threshold: float = 0.6) -> bool:
    sample = [s for s in series.dropna().astype(str).head(25) if s.strip()]
    if not sample:
        return False
    hits = 0
    for s in sample:
        if not _CDNA_SHAPE.search(s):
            continue          # bare integers etc. are not evidence
        try:
            parse_hgvs(s)
            hits += 1
        except CDNAError:
            pass
    return hits / len(sample) >= threshold


# --------------------------------------------------------------------------- #
#  Loading
# --------------------------------------------------------------------------- #

@dataclass
class VariantFile:
    """A parsed variant table."""
    variants: List[Variant] = field(default_factory=list)
    columns: Dict[str, Optional[str]] = field(default_factory=dict)
    dataframe: Optional["pd.DataFrame"] = None
    warnings: List[str] = field(default_factory=list)
    source: str = ""
    sheet: Optional[str] = None

    @property
    def categories(self) -> List[str]:
        seen: List[str] = []
        for v in self.variants:
            if v.category and v.category not in seen:
                seen.append(v.category)
        return seen

    def summary(self) -> str:
        used = ", ".join(f"{r}={c}" for r, c in self.columns.items() if c)
        return (f"{len(self.variants)} variant(s) from {self.source or 'table'}"
                + (f" [{used}]" if used else ""))


def read_table(source: Union[str, Path, io.IOBase], sheet: Optional[Any] = None
               ) -> "pd.DataFrame":
    """Read CSV/TSV/TXT/Excel into a DataFrame, sniffing the delimiter."""
    name = getattr(source, "name", None) or str(source)
    suffix = Path(str(name)).suffix.lower()

    if suffix in (".xlsx", ".xlsm", ".xls", ".xltx"):
        try:
            return pd.read_excel(source, sheet_name=sheet if sheet is not None else 0)
        except ImportError as e:  # pragma: no cover
            raise RuntimeError(
                "reading Excel needs openpyxl:  pip install openpyxl"
            ) from e
    # NB: do *not* use pandas' sep=None sniffer.  On a single-column file it
    # happily decides the delimiter is "c" and shreds "c.743G>A".  Count real
    # candidates in the header line instead.
    if suffix in (".tsv", ".tab"):
        sep = "\t"
    else:
        sep = _sniff_delimiter(source)
    return pd.read_csv(source, sep=sep, engine="python",
                       comment="#", skip_blank_lines=True)


def _sniff_delimiter(source: Union[str, Path, io.IOBase],
                     candidates: str = ",\t;|") -> str:
    """Pick the delimiter by counting candidates in the first non-comment line."""
    header = ""
    try:
        if hasattr(source, "read"):
            pos = source.tell()
            raw = source.read(8192)
            source.seek(pos)
            if isinstance(raw, bytes):
                raw = raw.decode("utf-8", "replace")
        else:
            with open(source, "r", encoding="utf-8-sig", errors="replace") as fh:
                raw = fh.read(8192)
        for line in raw.splitlines():
            if line.strip() and not line.lstrip().startswith("#"):
                header = line
                break
    except (OSError, UnicodeDecodeError):
        return ","
    counts = {c: header.count(c) for c in candidates}
    best = max(counts, key=lambda c: counts[c])
    return best if counts[best] > 0 else ","


def excel_sheet_names(source: Union[str, Path, io.IOBase]) -> List[str]:
    """Sheet names in an Excel workbook (empty list for CSV)."""
    name = getattr(source, "name", None) or str(source)
    if Path(str(name)).suffix.lower() not in (".xlsx", ".xlsm", ".xls", ".xltx"):
        return []
    try:
        return list(pd.ExcelFile(source).sheet_names)
    except Exception:  # pragma: no cover
        return []


def build_variants(df: "pd.DataFrame", columns: Dict[str, Optional[str]],
                   default_category: str = "") -> "VariantFile":
    """Turn a DataFrame plus a column mapping into Variant objects."""
    vf = VariantFile(columns=dict(columns), dataframe=df)

    if not columns.get("cdna") and not columns.get("genomic"):
        vf.warnings.append(
            "No cDNA or genomic position column found. Expected a column of "
            "HGVS c. positions such as 'c.743G>A' — check the header row, or "
            "pick the column manually."
        )
        return vf

    col_c = columns.get("cdna")
    col_g = columns.get("genomic")
    blank = skipped = 0

    for i, row in df.iterrows():
        raw_c = row[col_c] if col_c and col_c in df.columns else None
        raw_g = row[col_g] if col_g and col_g in df.columns else None
        if _empty(raw_c) and _empty(raw_g):
            blank += 1
            continue

        v = Variant(
            cdna="" if _empty(raw_c) else str(raw_c).strip(),
            label=_cell(row, columns.get("label"), df),
            category=_cell(row, columns.get("category"), df) or default_category,
            note=_cell(row, columns.get("note"), df),
            count=_int(_cell(row, columns.get("count"), df), 1),
            color=_cell(row, columns.get("color"), df) or None,
        )
        if not _empty(raw_g):
            g = _int(str(raw_g).replace(",", ""), None)
            # accept "chr17:7674220" as well as a bare number
            if g is None:
                m = re.search(r"(\d[\d,]*)\s*$", str(raw_g))
                g = _int(m.group(1).replace(",", ""), None) if m else None
            v.genomic = g
            if g is None and not v.cdna:
                v.error = f"row {i + 2}: cannot read position {raw_g!r}"
                skipped += 1
        if not v.label:
            v.label = v.cdna or (str(v.genomic) if v.genomic else "")
        vf.variants.append(v)

    if blank:
        vf.warnings.append(f"skipped {blank} row(s) with no position")
    if skipped:
        vf.warnings.append(f"{skipped} row(s) had an unreadable position")
    return vf


def load_variant_table(source: Union[str, Path, io.IOBase],
                       sheet: Optional[Any] = None,
                       columns: Optional[Dict[str, Optional[str]]] = None,
                       default_category: str = "") -> VariantFile:
    """
    Read a variant table and return placed-ready :class:`Variant` objects.

    Pass ``columns`` to override the automatic column detection, e.g.
    ``{"cdna": "HGVSc", "label": "Protein", "category": "ClinSig"}``.
    """
    df = read_table(source, sheet=sheet)
    df = df.dropna(how="all")
    df.columns = [str(c).strip() for c in df.columns]
    detected = detect_columns(df)
    if columns:
        detected.update({k: v for k, v in columns.items()})
    vf = build_variants(df, detected, default_category=default_category)
    vf.source = str(getattr(source, "name", None) or source)
    vf.sheet = sheet if isinstance(sheet, str) else None
    return vf


# --------------------------------------------------------------------------- #
#  small helpers
# --------------------------------------------------------------------------- #

def _empty(x: Any) -> bool:
    if x is None:
        return True
    if isinstance(x, float) and x != x:      # NaN
        return True
    return str(x).strip() == ""


def _cell(row: Any, col: Optional[str], df: "pd.DataFrame") -> str:
    if not col or col not in df.columns:
        return ""
    val = row[col]
    return "" if _empty(val) else str(val).strip()


def _int(x: Any, default: Any) -> Any:
    if _empty(x):
        return default
    try:
        return int(float(str(x).strip()))
    except (TypeError, ValueError):
        return default


def write_template(path: Union[str, Path]) -> Path:
    """Write an example variant file people can fill in."""
    rows = [
        {"cDNA": "c.524G>A", "Protein": "R175H",
         "Classification": "Pathogenic", "Count": 9,
         "Notes": "hotspot, DNA-binding domain"},
        {"cDNA": "c.659A>G", "Protein": "Y220C",
         "Classification": "Pathogenic", "Count": 5, "Notes": ""},
        {"cDNA": "c.743G>A", "Protein": "R248Q",
         "Classification": "Pathogenic", "Count": 12, "Notes": "hotspot"},
        {"cDNA": "c.818G>A", "Protein": "R273H",
         "Classification": "Pathogenic", "Count": 8, "Notes": "hotspot"},
        {"cDNA": "c.215C>G", "Protein": "P72R",
         "Classification": "Benign", "Count": 1, "Notes": "common polymorphism"},
        {"cDNA": "c.376-2A>G", "Protein": "splice acceptor",
         "Classification": "VUS", "Count": 1, "Notes": "intronic"},
        {"cDNA": "c.-28C>T", "Protein": "5'UTR",
         "Classification": "VUS", "Count": 1, "Notes": "upstream of ATG"},
        {"cDNA": "c.*38G>A", "Protein": "3'UTR",
         "Classification": "VUS", "Count": 2, "Notes": "after stop codon"},
    ]
    df = pd.DataFrame(rows)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() in (".xlsx", ".xlsm"):
        df.to_excel(path, index=False, sheet_name="variants")
    else:
        df.to_csv(path, index=False)
    return path


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--template":
        out = write_template(sys.argv[2] if len(sys.argv) > 2
                             else "variant_template.xlsx")
        print(f"wrote {out}")
    elif len(sys.argv) > 1:
        vf = load_variant_table(sys.argv[1])
        print(vf.summary())
        for w in vf.warnings:
            print("  !", w)
        for v in vf.variants[:20]:
            print(f"  {v.label:<18} {v.cdna:<14} {v.category:<14} n={v.count}")
    else:
        print(__doc__)
