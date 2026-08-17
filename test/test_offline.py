#!/usr/bin/env python3
"""Offline regression test: real UCSC payloads, no network.

Verifies transcript parsing, coordinate mapping and rendering against
genePred rows captured live from api.genome.ucsc.edu (hg38, 2025-08-13).
"""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import ucsc_gene_cartoon as ugc

OUT = Path(__file__).resolve().parent / "out"
OUT.mkdir(exist_ok=True)

# --- real rows captured from the UCSC REST API ---------------------------- #
TP53_MANE = {
    "name": "NM_000546.6", "chrom": "chr17", "strand": "-",
    "txStart": 7668420, "txEnd": 7687490,
    "cdsStart": 7669608, "cdsEnd": 7676594, "exonCount": 11,
    "exonStarts": "7668420,7670608,7673534,7673700,7674180,7674858,7675052,7675993,7676381,7676520,7687376,",
    "exonEnds": "7669690,7670715,7673608,7673837,7674290,7674971,7675236,7676272,7676403,7676622,7687490,",
    "name2": "TP53",
}
TP53_ISO2 = {
    "name": "NM_001126112.3", "chrom": "chr17", "strand": "-",
    "txStart": 7668420, "txEnd": 7687490,
    "cdsStart": 7669608, "cdsEnd": 7676594, "exonCount": 11,
    "exonStarts": "7668420,7670608,7673534,7673700,7674180,7674858,7675052,7675993,7676381,7676520,7687376,",
    "exonEnds": "7669690,7670715,7673608,7673837,7674290,7674971,7675236,7676272,7676403,7676619,7687490,",
    "name2": "TP53",
}
TP53_ISO3 = {
    "name": "NM_001126116.2", "chrom": "chr17", "strand": "-",
    "txStart": 7668420, "txEnd": 7675244,
    "cdsStart": 7673306, "cdsEnd": 7675215, "exonCount": 8,
    "exonStarts": "7668420,7670608,7673206,7673534,7673700,7674180,7674858,7675052,",
    "exonEnds": "7669690,7670715,7673339,7673608,7673837,7674290,7674971,7675244,",
    "name2": "TP53",
}
TP53_NC = {
    "name": "NR_176326.1", "chrom": "chr17", "strand": "-",
    "txStart": 7668420, "txEnd": 7687490,
    "cdsStart": 7687490, "cdsEnd": 7687490, "exonCount": 10,
    "exonStarts": "7668420,7670608,7673534,7673700,7674180,7675052,7675993,7676381,7676520,7687376,",
    "exonEnds": "7669690,7670715,7673608,7673837,7674290,7675236,7676272,7676403,7676622,7687490,",
    "name2": "TP53",
}
CFTR = {
    "name": "NM_000492.4", "chrom": "chr7", "strand": "+",
    "txStart": 117480024, "txEnd": 117668665,
    "cdsStart": 117480094, "cdsEnd": 117667108, "exonCount": 27,
    "exonStarts": "117480024,117504252,117509033,117530898,117534275,117535247,117536547,117540099,117542015,117548640,117559463,117587738,117590352,117591933,117594929,117602825,117603531,117606673,117610518,117611580,117614612,117627521,117642437,117652841,117664687,117665458,117666907,",
    "exonEnds": "117480147,117504363,117509142,117531114,117534365,117535411,117536673,117540346,117542108,117548823,117559655,117587833,117590439,117592657,117595058,117602863,117603782,117606753,117610669,117611808,117614713,117627770,117642593,117652931,117664860,117665564,117668665,",
    "name2": "CFTR",
}
XIST = {
    "name": "NR_001564.3", "chrom": "chrX", "strand": "-",
    "txStart": 73820655, "txEnd": 73852714,
    "cdsStart": 73852714, "cdsEnd": 73852714, "exonCount": 6,
    "exonStarts": "73820655,73829067,73831065,73833237,73837439,73841381,",
    "exonEnds": "73827984,73829231,73831274,73833374,73837503,73852714,",
    "name2": "XIST",
}

failures = []


def check(name, cond, detail=""):
    print(f"{'PASS' if cond else 'FAIL'}  {name}{'  ' + detail if detail else ''}")
    if not cond:
        failures.append(name)


# --- 1. parsing ----------------------------------------------------------- #
t = ugc.Transcript.from_genepred(TP53_MANE)
check("TP53 exon count", len(t.exons) == 11, f"got {len(t.exons)}")
check("TP53 strand/coding", t.strand == "-" and t.coding)
# 2512 = sum of genomic exon blocks in UCSC's alignment of NM_000546.6
# (the RefSeq mRNA record is 2591 nt; the difference is the poly-A tail).
check("TP53 exonic genomic length 2512 nt",
      t.spliced_length == 2512, f"got {t.spliced_length}")
check("TP53 CDS length divisible by 3",
      sum(e - s for s, e, k in t.segments() if k == "cds") % 3 == 0,
      f"cds={sum(e - s for s, e, k in t.segments() if k == 'cds')}")
check("TP53 CDS = 1182 nt (p53, 393 aa + stop)",
      sum(e - s for s, e, k in t.segments() if k == "cds") == 1182)
check("TP53 introns = exons - 1", len(t.introns) == 10)
check("TP53 minus-strand exon numbering: first block is exon 11",
      t.exon_number(0) == 11 and t.exon_number(10) == 1)

c = ugc.Transcript.from_genepred(CFTR)
check("CFTR exon count 27", len(c.exons) == 27, f"got {len(c.exons)}")
check("CFTR CDS = 4443 nt (CFTR, 1480 aa + stop)",
      sum(e - s for s, e, k in c.segments() if k == "cds") == 4443)
check("CFTR plus-strand numbering", c.exon_number(0) == 1)
check("CFTR span 188.6 kb", 188_000 < c.length < 189_000, f"{c.length:,} bp")

x = ugc.Transcript.from_genepred(XIST)
check("XIST is non-coding", not x.coding)
check("XIST segments all 'nc'", {k for _, _, k in x.segments()} == {"nc"})

# --- 2. UTR/CDS split integrity ------------------------------------------ #
for tx, label in ((t, "TP53"), (c, "CFTR"), (x, "XIST")):
    segs = tx.segments()
    check(f"{label} segments tile exons exactly",
          sum(e - s for s, e, _ in segs) == tx.spliced_length)
    check(f"{label} segments non-overlapping",
          all(segs[i][1] <= segs[i + 1][0] for i in range(len(segs) - 1)))
    if tx.coding:
        kinds = [k for _, _, k in segs]
        check(f"{label} has 5'UTR, CDS and 3'UTR",
              {"utr5", "cds", "utr3"} <= set(kinds))

# --- 3. coordinate mapper ------------------------------------------------- #
st = ugc.Style()
for mode in ("linear", "compress", "equal"):
    st.intron_mode = mode
    m = ugc.CoordinateMapper([c], st)
    xs = [m.x(p) for p in range(c.tx_start, c.tx_end, 500)]
    check(f"mapper[{mode}] monotonic", all(a <= b + 1e-12 for a, b in zip(xs, xs[1:])))
    check(f"mapper[{mode}] in [0,1]", min(xs) >= -1e-9 and max(xs) <= 1 + 1e-9)
    check(f"mapper[{mode}] endpoints", abs(m.x(c.tx_start)) < 1e-9 and abs(m.x(c.tx_end) - 1) < 1e-9)

st.intron_mode = "compress"
m = ugc.CoordinateMapper([c], st)
ex_w = sum(m.x(e) - m.x(s) for s, e in c.exons)
check("compress mode gives exons a real share of the axis (>25%)",
      ex_w > 0.25, f"exons occupy {ex_w:.0%} of width")
mlin = ugc.CoordinateMapper([c], ugc.Style(intron_mode="linear"))
ex_lin = sum(mlin.x(e) - mlin.x(s) for s, e in c.exons)
check("linear mode squashes CFTR exons (<5%)", ex_lin < 0.05,
      f"exons occupy {ex_lin:.1%} of width")

mflip = ugc.CoordinateMapper([t], st, flip=True)
check("flip puts minus-strand 5' end on the left",
      mflip.x(t.tx_end) < mflip.x(t.tx_start))

# --- 4. rendering --------------------------------------------------------- #
links = {"sequence": "https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr17;start=7668420;end=7687490",
         "browser": "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr17%3A7668421-7687490"}

cases = [
    ("tp53_default", [t], ugc.Style(), "TP53", "chr17", {}),
    ("tp53_isoforms", [t, ugc.Transcript.from_genepred(TP53_ISO2),
                       ugc.Transcript.from_genepred(TP53_ISO3),
                       ugc.Transcript.from_genepred(TP53_NC)],
     ugc.Style(intron_style="chevron_line"), "TP53", "chr17", {}),
    ("cftr_compress", [c], ugc.Style(figure_width_in=11, exon_number_size=5.5),
     "CFTR", "chr7", {}),
    ("cftr_linear", [c], ugc.Style(intron_mode="linear", show_exon_numbers=False),
     "CFTR", "chr7", {}),
    ("xist_noncoding", [x], ugc.Style(exon_shape="round"), "XIST", "chrX", {}),
    ("tp53_annotated", [t],
     ugc.Style(cds_color="#3D5A6C", utr_color="#C9D6DE", intron_style="angled",
               exon_number_position="inside", background="white"),
     "TP53", "chr17",
     {"annotations": [
         {"start": 7673535, "end": 7676272, "label": "DNA-binding domain",
          "style": "box", "color": "#E8A33D"},
         {"start": 7674220, "label": "R248Q", "style": "marker",
          "color": "#D1495B"}]}),
]

for name, txs, style, gene, chrom, kw in cases:
    try:
        cart = ugc.GeneCartoon(txs, style, gene, "hg38", "ncbiRefSeqSelect",
                               annotations=kw.get("annotations"), links=links)
        svg = OUT / f"{name}.svg"
        png = OUT / f"{name}.png"
        cart.render([str(svg), str(png)])
        body = svg.read_text()
        ok = (svg.stat().st_size > 4000 and png.stat().st_size > 4000
              and body.lstrip().startswith("<?xml") and "</svg>" in body)
        check(f"render {name}", ok,
              f"svg {svg.stat().st_size // 1024} kB, png {png.stat().st_size // 1024} kB")
        if name == "tp53_default":
            check("SVG keeps text as text (editable in Illustrator)",
                  "<text" in body and "font-family" in body)
            check("SVG carries the UCSC sequence/browser hyperlinks",
                  "api.genome.ucsc.edu/getData/sequence" in body
                  and "cgi-bin/hgTracks" in body)
            check("exon numbers 1..11 all present",
                  all(f">{n}<" in body for n in range(1, 12)))
    except Exception as e:  # noqa: BLE001
        check(f"render {name}", False, f"{type(e).__name__}: {e}")

# --- 5. helpers ----------------------------------------------------------- #
check("_fmt_bp", (ugc._fmt_bp(1000), ugc._fmt_bp(2_000_000), ugc._fmt_bp(750))
      == ("1 kb", "2 Mb", "750 bp"))
ticks = ugc._nice_ticks(7668420, 7687490, 5)
check("_nice_ticks inside range and ascending",
      ticks == sorted(ticks) and ticks[0] > 7668420 and ticks[-1] < 7687490,
      str(ticks))
check("style round-trips through JSON",
      ugc.Style.load(None) == ugc.Style()
      and json.loads(json.dumps(ugc.asdict(ugc.Style())))["cds_color"] == "#2C6FA6")

sel = ugc.select_transcripts(
    [ugc.Transcript.from_genepred(r) for r in (TP53_NC, TP53_ISO3, TP53_MANE)],
    "TP53", "longest-coding", None)
check("select_transcripts prefers a coding model",
      len(sel) == 1 and sel[0].coding, sel[0].name)
sel_all = ugc.select_transcripts(
    [ugc.Transcript.from_genepred(r) for r in (TP53_NC, TP53_ISO3, TP53_MANE)],
    "TP53", "all", None)
check("select_transcripts 'all' keeps everything", len(sel_all) == 3)
sel_id = ugc.select_transcripts(
    [ugc.Transcript.from_genepred(r) for r in (TP53_NC, TP53_MANE)],
    "TP53", "all", ["NM_000546"])
check("explicit accession lookup ignores version suffix",
      len(sel_id) == 1 and sel_id[0].name == "NM_000546.6")

# --- 6. end-to-end CLI, with the network layer stubbed --------------------- #
def fake_get(self, endpoint, params, use_cache=True):
    if endpoint == "search":
        return {"positionMatches": [{"trackName": "mane", "matches": [
            {"posName": "TP53 NM_000546.6", "position": "chr17:7668421-7687490"}]}]}
    if endpoint == "list/tracks":
        return {"hg38": {"ncbiRefSeqSelect": {"type": "genePred"},
                         "gap": {"type": "bed 3"}}}
    if endpoint == "getData/track":
        return {"ncbiRefSeqSelect": [dict(TP53_MANE)]}
    if endpoint == "getData/sequence":
        return {"dna": "acgt" * 100}
    raise AssertionError(endpoint)


ugc.UCSCClient._get = fake_get
import io                                        # noqa: E402
import contextlib                                # noqa: E402

buf = io.StringIO()
with contextlib.redirect_stderr(buf), contextlib.redirect_stdout(io.StringIO()):
    rc = ugc.main(["TP53", "-o", str(OUT / "cli_tp53.svg"),
                   "--save-sequence", str(OUT / "cli.fa"),
                   "--save-json", str(OUT / "cli.json")])
check("CLI: symbol -> figure exits 0", rc == 0)
check("CLI: writes SVG", (OUT / "cli_tp53.svg").stat().st_size > 4000)
cli = json.loads((OUT / "cli.json").read_text())
check("CLI: JSON dump has the right transcript",
      cli["transcripts"][0]["name"] == "NM_000546.6" and cli["track"] == "ncbiRefSeqSelect")
check("CLI: FASTA header is 1-based with strand",
      (OUT / "cli.fa").read_text().startswith(
          ">TP53 hg38 chr17:7668421-7687490 strand=-"))
check("CLI: sequence link points at the UCSC API",
      cli["links"]["sequence"].startswith("https://api.genome.ucsc.edu/getData/sequence?"))

with contextlib.redirect_stderr(io.StringIO()):
    rc_bad = ugc.main(["--region", "not-a-region", "-o", str(OUT / "x.svg")])
check("CLI: bad --region exits 1 without a traceback", rc_bad == 1)

buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    ugc.main(["--list-tracks"])
check("CLI: --list-tracks shows only transcript tracks",
      buf.getvalue().split() == ["ncbiRefSeqSelect"], buf.getvalue().strip())

# --- 7. notebook API ------------------------------------------------------- #
g = ugc.draw_gene("TP53", cds_color="#7A5195", exon_shape="round")
check("draw_gene returns a GeneFigure", isinstance(g, ugc.GeneFigure), repr(g))
svg = g._repr_svg_()
check("GeneFigure renders inline SVG",
      svg.lstrip().startswith("<?xml") and "</svg>" in svg and "<text" in svg,
      f"{len(svg)} chars")
check("style kwargs reach the drawing", "#7A5195" in svg.lower() or "#7a5195" in svg.lower())
row = g.table()[0]
check("GeneFigure.table summarises the model",
      row["transcript"] == "NM_000546.6" and row["exons"] == 11
      and row["cds_bp"] == 1182 and row["start"] == 7668421, str(row))
check("GeneFigure exposes UCSC links",
      g.sequence_url.startswith("https://api.genome.ucsc.edu/getData/sequence?")
      and "hgTracks" in g.browser_url)
check("GeneFigure.sequence() fetches DNA", g.sequence() == "ACGT" * 100)
check("GeneFigure.save writes files",
      [Path(p).stat().st_size > 4000
       for p in g.save(str(OUT / "nb.svg"), str(OUT / "nb.png"), quiet=True)]
      == [True, True])
try:
    ugc.draw_gene("TP53", cds_colour="red")
    check("unknown style option raises", False)
except TypeError as e:
    check("unknown style option raises TypeError", "cds_colour" in str(e))
check("draw_gene accepts a Style object",
      ugc.draw_gene("TP53", style=ugc.Style(intron_mode="linear")).cartoon.st.intron_mode
      == "linear")
check("draw_gene accepts region= instead of a symbol",
      ugc.draw_gene(region="chr17:7,668,421-7,687,490").transcripts[0].name
      == "NM_000546.6")

# --- 8. cDNA (HGVS) -> genomic mapping ------------------------------------- #
mt, mc = ugc.CDNAMapper(t), ugc.CDNAMapper(c)

# Published hg38 coordinates for well-characterised variants.  The two marked
# (ref) were confirmed against the UCSC reference base at that position:
# chr7:117,642,566 is G (matches c.3846G>A) and chr17:7,674,872 is T, the
# complement of the A in c.659A>G on the minus strand.
KNOWN = [
    ("TP53 R175H",  mt, "c.524G>A",   7675088, 175, 5),
    ("TP53 Y220C",  mt, "c.659A>G",   7674872, 220, 6),      # (ref)
    ("TP53 G245S",  mt, "c.733G>A",   7674230, 245, 7),
    ("TP53 R248Q",  mt, "c.743G>A",   7674220, 248, 7),
    ("TP53 R273H",  mt, "c.818G>A",   7673802, 273, 8),
    ("TP53 R282W",  mt, "c.844C>T",   7673776, 282, 8),
    ("CFTR R117H",  mc, "c.350G>A",   117530975, 117, 4),
    ("CFTR G551D",  mc, "c.1652G>A",  117587806, 551, 12),
    ("CFTR W1282X", mc, "c.3846G>A",  117642566, 1282, 23),  # (ref)
]
for name, mp, hgvs, want_g, want_p, want_ex in KNOWN:
    got = mp.to_genomic(hgvs)
    check(f"{name} {hgvs} maps to the published coordinate",
          got == want_g and mp.protein_position(hgvs) == want_p
          and mp.exon_of(hgvs) == want_ex,
          f"got {got} / codon {mp.protein_position(hgvs)} / exon {mp.exon_of(hgvs)}")

check("TP53 c.1 is the ATG at chr17:7,676,594", mt.to_genomic("c.1") == 7676594)
check("CFTR c.1 is the ATG at chr7:117,480,095", mc.to_genomic("c.1") == 117480095)
check("c.-1 is one base 5' of the ATG (minus strand -> higher coordinate)",
      mt.to_genomic("c.-1") == 7676595)
check("CFTR c.-1 sits one base lower on the plus strand",
      mc.to_genomic("c.-1") == 117480094)
check("last CDS base maps to the genomic CDS end",
      mc.to_genomic("c.4443") == c.cds_end and mt.to_genomic("c.1182") == t.cds_start + 1)
check("c.* numbering starts after the stop codon",
      mc.to_genomic("c.*1") == c.cds_end + 1 and mt.to_genomic("c.*1") == t.cds_start)
check("intronic offsets step off the exon edge in transcript direction",
      mt.to_genomic("c.375+1") == mt.to_genomic("c.375") - 1
      and mc.to_genomic("c.1584+1") == mc.to_genomic("c.1584") + 1)
check("intronic positions report no exon", mt.exon_of("c.375+1") is None)
check("round trip: genomic -> tx -> genomic",
      all(mc.tx_to_genomic(mc.genomic_to_tx(g)) == g
          for g in range(117480025, 117480148)))
check("every coding base maps back to a codon",
      mc.protein_position("c.4443") == 1481 and mc.protein_position("c.3") == 1)
check("non-coding transcript uses n. numbering from the 5' end",
      ugc.CDNAMapper(x).to_genomic("c.1") == x.tx_end)

# parsing
p = ugc.parse_hgvs("NM_000546.6:c.376-2A>G")
check("parser handles prefix, offset and ref>alt",
      (p.base, p.offset, p.utr3) == (376, -2, False), str(p))
check("parser handles c.*38", ugc.parse_hgvs("c.*38G>A").utr3)
check("parser handles c.-124", ugc.parse_hgvs("c.-124C>T").base == -124)
check("parser handles a bare integer", ugc.parse_hgvs("743").base == 743)
check("parser takes the first position of a range",
      ugc.parse_hgvs("c.1521_1523del").base == 1521)
for bad, why in [("p.Arg248Gln", "protein"), ("", "empty"), ("rubbish", "junk")]:
    try:
        ugc.parse_hgvs(bad)
        check(f"parser rejects {why}", False)
    except ugc.CDNAError as e:
        check(f"parser rejects {why} with a useful message",
              len(str(e)) > 15, str(e)[:70])
try:
    mt.to_genomic("c.99999")
    check("position past the 3' end is rejected", False)
except ugc.CDNAError as e:
    check("position past the 3' end is rejected", "3' end" in str(e), str(e)[:60])
try:
    mt.to_genomic("c.-99999")
    check("promoter position upstream of the transcript is rejected", False)
except ugc.CDNAError as e:
    check("promoter position gets an explanation, not just a rejection",
          "upstream" in str(e) and "5'UTR" in str(e), str(e)[:70])

# placement + rendering
vlist = [ugc.Variant(label=l, cdna=cd, category=k, count=n) for l, cd, k, n in [
    ("R175H", "c.524G>A", "Pathogenic", 9), ("Y220C", "c.659A>G", "Pathogenic", 5),
    ("G245S", "c.733G>A", "Pathogenic", 4), ("R248Q", "c.743G>A", "Pathogenic", 12),
    ("R248W", "c.742C>T", "Pathogenic", 7), ("R273H", "c.818G>A", "Pathogenic", 8),
    ("R273C", "c.817C>T", "Pathogenic", 4), ("R282W", "c.844C>T", "Pathogenic", 3),
    ("P72R", "c.215C>G", "Benign", 1), ("splice", "c.376-2A>G", "VUS", 1),
    ("5'UTR", "c.-28C>T", "VUS", 1), ("3'UTR", "c.*38G>A", "VUS", 2),
    ("bad", "p.Arg248Gln", "", 1)]]
placed, failed = ugc.place_variants(vlist, t)
check("place_variants places the valid rows and isolates the bad one",
      len(placed) == 12 and len(failed) == 1 and failed[0].label == "bad")
check("failed variant explains itself", "protein" in failed[0].error.lower())

vc = ugc.GeneCartoon([t], ugc.Style(figure_width_in=10), "TP53", "hg38",
                     "ncbiRefSeqSelect", links=links, variants=placed)
vc.render([str(OUT / "tp53_variants.svg"), str(OUT / "tp53_variants.png")],
          quiet=True)
vbody = (OUT / "tp53_variants.svg").read_text()
check("variant figure renders", (OUT / "tp53_variants.svg").stat().st_size > 8000)
check("all variant labels appear in the SVG",
      all(f">{v.label}<" in vbody for v in placed if v.label.isalnum()),
      )
check("categories get distinct colours",
      len({vc.variant_color(v) for v in placed}) == 3)
genomic_order = [v.label for v, _ in
                 sorted(((v, vc.map.x(v.genomic - 1)) for v in placed),
                        key=lambda p: p[1])]

# default is the stacked style: no stems, markers on their true position
layout = vc._variant_layout()
check("stacked is the default variant style", ugc.Style().variant_style == "stacked")
check("stacked keeps every marker on its true position",
      all(abs(anchor - head) < 1e-12 for _, anchor, head, _ in layout))
check("stacked piles crowded variants into tiers",
      max(tier for _, _, _, tier in layout) > 0,
      f"{max(tier for _, _, _, tier in layout) + 1} tiers")
check("stacked reuses the bottom tier for well-separated variants",
      sum(1 for _, _, _, tier in layout if tier == 0) >= 5,
      f"{sum(1 for _, _, _, tier in layout if tier == 0)} on tier 0")
check("stacked keeps variants in genomic order",
      [v.label for v, _, _, _ in layout] == genomic_order)
# NB: matplotlib emits Line2D as <path>, not <line> -- counting tags would
# pass vacuously.  Count strokes in the stem colour instead.
STEM = ugc.Style().variant_stem_color.lower().lstrip("#")
stacked_svg = ugc.GeneCartoon([t], ugc.Style(variant_style="stacked"), "TP53",
                              "hg38", "x", variants=placed).to_svg().lower()
lolli_svg = ugc.GeneCartoon([t], ugc.Style(variant_style="lollipop"), "TP53",
                            "hg38", "x", variants=placed).to_svg().lower()
check("stacked draws no stems at all",
      stacked_svg.count(STEM) == 0,
      f"{stacked_svg.count(STEM)} stem-coloured strokes")
check("lollipop does draw stems (so the check above is meaningful)",
      lolli_svg.count(STEM) >= len(placed),
      f"{lolli_svg.count(STEM)} stem-coloured strokes for {len(placed)} variants")
check("the two styles produce genuinely different output",
      stacked_svg != lolli_svg)

# lollipop remains available
lolli = ugc.GeneCartoon([t], ugc.Style(variant_style="lollipop"), "TP53",
                        "hg38", "x", variants=placed)
lolli.figure()
llayout = lolli._variant_layout()
check("lollipop spread separates every head",
      all(b[2] - a[2] > 1e-6 for a, b in zip(llayout, llayout[1:])))
check("lollipop keeps heads in genomic order",
      [v.label for v, _, _, _ in llayout] == genomic_order)
check("lollipop stack mode tiers colliding heads",
      max(tier for _, _, _, tier in ugc.GeneCartoon(
          [t], ugc.Style(variant_style="lollipop", variant_collision="stack"),
          "TP53", "hg38", "x", variants=placed)._variant_layout()) > 0)

# --- 9. variant file loading ----------------------------------------------- #
try:
    import pandas as pd
    import variants as varmod

    def _csv(df, sep=","):
        p = OUT / "vars.csv"
        df.to_csv(p, index=False, sep=sep)
        return p

    vf = varmod.load_variant_table(_csv(pd.DataFrame({
        "Sample ID": ["S1", "S2", "S3"],
        "HGVSc": ["NM_000546.6:c.743G>A", "c.524G>A", "c.818G>A"],
        "HGVSp": ["R248Q", "R175H", "R273H"],
        "Clinical Significance": ["Pathogenic", "Pathogenic", "VUS"],
        "No. of cases": [12, 9, 8]})))
    check("loader detects messy real-world headers",
          (vf.columns["cdna"], vf.columns["label"], vf.columns["category"],
           vf.columns["count"])
          == ("HGVSc", "HGVSp", "Clinical Significance", "No. of cases"),
          str({k: v for k, v in vf.columns.items() if v}))
    check("loader reads counts", [v.count for v in vf.variants] == [12, 9, 8])

    vf1 = varmod.load_variant_table(_csv(pd.DataFrame(
        {"cDNA": ["c.743G>A", "c.524G>A", "c.376-2A>G"]})))
    check("single-column CSV is not shredded by delimiter sniffing",
          [v.cdna for v in vf1.variants]
          == ["c.743G>A", "c.524G>A", "c.376-2A>G"],
          str([v.cdna for v in vf1.variants]))

    vf2 = varmod.load_variant_table(_csv(pd.DataFrame(
        {"id": [1, 2], "tissue": ["lung", "breast"]})))
    check("a numeric ID column is not mistaken for cDNA positions",
          not vf2.variants and vf2.warnings)

    vf3 = varmod.load_variant_table(_csv(pd.DataFrame(
        {"col_a": ["x", "y"], "col_b": ["c.743G>A", "c.524G>A"]})))
    check("unlabelled cDNA column is found by content", vf3.columns["cdna"] == "col_b")

    vf4 = varmod.load_variant_table(_csv(pd.DataFrame(
        {"Variant": ["c.743G>A"], "hg38": ["chr17:7,674,220"]}), sep=";"))
    check("semicolon files and chr:pos genomic columns both work",
          vf4.variants and vf4.variants[0].genomic == 7674220)

    xl = OUT / "vars.xlsx"
    varmod.write_template(xl)
    vf5 = varmod.load_variant_table(xl)
    check("Excel round-trips through the template writer",
          len(vf5.variants) == 8 and vf5.columns["cdna"] == "cDNA")
    pl5, fa5 = ugc.place_variants(vf5.variants, t)
    check("every template variant places on TP53", len(pl5) == 8 and not fa5)

    # header styles seen in real curated variant databases
    for headers, want in [
        (["HGVS_cDNA", "Location", "Domain", "Mutation_Type"],
         {"cdna": "HGVS_cDNA", "category": "Mutation_Type", "domain": "Domain"}),
        (["cDNA", "Variant type"], {"cdna": "cDNA", "category": "Variant type"}),
        (["Mutation CDS", "Mutation AA", "Mutation Description", "Count"],
         {"cdna": "Mutation CDS", "label": "Mutation AA", "count": "Count"}),
    ]:
        d = pd.DataFrame({h: (["c.743G>A", "c.524G>A"] if "cdna" in h.lower()
                              or "cds" in h.lower() else ["x", "y"])
                          for h in headers})
        got = varmod.detect_columns(d)
        check(f"detects {headers[0]!r}-style headers",
              all(got.get(k) == v for k, v in want.items()),
              str({k: v for k, v in got.items() if v}))

    # NB: 'Location' here holds cDNA positions, not genomic ones -- claiming it
    # as the genomic column would silently place every variant in the wrong place
    d = pd.DataFrame({"HGVS_cDNA": ["c.743G>A"], "Location": [743]})
    check("a cDNA 'Location' column is not mistaken for a genomic coordinate",
          varmod.detect_columns(d).get("genomic") is None)

except ImportError:                                    # pragma: no cover
    print("SKIP  variant file tests (pandas not available)")


# --- 10. category filtering and per-category colour ------------------------ #
cat_vars = [ugc.Variant(label=f"v{i}", cdna=cd, category=cat)
            for i, (cd, cat) in enumerate([
                ("c.524G>A", "Missense"), ("c.743G>A", "Missense"),
                ("c.818G>A", "Nonsense"), ("c.844C>T", "Nonsense"),
                ("c.376-2A>G", "Splice site change")])]
cat_placed, _ = ugc.place_variants(cat_vars, t)

subset = [v for v in cat_placed if v.category != "Nonsense"]
cs = ugc.GeneCartoon([t], ugc.Style(), "TP53", "hg38", "x", variants=subset)
check("filtering a category removes it from the legend too",
      "Nonsense" not in cs._variant_colors and len(cs._variant_colors) == 2,
      str(sorted(cs._variant_colors)))

for v in cat_placed:
    v.color = {"Missense": "#111111", "Nonsense": "#222222",
               "Splice site change": "#333333"}[v.category]
cc = ugc.GeneCartoon([t], ugc.Style(), "TP53", "hg38", "x", variants=cat_placed)
check("an explicit colour drives the markers",
      {cc.variant_color(v) for v in cat_placed}
      == {"#111111", "#222222", "#333333"})
check("...and the legend agrees with the markers",
      all(cc._variant_colors[v.category] == cc.variant_color(v)
          for v in cat_placed),
      str(cc._variant_colors))

# --- 10b. lanes: one labelled row per category ----------------------------- #
lane_style = ugc.Style(variant_style="lanes", figure_width_in=11)
ln = ugc.GeneCartoon([t], lane_style, "TP53", "hg38", "x", variants=cat_placed)
lane_svg = ln.to_svg()
check("lanes gives every category its own row",
      [c for c, _, _ in ln._lane_rows]
      == list(ln._variant_colors),
      str([(c, n) for c, _, n in ln._lane_rows]))
check("lanes keeps each variant inside its own lane",
      all(any(base <= tier < base + n
              for cat, base, n in ln._lane_rows
              if cat == v.category)
          for v, _, _, tier in ln._variant_layout()))
check("lanes writes the category name on the left",
      all(f">{c}<" in lane_svg for c in ln._variant_colors))
check("lanes suppresses the now-redundant colour key",
      lane_svg.count(">Missense<") == 1, "category named once, not twice")
check("lanes widens the left margin for long names",
      ugc.GeneCartoon([t], ugc.Style(variant_style="lanes"), "TP53", "hg38",
                      "x", variants=[ugc.Variant(
                          label="v", cdna="c.1", genomic=7676594,
                          category="Large structural change (>50 bp)")]
                      ).figure().axes[0].get_xlim()[0] < -0.2)
check("lanes and stacked differ", lane_svg != ugc.GeneCartoon(
    [t], ugc.Style(figure_width_in=11), "TP53", "hg38", "x",
    variants=cat_placed).to_svg())

# --- 10bb. focusing on an exon --------------------------------------------- #
# Exon numbering is biological, so on the minus-strand TP53 exon 1 must be the
# *highest* genomic block and exon 11 the lowest -- the classic thing to get
# backwards.
check("exon_span follows biological numbering on the minus strand",
      t.exon_span(1) == t.exons[-1] and t.exon_span(11) == t.exons[0],
      f"exon1={t.exon_span(1)} exon11={t.exon_span(11)}")
check("exon_span follows genomic order on the plus strand",
      c.exon_span(1) == c.exons[0] and c.exon_span(27) == c.exons[-1])
check("exon_span rejects numbers outside the transcript",
      t.exon_span(0) is None and t.exon_span(99) is None)

w1 = t.exon_window(5, flank=0)
check("a single-exon window is exactly that exon", w1 == t.exon_span(5))
check("flank widens the window symmetrically",
      t.exon_window(5, flank=250) == (w1[0] - 250, w1[1] + 250))
w2 = t.exon_window(4, 6, flank=0)
check("an exon range spans all of them",
      w2[0] == min(t.exon_span(n)[0] for n in (4, 5, 6))
      and w2[1] == max(t.exon_span(n)[1] for n in (4, 5, 6)))
check("a range reads the same either way round",
      t.exon_window(6, 4) == t.exon_window(4, 6))
try:
    t.exon_window(99)
    check("an impossible exon is refused", False)
except ValueError as e:
    check("an impossible exon is refused with a helpful message",
          "1-11" in str(e), str(e)[:60])

for spec, want in [(5, "exon 5"), ("5", "exon 5"), ("exon 5", "exon 5"),
                   ("4-6", "exons 4–6"), ((4, 6), "exons 4–6"),
                   (None, ""), ("whole gene", "")]:
    b, lbl = ugc.resolve_focus(t, spec, flank=0)
    check(f"resolve_focus({spec!r}) -> {want or 'whole gene'}",
          lbl == want and (b is None) == (want == ""))

foc = ugc.GeneCartoon([t], ugc.Style(), "TP53", "hg38", "x",
                      variants=cat_placed, focus=t.exon_window(7, flank=200),
                      focus_label="exon 7")
check("focus keeps only the variants inside the window",
      len(foc.variants) + len(foc.offscreen_variants) == len(cat_placed)
      and all(foc.map.g_start < v.genomic <= foc.map.g_end
              for v in foc.variants),
      f"{len(foc.variants)} in, {len(foc.offscreen_variants)} out")
foc_svg = foc.to_svg()
check("focus names the region in the subtitle",
      "exon 7" in foc_svg and "exon 7" not in ugc.GeneCartoon(
          [t], ugc.Style(), "TP53", "hg38", "x").to_svg())
check("focus numbers only the exons on screen",
      foc_svg.count(">7<") >= 1 and ">1<" not in foc_svg,
      "distant exon numbers are not drawn")
check("the whole-gene view still shows every exon number",
      all(f">{n}<" in ugc.GeneCartoon([t], ugc.Style(), "TP53", "hg38",
                                      "x").to_svg() for n in range(1, 12)))

part = ugc.CoordinateMapper([t], ugc.Style(), bounds=t.exon_window(7, flank=100))
check("a bounded mapper spans only its window",
      (part.g_start, part.g_end) == t.exon_window(7, flank=100))
check("clip trims a feature that straddles the edge",
      part.clip(part.g_start - 500, part.g_start + 50)
      == (part.g_start, part.g_start + 50))
check("clip drops a feature entirely outside",
      part.clip(part.g_start - 900, part.g_start - 400) is None)
check("visible() agrees with clip()",
      all((part.clip(a, b) is not None) == part.visible(a, b)
          for a, b in [(part.g_start - 9, part.g_start - 1),
                       (part.g_start - 9, part.g_start + 9),
                       (part.g_start + 5, part.g_end - 5),
                       (part.g_end + 1, part.g_end + 9)]))
check("a focused figure is narrower in bp than the whole gene",
      (part.g_end - part.g_start) < (t.tx_end - t.tx_start))

# --- 10bc. density: shrinking markers must shrink the figure --------------- #
# A dense synthetic cohort: enough collisions to force many rows.
dense = []
for i in range(160):
    dense.append(ugc.Variant(label="", cdna=f"c.{200 + (i % 40)}",
                             category=["A", "B", "C"][i % 3]))
dense, _ = ugc.place_variants(dense, t)


def fig_height(**kw):
    st_ = ugc.Style(variant_style="lanes", figure_width_in=11, **kw)
    cart = ugc.GeneCartoon([t], st_, "TP53", "hg38", "x", variants=dense)
    f = cart.figure()
    h = f.get_size_inches()[1]
    import matplotlib.pyplot as _plt
    _plt.close(f)
    return h, max(tr for _, _, _, tr in cart._variant_layout()) + 1


h_big, rows_big = fig_height(variant_head_size=42)
h_small, rows_small = fig_height(variant_head_size=12)
check("shrinking the markers shortens the figure",
      h_small < h_big * 0.9,
      f"{h_big:.1f}in -> {h_small:.1f}in")
check("row pitch tracks marker size",
      ugc.GeneCartoon([t], ugc.Style(variant_head_size=12), "TP53", "hg38",
                      "x", variants=dense)._stack_gap()
      < ugc.GeneCartoon([t], ugc.Style(variant_head_size=42), "TP53", "hg38",
                        "x", variants=dense)._stack_gap())
h_tight, _ = fig_height(variant_head_size=12, variant_row_spacing=1.0)
check("row spacing multiplier shortens it further",
      h_tight < h_small, f"{h_small:.1f}in -> {h_tight:.1f}in")
_, rows_packed = fig_height(variant_head_size=12, variant_pack=0.7)
check("tighter horizontal packing needs fewer rows",
      rows_packed < rows_small, f"{rows_small} -> {rows_packed} rows")
check("an explicit stack gap still overrides the automatic one",
      abs(ugc.GeneCartoon([t], ugc.Style(variant_stack_gap=0.5), "TP53",
                          "hg38", "x", variants=dense)._stack_gap() - 0.5)
      < 1e-9)
check("no variant is dropped when compacting",
      all(len(ugc.GeneCartoon([t], ugc.Style(variant_style="lanes", **kw),
                              "TP53", "hg38", "x", variants=dense).variants)
          == len(dense)
          for kw in ({}, {"variant_head_size": 8},
                     {"variant_pack": 0.6, "variant_row_spacing": 0.9})),
      "compacting is presentation only")

# --- 10c. SVG ids for interactivity ---------------------------------------- #
gid_svg = ugc.GeneCartoon([t], ugc.Style(), "TP53", "hg38", "x",
                          variants=cat_placed)
svg_txt = gid_svg.to_svg()
gids = re.findall(r'id="(variant-\d+)"', svg_txt)
tips = gid_svg.variant_tooltips()
check("every drawn variant gets an id in the SVG",
      len(gids) == len(cat_placed) and len(set(gids)) == len(gids),
      f"{len(gids)} ids for {len(cat_placed)} variants")
check("ids line up with the tooltip records",
      {t_["gid"] for t_ in tips} == set(gids))
check("tooltip records carry what a reader needs",
      all(t_["label"] and t_["position"] and t_["color"] for t_ in tips))
check("tooltip colour matches the drawn marker",
      all(t_["color"] == gid_svg.variant_color(v)
          for t_, v in zip(tips, gid_svg.variants)))
check("ids are stable across re-renders",
      re.findall(r'id="(variant-\d+)"', gid_svg.to_svg()) == gids)

many = [ugc.Variant(label=f"m{i}", cdna="c.524G>A",
                    category=f"a rather long category name {i}")
        for i in range(8)]
mp, _ = ugc.place_variants(many, t)
narrow = ugc.GeneCartoon([t], ugc.Style(figure_width_in=6.0), "TP53", "hg38",
                         "x", variants=mp)
narrow.figure()
items, lines = narrow._variant_legend_layout()
check("a long category legend wraps instead of running off the page",
      lines > 1 and all(x < 1.0 for _, x, _ in items),
      f"{lines} rows, max x {max(x for _, x, _ in items):.2f}")
check("wrapped legend rows are reserved in the layout",
      narrow.figure().get_size_inches()[1]
      > ugc.GeneCartoon([t], ugc.Style(figure_width_in=6.0), "TP53", "hg38",
                        "x", variants=mp[:1]).figure().get_size_inches()[1])
# --- 11. every control, across its full range ------------------------------ #
# Renders ~100 configurations and inspects the drawn result for overlapping
# text, artists escaping the canvas and markers colliding with the gene.
try:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import sweep as _sweep
    n_cfg, bad = _sweep.sweep()
    check(f"all {n_cfg} control settings render cleanly", not bad,
          "; ".join(f"{name}: {p[0]}" for (_, name), p in list(bad.items())[:3])
          or f"{n_cfg} configurations")
except ImportError:                                    # pragma: no cover
    print("SKIP  control sweep (sweep.py not importable)")

print()
print(f"{'ALL TESTS PASSED' if not failures else str(len(failures)) + ' FAILURE(S): ' + ', '.join(failures)}")
sys.exit(1 if failures else 0)
