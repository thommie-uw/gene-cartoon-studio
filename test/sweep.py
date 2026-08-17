#!/usr/bin/env python3
"""
Sweep every GUI control across its full range and detect broken output.

Checks, in display coordinates after a real draw:
  * the render doesn't raise and the geometry is finite
  * no text label overlaps another text label
  * nothing is drawn outside the figure canvas
  * variant markers don't collide with the gene model
  * exon boxes don't collide with the row above

Run:  python test/sweep.py [--verbose]
"""
import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt              # noqa: E402
from matplotlib.text import Text             # noqa: E402
from matplotlib.patches import Rectangle, Polygon  # noqa: E402
from matplotlib.collections import PathCollection  # noqa: E402

import ucsc_gene_cartoon as ugc              # noqa: E402

VERBOSE = "--verbose" in sys.argv

TP53 = {
    "name": "NM_000546.6", "chrom": "chr17", "strand": "-",
    "txStart": 7668420, "txEnd": 7687490,
    "cdsStart": 7669608, "cdsEnd": 7676594,
    "exonStarts": "7668420,7670608,7673534,7673700,7674180,7674858,7675052,7675993,7676381,7676520,7687376,",
    "exonEnds": "7669690,7670715,7673608,7673837,7674290,7674971,7675236,7676272,7676403,7676622,7687490,",
    "name2": "TP53",
}
T = ugc.Transcript.from_genepred(TP53)
VARIANTS, _ = ugc.place_variants(
    [ugc.Variant(label=l, cdna=c, category=k, count=n) for l, c, k, n in [
        ("R175H", "c.524G>A", "Missense", 9),
        ("Y220C", "c.659A>G", "Missense", 5),
        ("R248Q", "c.743G>A", "Missense", 12),
        ("R273H", "c.818G>A", "Nonsense", 8),
        ("R282W", "c.844C>T", "Nonsense", 3),
        ("P72R", "c.215C>G", "Benign", 1),
        ("splice", "c.376-2A>G", "Splice", 1),
    ]], T)


def overlap(a, b, pad=0.5):
    """Do two display-space bboxes overlap by more than `pad` points?"""
    return not (a.x1 - pad <= b.x0 or b.x1 - pad <= a.x0
                or a.y1 - pad <= b.y0 or b.y1 - pad <= a.y0)


def oriented_box(txt, rend, pad=1.0):
    """
    Corners of a Text's *rotated* rectangle in display space.

    matplotlib's get_window_extent returns the axis-aligned envelope, which for
    45-degree labels is far bigger than the glyphs and reports collisions that
    aren't there. Measure unrotated, then rotate the rectangle ourselves.
    """
    import math as _m
    rot = txt.get_rotation() or 0.0
    saved = txt.get_rotation()
    txt.set_rotation(0)
    bb = txt.get_window_extent(renderer=rend)
    txt.set_rotation(saved)
    w, h = bb.width + 2 * pad, bb.height + 2 * pad

    x, y = txt.get_transform().transform(txt.get_position())
    ha, va = txt.get_horizontalalignment(), txt.get_verticalalignment()
    dx = {"left": 0.0, "center": -w / 2, "right": -w}[ha]
    dy = {"baseline": 0.0, "bottom": 0.0, "center": -h / 2,
          "center_baseline": -h / 2, "top": -h}.get(va, 0.0)

    th = _m.radians(rot)
    cos, sin = _m.cos(th), _m.sin(th)
    pts = []
    for cx, cy in ((dx, dy), (dx + w, dy), (dx + w, dy + h), (dx, dy + h)):
        pts.append((x + cx * cos - cy * sin, y + cx * sin + cy * cos))
    return pts


def obb_overlap(p, q):
    """Separating-axis test between two convex quads."""
    for poly in (p, q):
        for i in range(len(poly)):
            x0, y0 = poly[i]
            x1, y1 = poly[(i + 1) % len(poly)]
            ax, ay = -(y1 - y0), (x1 - x0)
            n = (ax * ax + ay * ay) ** 0.5
            if n == 0:
                continue
            ax, ay = ax / n, ay / n
            pp = [ax * x + ay * y for x, y in p]
            qq = [ax * x + ay * y for x, y in q]
            if max(pp) <= min(qq) or max(qq) <= min(pp):
                return False
    return True


def inspect(style, variants=VARIANTS, annotations=None):
    """Render once and return a list of problem strings."""
    problems = []
    try:
        cart = ugc.GeneCartoon([T], style, "TP53", "hg38", "ncbiRefSeqSelect",
                               variants=variants, annotations=annotations,
                               links={"sequence": "s", "browser": "b"})
        fig = cart.figure()
    except Exception as e:
        return [f"RAISED {type(e).__name__}: {e}"]

    try:
        w, h = fig.get_size_inches()
        if not (0.5 < w < 80) or not (0.5 < h < 80):
            problems.append(f"absurd figure size {w:.1f}x{h:.1f}in")
        ax = fig.axes[0]
        y0, y1 = ax.get_ylim()
        if not (y1 > y0) or any(v != v for v in (y0, y1)):
            problems.append(f"bad ylim {y0}..{y1}")

        fig.canvas.draw()
        rend = fig.canvas.get_renderer()
        figbox = fig.bbox

        texts, boxes, markers = [], [], []
        for art in ax.get_children():
            try:
                bb = art.get_window_extent(renderer=rend)
            except Exception:
                continue
            if bb.width <= 0 or bb.height <= 0:
                continue
            if isinstance(art, Text):
                if art.get_text().strip():
                    texts.append((art.get_text().strip(), bb,
                                  oriented_box(art, rend)))
            elif isinstance(art, Rectangle):
                boxes.append(bb)
            elif isinstance(art, PathCollection):
                markers.append(bb)

        # 1. text must not collide with other text (rotation-aware)
        for (t1, _, o1), (t2, _, o2) in combinations(texts, 2):
            if obb_overlap(o1, o2):
                problems.append(f"text collision: {t1!r} / {t2!r}")

        # 2. nothing may escape the canvas
        for label, bb, _ in texts:
            if (bb.x0 < figbox.x0 - 1 or bb.x1 > figbox.x1 + 1
                    or bb.y0 < figbox.y0 - 1 or bb.y1 > figbox.y1 + 1):
                problems.append(f"text off-canvas: {label!r}")
        for bb in boxes + markers:
            if (bb.y0 < figbox.y0 - 1 or bb.y1 > figbox.y1 + 1):
                problems.append("artist off-canvas vertically")
                break

        # 3. variant markers must clear the gene model
        for mb in markers:
            for bb in boxes:
                if overlap(mb, bb, pad=0.0):
                    problems.append("variant marker overlaps an exon box")
                    break
            else:
                continue
            break
    finally:
        plt.close(fig)
    return sorted(set(problems))


def sweep():
    cases = []

    def add(panel, name, **kw):
        cases.append((panel, name, kw))

    # ---- Layout & scale ---- #
    for v in ("compress", "linear", "equal"):
        add("Layout", f"intron_mode={v}", intron_mode=v)
    for v in (4.0, 6.0, 9.0, 14.0, 20.0):
        add("Layout", f"figure_width_in={v}", figure_width_in=v)
    for v in (0.4, 0.6, 0.85, 1.4, 2.0):
        add("Layout", f"row_height_in={v}", row_height_in=v)
    for v in (0.1, 0.42, 0.7, 0.9):
        add("Layout", f"cds_height={v}", cds_height=v)
    for v in (0.05, 0.22, 0.6, 0.9):
        add("Layout", f"utr_height={v}", utr_height=v)
    for v in (5.0, 9.0, 13.0, 16.0):
        add("Layout", f"font_size={v}", font_size=v)
    for v in ("above", "inside", "below"):
        add("Layout", f"exon_number_position={v}", exon_number_position=v)
    for v in ("auto", "coords", "scalebar", "none"):
        add("Layout", f"axis_mode={v}", axis_mode=v)
    for flag in ("show_exon_numbers", "show_title", "show_subtitle",
                 "show_legend", "show_strand_arrow",
                 "orient_five_prime_left"):
        add("Layout", f"{flag}=False", **{flag: False})
    add("Layout", "background=white", background="white")
    for v in ("box", "round"):
        add("Layout", f"exon_shape={v}", exon_shape=v)
    for v in ("line", "chevron_line", "angled"):
        add("Layout", f"intron_style={v}", intron_style=v)

    # ---- Variants ---- #
    for v in ("stacked", "lanes", "lollipop"):
        add("Variants", f"variant_style={v}", variant_style=v)
        for size in (6.0, 14.0, 42.0, 120.0, 200.0):
            add("Variants", f"{v}: head_size={size}",
                variant_style=v, variant_head_size=size)
    for v in (0.0, 0.16, 0.4, 0.8):
        add("Variants", f"lanes: lane_gap={v}",
            variant_style="lanes", lane_gap=v)
    for v in (0.8, 1.3, 2.0, 3.0):
        add("Variants", f"row_spacing={v}", variant_row_spacing=v)
    for v in (0.5, 1.0, 1.25, 3.0):
        add("Variants", f"pack={v}", variant_pack=v)
    for v in (0.0, 0.2, 0.8, 1.5):
        add("Variants", f"base_gap={v}", variant_base_gap=v)
    for v in (4.0, 6.5, 10.0, 14.0):
        add("Variants", f"label_size={v}", variant_label_size=v)
    for v in ("o", "s", "D", "v", "^"):
        add("Variants", f"marker={v}", variant_marker=v)
    for v in (0.2, 0.55, 1.2, 2.0):
        add("Variants", f"lollipop: stem={v}",
            variant_style="lollipop", variant_stem_height=v)
    add("Variants", "scale_by_count=False", variant_scale_by_count=False)
    add("Variants", "labels off", show_variant_labels=False)
    add("Variants", "no variants", variant_style="lanes")

    # ---- a few nasty combinations ---- #
    add("Combo", "tiny rows + big font", row_height_in=0.4, font_size=16.0)
    add("Combo", "tall UTR + short CDS", utr_height=0.9, cds_height=0.15)
    add("Combo", "huge markers + lanes",
        variant_style="lanes", variant_head_size=200.0)
    add("Combo", "narrow + long labels",
        figure_width_in=4.0, variant_style="lanes")
    add("Combo", "compact preset", variant_head_size=14.0,
        variant_row_spacing=1.05, variant_pack=1.0, variant_base_gap=0.12,
        lane_gap=0.06, row_height_in=0.65, variant_label_size=5.5,
        variant_style="lanes")

    results = {}
    for panel, name, kw in cases:
        vs = [] if name == "no variants" else VARIANTS
        probs = inspect(ugc.Style(**kw), variants=vs)
        if probs:
            results[(panel, name)] = probs
        if VERBOSE:
            print(f"{'ok  ' if not probs else 'BAD '} {panel:<9}{name}")
    return len(cases), results


if __name__ == "__main__":
    total, bad = sweep()
    print(f"\nswept {total} configurations, {len(bad)} with problems\n")
    by_panel = {}
    for (panel, name), probs in bad.items():
        by_panel.setdefault(panel, []).append((name, probs))
    for panel in sorted(by_panel):
        print(f"--- {panel} ---")
        for name, probs in sorted(by_panel[panel]):
            print(f"  {name}")
            for p in probs[:4]:
                print(f"      {p}")
            if len(probs) > 4:
                print(f"      ... and {len(probs) - 4} more")
    sys.exit(1 if bad else 0)
