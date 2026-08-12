#!/usr/bin/env python3
"""
Gene Cartoon Studio -- a point-and-click front end for ucsc_gene_cartoon.

Start it with:

    python run_gui.py                 (or)      streamlit run app.py

Type a gene symbol, adjust the controls, upload a variant spreadsheet,
download a figure.  No code required.
"""

from __future__ import annotations

import io
import json
import re
import traceback
from dataclasses import asdict, fields
from pathlib import Path
from typing import Any, Dict, List, Optional

import streamlit as st
import streamlit.components.v1 as components

import ucsc_gene_cartoon as ugc
from ucsc_gene_cartoon import (
    GENE_TRACKS, CDNAMapper, GeneCartoon, Style, Transcript, UCSCClient,
    UCSCError, Variant, place_variants, select_transcripts,
)

try:
    import variants as varmod
    import pandas as pd
    HAVE_TABLES = True
except ImportError:                                    # pragma: no cover
    HAVE_TABLES = False

#: this app needs at least this version of ucsc_gene_cartoon.py
ENGINE_MIN = (1, 5, 0)

ASSEMBLIES = ["hg38", "hg19", "mm39", "mm10", "rn7", "danRer11",
              "dm6", "ce11", "sacCer3", "galGal6", "susScr11", "bosTau9"]

#: one-click density settings. Neither hides a variant -- they only change
#: how much room each one is given.
COMPACT_PRESET = {
    "variant_head_size": 14.0, "variant_row_spacing": 1.05,
    "variant_pack": 1.0, "variant_base_gap": 0.12, "lane_gap": 0.06,
    "row_height_in": 0.65, "variant_label_size": 5.5,
}
ROOMY_PRESET = {
    "variant_head_size": 42.0, "variant_row_spacing": 1.30,
    "variant_pack": 1.25, "variant_base_gap": 0.20, "lane_gap": 0.16,
    "row_height_in": 0.85, "variant_label_size": 6.5,
}

FORMATS = [("SVG (vector, editable)", "svg", "image/svg+xml"),
           ("PDF (vector)", "pdf", "application/pdf"),
           ("PNG", "png", "image/png"),
           ("TIFF", "tif", "image/tiff")]


# --------------------------------------------------------------------------- #
#  Data fetching (cached so the sliders stay snappy)
# --------------------------------------------------------------------------- #

@st.cache_data(show_spinner=False, ttl=3600)
def fetch_models(gene: str, genome: str, track: str):
    """Resolve a symbol and pull transcript models. Returns (txs, links, track)."""
    client = UCSCClient(genome, quiet=True)
    chrom, start, end = client.locate_symbol(gene)

    if track and track != "auto":
        candidates = [track]
    else:
        available = set(client.gene_pred_tracks())
        candidates = [t for t in GENE_TRACKS if t in available] or ["ncbiRefSeq"]

    txs, used = [], candidates[0]
    for tr in candidates:
        found = client.transcripts_in_region(tr, chrom, start, end)
        found = [t for t in found if t.gene.upper() == gene.upper()] or found
        if found:
            txs, used = found, tr
            break
    if not txs:
        raise UCSCError(f"No transcript models found for {gene} in {genome}.")

    g0 = min(t.tx_start for t in txs)
    g1 = max(t.tx_end for t in txs)
    links = {"sequence": client.sequence_api_url(chrom, g0, g1),
             "browser": client.browser_url(chrom, g0, g1)}
    return txs, links, used


@st.cache_data(show_spinner=False, ttl=3600)
def available_tracks(genome: str) -> List[str]:
    try:
        return UCSCClient(genome, quiet=True).gene_pred_tracks()
    except UCSCError:
        return []


# --------------------------------------------------------------------------- #
#  Helpers
# --------------------------------------------------------------------------- #

def show_figure(cartoon: GeneCartoon, interactive: bool = True) -> None:
    """Preview the cartoon, scaled to the page."""
    if interactive and cartoon.variants:
        try:
            show_interactive(cartoon)
            return
        except Exception:                              # pragma: no cover
            st.caption("Interactive view unavailable; showing a static image.")
    png = cartoon.to_bytes("png")
    try:
        st.image(png, use_container_width=True)
    except TypeError:                                  # older Streamlit
        st.image(png, use_column_width=True)


def show_interactive(cartoon: GeneCartoon) -> None:
    """
    Render the figure as live SVG with hover tooltips and click-to-pin.

    This is the *same* SVG offered for download -- the markers simply carry
    ids, and a little JavaScript attaches behaviour to them. Nothing about the
    published figure changes.
    """
    svg = cartoon.to_svg()
    svg = svg[svg.index("<svg"):]

    m = re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', svg)
    aspect = (float(m.group(2)) / float(m.group(1))) if m else 0.4
    # make it fluid inside the iframe
    svg = re.sub(r'(<svg[^>]*?)\s+width="[\d.]+pt"\s+height="[\d.]+pt"',
                 r'\1 width="100%"', svg, count=1)

    tips = {t["gid"]: t for t in cartoon.variant_tooltips()}
    width_px = 1150
    height = int(width_px * aspect) + 190

    html = _INTERACTIVE_HTML.replace("__SVG__", svg).replace(
        "__DATA__", json.dumps(tips))
    components.html(html, height=height, scrolling=False)


_INTERACTIVE_HTML = """
<style>
  body { margin:0; font-family: -apple-system, "Segoe UI", Roboto, sans-serif; }
  #fig { position:relative; width:100%; }
  #fig svg { width:100%; height:auto; display:block; }
  g[id^="variant-"] { cursor:pointer; }
  g[id^="variant-"] use { transition: stroke-width .08s ease; }
  g[id^="variant-"].hot use,
  g[id^="variant-"].pinned use { stroke:#111 !important; stroke-width:1.6 !important; }
  #tip {
    position:fixed; display:none; z-index:9999; pointer-events:none;
    background:#1f2933; color:#fff; padding:7px 10px; border-radius:6px;
    font-size:12px; line-height:1.45; max-width:280px;
    box-shadow:0 4px 14px rgba(0,0,0,.28);
  }
  #tip b { font-size:13px; }
  #tip .sub { color:#c3ccd6; }
  #pins { margin:10px 2px 0; font-size:13px; color:#31333F; }
  #pins .hdr { color:#6b7480; margin-bottom:6px; }
  .chip {
    display:inline-flex; align-items:center; gap:6px; cursor:pointer;
    border:1px solid #d7dbe0; border-radius:14px; padding:3px 10px;
    margin:0 6px 6px 0; background:#fff;
  }
  .chip:hover { background:#f4f6f8; }
  .dot { width:9px; height:9px; border-radius:50%; display:inline-block; }
  .x { color:#98a2ad; font-weight:700; }
</style>

<div id="fig">__SVG__</div>
<div id="tip"></div>
<div id="pins"></div>

<script>
const DATA = __DATA__;
const tip  = document.getElementById('tip');
const pins = document.getElementById('pins');
const pinned = new Set();

function rowsFor(d) {
  const bits = [];
  if (d.cdna && d.cdna !== d.label) bits.push(d.cdna);
  if (d.category) bits.push(d.category);
  const loc = [d.exon, d.codon].filter(Boolean).join(' &middot; ');
  let html = '<b>' + d.label + '</b>';
  if (bits.length) html += '<br><span class="sub">' + bits.join('<br>') + '</span>';
  if (loc)      html += '<br><span class="sub">' + loc + '</span>';
  if (d.position) html += '<br><span class="sub">' + d.position + '</span>';
  if (d.count > 1) html += '<br><span class="sub">seen ' + d.count + '&times;</span>';
  if (d.note)   html += '<br><span class="sub">' + d.note + '</span>';
  return html;
}

function renderPins() {
  if (!pinned.size) { pins.innerHTML = ''; return; }
  let h = '<div class="hdr">Pinned (click a chip to remove)</div>';
  pinned.forEach(gid => {
    const d = DATA[gid];
    h += '<span class="chip" data-gid="' + gid + '">'
       + '<span class="dot" style="background:' + d.color + '"></span>'
       + d.label + (d.category ? ' <span class="sub">' + d.category + '</span>' : '')
       + ' <span class="x">&times;</span></span>';
  });
  pins.innerHTML = h;
  pins.querySelectorAll('.chip').forEach(ch => {
    ch.onclick = () => { toggle(ch.dataset.gid); };
  });
}

function toggle(gid) {
  const el = document.getElementById(gid);
  if (pinned.has(gid)) { pinned.delete(gid); el && el.classList.remove('pinned'); }
  else { pinned.add(gid); el && el.classList.add('pinned'); }
  renderPins();
}

document.querySelectorAll('g[id^="variant-"]').forEach(el => {
  const d = DATA[el.id];
  if (!d) return;
  el.addEventListener('mouseenter', e => {
    el.classList.add('hot');
    tip.innerHTML = rowsFor(d);
    tip.style.display = 'block';
  });
  el.addEventListener('mousemove', e => {
    const pad = 14;
    let x = e.clientX + pad, y = e.clientY + pad;
    const r = tip.getBoundingClientRect();
    if (x + r.width  > window.innerWidth)  x = e.clientX - r.width  - pad;
    if (y + r.height > window.innerHeight) y = e.clientY - r.height - pad;
    tip.style.left = x + 'px';
    tip.style.top  = y + 'px';
  });
  el.addEventListener('mouseleave', () => {
    el.classList.remove('hot');
    tip.style.display = 'none';
  });
  el.addEventListener('click', () => toggle(el.id));
});
</script>
"""


def style_from_state() -> Style:
    """Build a Style from whatever is in st.session_state."""
    st_obj = Style()
    for f in fields(Style):
        key = f"sty_{f.name}"
        if key in st.session_state:
            setattr(st_obj, f.name, st.session_state[key])
    return st_obj


#: where a queued preset waits between reruns
PENDING_STYLE = "_pending_style"


def queue_style(preset: Dict[str, Any]) -> None:
    """
    Ask for a preset to be applied on the next run.

    Streamlit forbids writing to a widget's key once that widget has been
    created in the current run, and the preset buttons necessarily sit below
    the sliders they want to change. So stash the values under a plain key and
    let :func:`apply_pending_style` install them at the top of the next run,
    before any widget exists.
    """
    st.session_state[PENDING_STYLE] = dict(preset)


def apply_pending_style() -> None:
    """Install a queued preset. Must run before any widget is created."""
    preset = st.session_state.pop(PENDING_STYLE, None)
    if not preset:
        return
    known = {f.name for f in fields(Style)}
    for k, v in preset.items():
        if k in known:
            st.session_state[f"sty_{k}"] = v


def seed_state(preset: Dict[str, Any]) -> None:
    """Backwards-compatible alias: queue a preset for the next run."""
    queue_style(preset)


def annotations_from_editor(df) -> List[Dict[str, Any]]:
    out = []
    if df is None:
        return out
    for _, row in df.iterrows():
        try:
            start = int(row["start"])
        except (TypeError, ValueError):
            continue
        try:
            end = int(row["end"])
        except (TypeError, ValueError):
            end = start
        out.append({
            "start": start, "end": end,
            "label": ("" if _isblank(row.get("label")) else str(row["label"])),
            "color": ("#E8A33D" if _isblank(row.get("color"))
                      else str(row["color"])),
            "style": ("box" if _isblank(row.get("style")) else str(row["style"])),
            "row": 0 if _isblank(row.get("row")) else int(row["row"]),
            "coords": "browser",
        })
    return out


def _isblank(v: Any) -> bool:
    if v is None:
        return True
    if isinstance(v, float) and v != v:
        return True
    return str(v).strip() == ""


# --------------------------------------------------------------------------- #
#  App
# --------------------------------------------------------------------------- #

def check_engine_version() -> None:
    """
    Refuse to run against a stale ucsc_gene_cartoon.py.

    Updating one file but not the other is easy to do when uploading through
    the GitHub web interface, and the result is baffling: new controls appear
    in the sidebar but silently do nothing, because the old Style dataclass
    has no field for them. Better to say so plainly.
    """
    raw = getattr(ugc, "__version__", None)
    try:
        got = tuple(int(p) for p in str(raw).split(".")[:3])
    except (ValueError, AttributeError):
        got = (0, 0, 0)
    if got >= ENGINE_MIN:
        return

    want = ".".join(str(p) for p in ENGINE_MIN)
    loaded_from = getattr(ugc, "__file__", "unknown location")
    app_dir = str(Path(__file__).resolve().parent)
    same_folder = str(Path(loaded_from).resolve().parent) == app_dir

    st.error(
        f"**`ucsc_gene_cartoon.py` is out of date.** This app needs version "
        f"{want} or newer, but the copy it loaded reports "
        f"`{raw if raw else 'no version at all'}`."
    )
    st.markdown(
        f"""
**Python loaded this file:**

```
{loaded_from}
```

**app.py is in this folder:**

```
{app_dir}
```
"""
    )
    if not same_folder:
        st.warning(
            "Those are different folders — so there are **two copies** of "
            "`ucsc_gene_cartoon.py` and the wrong one is winning. Delete the "
            "one that isn't next to `app.py`."
        )
    else:
        st.markdown(
            """
Both are in the same folder, so the file simply hasn't been replaced yet.
**Check it directly on GitHub:** open `ucsc_gene_cartoon.py` in your
repository and look at about line 30. The current version has:

```python
__version__ = "1.5.1"
```

If that line isn't there, the upload didn't land. Common reasons:

- the files went into a **subfolder** instead of the top level of the repo
- the commit went to a **different branch** than the one Streamlit deploys
- **Commit changes** at the bottom of the upload page was never clicked

Re-upload via **Add file → Upload files**, drop the file at the top level,
and commit. Then, if it still shows old, click **Manage app** at the bottom
right of this page and choose **Reboot app** to force a clean rebuild.
"""
        )
    st.stop()


def main() -> None:
    st.set_page_config(page_title="Gene Cartoon Studio",
                       page_icon="🧬", layout="wide")
    check_engine_version()
    # before any widget is instantiated, or Streamlit refuses the assignment
    apply_pending_style()
    st.title("Gene Cartoon Studio")
    st.caption("Publication-ready gene diagrams from live UCSC Genome Browser "
               "data. Type a gene, adjust, download.")

    with st.expander("New here? Start with this", expanded=False):
        st.markdown(
            """
            1. **Type a gene symbol** in the left-hand panel — `TP53`, `CFTR`,
               `BRCA1`. The diagram appears below within a second or two.
            2. **Change how it looks** using *Appearance* on the left: colours,
               exon shape, how introns are drawn, what's labelled.
            3. **Add your variants** in the *Variants* tab — upload a
               spreadsheet of cDNA positions (there's a downloadable example)
               and they appear as lollipops above the gene.
            4. **Download** with the buttons under the figure. Use **SVG** for
               journal submission: it stays sharp at any size and the text
               remains editable in Illustrator or Inkscape.

            Nothing you do here can break anything — change a control, look at
            the figure, change it back.
            """
        )

    # ---------------- gene selection ---------------- #
    with st.sidebar:
        st.header("Gene")
        gene = st.text_input("Gene symbol", value="TP53",
                             help="e.g. TP53, CFTR, BRCA1").strip()
        genome = st.selectbox("Assembly", ASSEMBLIES, index=0)
        tracks = ["auto"] + available_tracks(genome)
        track = st.selectbox(
            "Annotation track", tracks, index=0,
            help="'auto' prefers MANE Select (one clean transcript per gene).")

    if not gene:
        st.info("Enter a gene symbol in the sidebar to begin.")
        return

    try:
        with st.spinner(f"Fetching {gene} from UCSC…"):
            txs, links, used_track = fetch_models(gene, genome, track)
    except UCSCError as e:
        st.error(str(e))
        st.stop()
        return
    except Exception as e:                             # pragma: no cover
        st.error(f"Unexpected problem talking to UCSC: {e}")
        st.code(traceback.format_exc())
        st.stop()
        return

    # ---------------- transcript picking ---------------- #
    with st.sidebar:
        st.header("Transcripts")
        default = select_transcripts(txs, gene, "longest-coding", None)
        by_name = {t.name: t for t in txs}
        order = sorted(by_name, key=lambda n: (-by_name[n].spliced_length, n))
        chosen = st.multiselect(
            f"Isoforms in {used_track}", order,
            default=[t.name for t in default],
            format_func=lambda n: (
                f"{n} — {len(by_name[n].exons)} exons, "
                f"{by_name[n].spliced_length:,} nt"
                f"{'' if by_name[n].coding else ' (non-coding)'}"),
            help="Tick more than one to stack isoforms.")
        drawn = [by_name[n] for n in chosen] or default
        primary = drawn[0]
        st.caption(f"cDNA positions are mapped against **{primary.name}**.")

        st.header("Focus")
        n_ex = len(primary.exons)
        scope = st.radio("Show", ["Whole gene", "Single exon", "Exon range"],
                         horizontal=False, key="focus_scope")
        exon_spec = None
        if scope == "Single exon":
            exon_spec = st.selectbox(
                "Exon", list(range(1, n_ex + 1)), key="focus_exon",
                format_func=lambda n: (
                    f"exon {n}  ({primary.exon_span(n)[1] - primary.exon_span(n)[0]:,} bp)"))
        elif scope == "Exon range":
            lo, hi = st.select_slider(
                "Exons", options=list(range(1, n_ex + 1)),
                value=(1, min(3, n_ex)), key="focus_range")
            exon_spec = (lo, hi)
        flank = 200
        if exon_spec is not None:
            flank = st.slider("Flanking intron (bp)", 0, 5000, 200, 50,
                              key="focus_flank",
                              help="How much intron to show either side.")

    # ---------------- style controls ---------------- #
    with st.sidebar:
        st.header("Appearance")

        with st.expander("Colours & shapes", expanded=True):
            st.color_picker("CDS", key="sty_cds_color", value="#2C6FA6")
            st.color_picker("CDS outline", key="sty_cds_edge", value="#123A5C")
            st.color_picker("UTR", key="sty_utr_color", value="#BBD5E8")
            st.color_picker("Non-coding exon", key="sty_noncoding_color",
                            value="#B7B7B7")
            st.color_picker("Intron", key="sty_intron_color", value="#4A4A4A")
            st.selectbox("Exon shape", ["box", "round"], key="sty_exon_shape")
            st.slider("Corner rounding", 0.0, 0.5, key="sty_corner_radius",
                      value=0.30, step=0.05)
            st.selectbox("Intron style", ["line", "chevron_line", "angled"],
                         key="sty_intron_style")

        with st.expander("Layout & scale"):
            st.selectbox(
                "Intron scale", ["compress", "linear", "equal"],
                key="sty_intron_mode",
                help="compress = fixed-width introns (usual for figures); "
                     "linear = true genomic scale.")
            st.slider("Figure width (in)", 4.0, 20.0, key="sty_figure_width_in",
                      value=9.0, step=0.5)
            st.slider("Row height (in)", 0.4, 2.0, key="sty_row_height_in",
                      value=0.85, step=0.05)
            st.slider("CDS height", 0.1, 0.9, key="sty_cds_height",
                      value=0.42, step=0.02)
            st.slider("UTR height", 0.05, 0.9, key="sty_utr_height",
                      value=0.22, step=0.02)
            st.slider("Base font size", 5.0, 16.0, key="sty_font_size",
                      value=9.0, step=0.5)
            st.checkbox("Exon numbers", key="sty_show_exon_numbers", value=True)
            st.selectbox("Number position", ["above", "inside", "below"],
                         key="sty_exon_number_position")
            st.checkbox("Title", key="sty_show_title", value=True)
            st.checkbox("Locus subtitle", key="sty_show_subtitle", value=True)
            st.checkbox("Legend", key="sty_show_legend", value=True)
            st.checkbox("5'→3' arrow", key="sty_show_strand_arrow", value=True)
            st.selectbox("Axis", ["auto", "coords", "scalebar", "none"],
                         key="sty_axis_mode")
            st.checkbox("Put 5' end on the left", value=True,
                        key="sty_orient_five_prime_left",
                        help="Flips minus-strand genes so they read 5'→3'.")
            st.selectbox("Background", ["none", "white"], key="sty_background",
                         help="'none' = transparent.")

        with st.expander("Variants"):
            cc1, cc2 = st.columns(2)
            with cc1:
                if st.button("Compact", use_container_width=True,
                             help="Small markers, tight rows — for large "
                                  "cohorts. Nothing is hidden."):
                    queue_style(COMPACT_PRESET)
                    st.rerun()
            with cc2:
                if st.button("Roomy", use_container_width=True,
                             help="Back to the spacious defaults."):
                    queue_style(ROOMY_PRESET)
                    st.rerun()
            st.selectbox(
                "Display", ["stacked", "lanes", "lollipop"],
                key="sty_variant_style",
                help="stacked = markers pile up above their position. "
                     "lanes = one labelled row per category, so mutation "
                     "types can be compared along the gene. "
                     "lollipop = classic stems.")
            st.slider("Gap between lanes", 0.0, 0.8, key="sty_lane_gap",
                      value=0.16, step=0.02,
                      help="Only used by the 'lanes' display.")
            st.selectbox("Marker", ["o", "s", "D", "v", "^"],
                         key="sty_variant_marker",
                         format_func=lambda m: {"o": "circle", "s": "square",
                                                "D": "diamond",
                                                "v": "triangle down",
                                                "^": "triangle up"}[m])
            st.slider("Marker size", 6.0, 200.0, key="sty_variant_head_size",
                      value=42.0, step=2.0,
                      help="Row spacing follows this, so turning it down "
                           "genuinely shortens the figure.")
            st.slider("Row spacing (× marker)", 0.8, 3.0,
                      key="sty_variant_row_spacing", value=1.30, step=0.05,
                      help="Vertical pitch between rows, as a multiple of the "
                           "marker diameter. 1.0 = markers just touching.")
            st.slider("Horizontal packing (× marker)", 0.5, 3.0,
                      key="sty_variant_pack", value=1.25, step=0.05,
                      help="How much clearance a marker claims sideways. "
                           "Lower packs more per row, so fewer rows are "
                           "needed; below 1.0 markers may touch.")
            st.slider("Gap above gene", 0.0, 1.5, key="sty_variant_base_gap",
                      value=0.20, step=0.05)
            st.checkbox("Scale marker by count", value=True,
                        key="sty_variant_scale_by_count")
            st.checkbox("Show variant labels", value=True,
                        key="sty_show_variant_labels")
            st.slider("Label size", 4.0, 14.0, key="sty_variant_label_size",
                      value=6.5, step=0.5)
            st.slider("Stem height (lollipop only)", 0.2, 2.0,
                      key="sty_variant_stem_height", value=0.55, step=0.05)

    style = style_from_state()

    # ---------------- data panels ---------------- #
    tab_variants, tab_annot, tab_data, tab_preset = st.tabs(
        ["Variants", "Annotations", "Transcript data", "Style presets"])

    # -- variants -- #
    variant_objs: List[Variant] = []
    failed: List[Variant] = []

    with tab_variants:
        if not HAVE_TABLES:
            st.warning("Install pandas and openpyxl to import variant files: "
                       "`pip install pandas openpyxl`")
        else:
            st.markdown(
                "Upload a **CSV or Excel** file of variants with cDNA (HGVS "
                "`c.`) positions. Columns are detected automatically — you can "
                "override them below.")
            c1, c2 = st.columns([3, 1])
            with c1:
                up = st.file_uploader("Variant table",
                                      type=["csv", "tsv", "txt", "xlsx", "xls"])
            with c2:
                st.download_button(
                    "Example file", data=_template_bytes(),
                    file_name="variant_template.xlsx",
                    mime="application/vnd.openxmlformats-officedocument."
                         "spreadsheetml.sheet",
                    help="A filled-in example showing the expected columns.")

            if up is not None:
                variant_objs, failed = _variant_panel(up, primary)

    # -- annotations -- #
    with tab_annot:
        st.markdown(
            "Domains, regions and markers drawn above the gene. Positions are "
            "**genomic, 1-based** (as shown in the UCSC browser).")
        if HAVE_TABLES:
            blank = pd.DataFrame([{"start": None, "end": None, "label": "",
                                   "color": "#E8A33D", "style": "box",
                                   "row": 0}])
            edited = st.data_editor(
                blank, num_rows="dynamic", key="annot_editor",
                use_container_width=True,
                column_config={
                    "style": st.column_config.SelectboxColumn(
                        "style", options=["box", "bracket", "marker"]),
                    "row": st.column_config.NumberColumn(
                        "row", min_value=0, max_value=4, step=1),
                })
            annotations = annotations_from_editor(edited)
            st.caption(f"{primary.chrom}:{primary.tx_start + 1:,}–"
                       f"{primary.tx_end:,} is the drawable range.")
        else:
            annotations = []
            st.info("Install pandas to use the annotation editor.")

    # -- transcript data -- #
    with tab_data:
        rows = [{"transcript": t.name, "exons": len(t.exons),
                 "exonic nt": t.spliced_length,
                 "CDS nt": sum(e - s for s, e, k in t.segments() if k == "cds"),
                 "start": t.tx_start + 1, "end": t.tx_end,
                 "strand": t.strand, "coding": t.coding} for t in drawn]
        if HAVE_TABLES:
            st.dataframe(pd.DataFrame(rows), use_container_width=True,
                         hide_index=True)
        else:
            st.write(rows)
        st.markdown(f"**Sequence (UCSC API):** {links['sequence']}")
        st.markdown(f"**Browser view:** {links['browser']}")

    # -- presets -- #
    with tab_preset:
        st.markdown("Save the current look as a house style, and load it back "
                    "later so every figure in a paper matches.")
        c1, c2 = st.columns(2)
        with c1:
            st.download_button(
                "Download style JSON",
                data=json.dumps(asdict(style), indent=2),
                file_name="lab_style.json", mime="application/json")
        with c2:
            preset = st.file_uploader("Load style JSON", type=["json"],
                                      key="style_upload")
            if preset is not None and st.button("Apply loaded style"):
                try:
                    queue_style(json.load(preset))
                except (json.JSONDecodeError, ValueError) as e:
                    st.error(f"Not a valid style file: {e}")
                else:
                    st.rerun()

    # ---------------- render ---------------- #
    try:
        bounds, focus_label = ugc.resolve_focus(primary, exon_spec, flank)
        cartoon = GeneCartoon(drawn, style, gene, genome, used_track,
                              annotations=annotations, links=links,
                              variants=variant_objs, focus=bounds,
                              focus_label=focus_label)
        if bounds:
            n_off = len(cartoon.offscreen_variants)
            st.info(
                f"Showing **{focus_label}** of {primary.name} — "
                f"{primary.chrom}:{bounds[0] + 1:,}–{bounds[1]:,}"
                + (f". {len(cartoon.variants)} variant(s) here, "
                   f"{n_off} elsewhere in the gene." if variant_objs else "."))
        if variant_objs:
            interactive = st.checkbox(
                "Interactive figure — hover a variant for details, click to pin",
                value=True, key="interactive_preview")
        else:
            interactive = False
        show_figure(cartoon, interactive=interactive)
    except Exception as e:                             # pragma: no cover
        st.error(f"Could not draw the figure: {e}")
        st.code(traceback.format_exc())
        return

    st.subheader("Download")
    cols = st.columns(len(FORMATS) + 1)
    for col, (label, ext, mime) in zip(cols, FORMATS):
        with col:
            try:
                data = cartoon.to_bytes("tiff" if ext == "tif" else ext)
            except Exception as e:                     # pragma: no cover
                st.caption(f"{label}: {e}")
                continue
            st.download_button(label, data=data,
                               file_name=f"{gene}_{genome}.{ext}", mime=mime,
                               use_container_width=True)
    with cols[-1]:
        st.number_input("DPI", 72, 1200, key="sty_dpi", value=400, step=50,
                        help="Applies to PNG and TIFF.")

    if failed:
        st.warning(f"{len(failed)} variant(s) could not be placed — see the "
                   "Variants tab.")

    st.divider()
    st.caption(
        f"Gene models from the UCSC Genome Browser REST API "
        f"({used_track}, {genome}). Please cite UCSC if you publish a figure "
        "made with this tool."
    )
    with st.expander("What happens to an uploaded file?"):
        st.markdown(
            """
            An uploaded spreadsheet is held in memory for your session only.
            This app never writes it to disk, and other users of the app
            cannot see it.

            It is still read by the server this app runs on. **If your file
            contains identifiable patient data, do not upload it to a public
            deployment** — run the app locally instead
            (`python run_gui.py`), where nothing leaves your machine except
            the gene lookup to UCSC.
            """
        )


def category_controls(placed: List[Variant],
                      column: Optional[str]) -> List[Variant]:
    """
    Per-category show/hide and colour.

    Returns only the variants the user has left switched on, with an explicit
    colour stamped on each so both the markers and the legend use it.
    """
    cats: List[str] = []
    for v in placed:
        c = v.category or "(no category)"
        if c not in cats:
            cats.append(c)
    if not cats:
        return placed
    cats.sort(key=lambda c: (-sum(1 for v in placed
                                  if (v.category or "(no category)") == c), c))

    palette = Style().variant_palette
    label = column or "category"
    st.markdown(f"**Show / colour by `{label}`** — untick to hide a group.")

    c1, c2, _ = st.columns([1, 1, 3])
    with c1:
        if st.button("Select all", use_container_width=True):
            for c in cats:
                st.session_state[f"cat_on_{c}"] = True
    with c2:
        if st.button("Select none", use_container_width=True):
            for c in cats:
                st.session_state[f"cat_on_{c}"] = False

    chosen: Dict[str, str] = {}
    for i, cat in enumerate(cats):
        n = sum(1 for v in placed if (v.category or "(no category)") == cat)
        row = st.columns([0.5, 4, 1.2])
        with row[0]:
            on = st.checkbox(" ", value=True, key=f"cat_on_{cat}",
                             label_visibility="collapsed")
        with row[1]:
            st.markdown(f"{cat} &nbsp;<span style='color:#888'>({n})</span>",
                        unsafe_allow_html=True)
        with row[2]:
            colour = st.color_picker(
                " ", value=palette[i % len(palette)], key=f"cat_col_{cat}",
                label_visibility="collapsed")
        if on:
            chosen[cat] = colour

    kept = []
    for v in placed:
        cat = v.category or "(no category)"
        if cat in chosen:
            v.color = chosen[cat]
            kept.append(v)

    hidden = len(placed) - len(kept)
    if hidden:
        st.caption(f"{hidden} variant(s) hidden by the filter above.")
    if not kept:
        st.warning("Every group is hidden, so no variants will be drawn.")
    return kept


def _variant_panel(upload, transcript: Transcript):
    """Column mapping + placement report. Returns (placed, failed)."""
    sheets = varmod.excel_sheet_names(upload)
    sheet = None
    if len(sheets) > 1:
        sheet = st.selectbox("Worksheet", sheets)

    try:
        df = varmod.read_table(upload, sheet=sheet)
    except Exception as e:
        st.error(f"Could not read that file: {e}")
        return [], []

    df = df.dropna(how="all")
    df.columns = [str(c).strip() for c in df.columns]
    if df.empty:
        st.error("That file has no rows.")
        return [], []

    detected = varmod.detect_columns(df)
    st.markdown("**Columns** — detected automatically; change if wrong.")
    choices = ["(none)"] + list(df.columns)
    cols = st.columns(4)
    mapping: Dict[str, Optional[str]] = {}
    for col_box, role, label in zip(
            cols * 2,
            ["cdna", "label", "category", "count"],
            ["cDNA position", "Label", "Category / class", "Count"]):
        with col_box:
            cur = detected.get(role)
            idx = choices.index(cur) if cur in choices else 0
            pick = st.selectbox(label, choices, index=idx, key=f"map_{role}")
            mapping[role] = None if pick == "(none)" else pick
    for role in ("genomic", "color", "note"):
        mapping[role] = detected.get(role)

    if not mapping["cdna"] and not mapping["genomic"]:
        st.error("Pick the column holding cDNA positions (e.g. `c.743G>A`).")
        return [], []

    vf = varmod.build_variants(df, mapping)
    for w in vf.warnings:
        st.warning(w)

    placed, failed = place_variants(vf.variants, transcript)
    st.success(f"Placed {len(placed)} of {len(vf.variants)} variant(s) on "
               f"{transcript.name}.")

    all_placed = list(placed)
    placed = category_controls(placed, mapping.get("category"))
    shown = {id(v) for v in placed}

    mapper = CDNAMapper(transcript)
    report = []
    for v in all_placed:
        report.append({"label": v.label, "cDNA": v.cdna,
                       f"{transcript.chrom} (1-based)": v.genomic,
                       "exon": v.exon, "codon": v.protein,
                       "category": v.category, "count": v.count,
                       "status": "shown" if id(v) in shown else "hidden"})
    for v in failed:
        report.append({"label": v.label, "cDNA": v.cdna,
                       f"{transcript.chrom} (1-based)": None,
                       "exon": None, "codon": None,
                       "category": v.category, "count": v.count,
                       "status": v.error})
    st.dataframe(pd.DataFrame(report), use_container_width=True,
                 hide_index=True)
    st.caption(
        f"{transcript.name}: c.1 is at {transcript.chrom}:"
        f"{mapper.to_genomic('c.1'):,} and the transcript is "
        f"{mapper.tx_length:,} nt. Negative positions (c.-30) sit in the "
        "5'UTR, c.* positions after the stop codon, and c.123+4 is intronic.")

    if failed:
        st.info("Rows that didn't place are usually protein (p.) notation, "
                "numbering from a different transcript, or a typo. Everything "
                "else still drew fine.")
    return placed, failed


@st.cache_data(show_spinner=False)
def _template_bytes() -> bytes:
    buf = io.BytesIO()
    rows = [
        ("c.524G>A", "R175H", "Pathogenic", 9, "hotspot"),
        ("c.659A>G", "Y220C", "Pathogenic", 5, ""),
        ("c.743G>A", "R248Q", "Pathogenic", 12, "hotspot"),
        ("c.818G>A", "R273H", "Pathogenic", 8, "hotspot"),
        ("c.215C>G", "P72R", "Benign", 1, "common polymorphism"),
        ("c.376-2A>G", "splice acceptor", "VUS", 1, "intronic"),
        ("c.-28C>T", "5'UTR", "VUS", 1, "upstream of ATG"),
        ("c.*38G>A", "3'UTR", "VUS", 2, "after stop codon"),
    ]
    df = pd.DataFrame(rows, columns=["cDNA", "Protein", "Classification",
                                     "Count", "Notes"])
    df.to_excel(buf, index=False, sheet_name="variants")
    return buf.getvalue()


if __name__ == "__main__":
    main()
