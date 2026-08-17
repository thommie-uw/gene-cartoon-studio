#!/usr/bin/env python3
"""
Headless test for app.py.

Streamlit can't be installed in this environment, so we stand in a mock that
implements the widget API and returns the defaults each widget was given.
Everything that isn't a widget -- fetching, column mapping, cDNA placement,
figure rendering, download bytes -- runs for real.

This catches wrong Streamlit call signatures, bad state keys, and any logic
error in the app, without a browser.
"""

import io
import sys
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

failures = []
calls = {}


def check(name, cond, detail=""):
    print(f"{'PASS' if cond else 'FAIL'}  {name}{'  ' + detail if detail else ''}")
    if not cond:
        failures.append(name)


# --------------------------------------------------------------------------- #
#  Mock streamlit
# --------------------------------------------------------------------------- #

class _Ctx:
    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


class StreamlitAPIException(Exception):
    """Stand-in for the real thing."""


class RerunException(Exception):
    """Raised by the mock's st.rerun(), like Streamlit's own control flow."""


class _SessionState(dict):
    """
    Session state that enforces Streamlit's real rule:

    once a widget with key K has been created during this run, assigning to
    st.session_state[K] raises. Without this the mock silently accepts writes
    that blow up in production -- which is exactly how the preset buttons
    shipped broken.
    """

    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        dict.__setattr__(self, "_widget_keys", set())

    def _claim(self, key):
        if key:
            self._widget_keys.add(key)

    def __setitem__(self, k, v):
        if k in self._widget_keys:
            raise StreamlitAPIException(
                f"st.session_state.{k} cannot be modified after the widget "
                f"with key {k} is instantiated."
            )
        super().__setitem__(k, v)

    def __getattr__(self, k):
        try:
            return self[k]
        except KeyError as e:
            raise AttributeError(k) from e

    def __setattr__(self, k, v):
        self[k] = v


class _ColumnConfig:
    def SelectboxColumn(self, *a, **k):
        return {"kind": "selectbox", "args": a, "kwargs": k}

    def NumberColumn(self, *a, **k):
        return {"kind": "number", "args": a, "kwargs": k}

    def TextColumn(self, *a, **k):
        return {"kind": "text", "args": a, "kwargs": k}


class MockStreamlit(types.ModuleType):
    """Records calls; widgets return their default/first option.

    ``overrides`` simulates the user changing a control: {widget_key: value}.
    """

    def __init__(self, uploads=None, overrides=None, buttons=None):
        super().__init__("streamlit")
        self.session_state = _SessionState()
        self.sidebar = _Ctx()
        self.column_config = _ColumnConfig()
        self._uploads = uploads or {}
        self._overrides = overrides or {}
        self._buttons = buttons or {}
        self.errors, self.warnings_, self.downloads, self.images = [], [], [], []
        self.embeds = []
        self.stopped = False
        self.rerun_requested = False
        comp = types.ModuleType("streamlit.components.v1")
        comp.html = lambda h, height=None, **k: self.embeds.append((h, height))
        self.components = types.SimpleNamespace(v1=comp)
        sys.modules["streamlit.components"] = types.ModuleType(
            "streamlit.components")
        sys.modules["streamlit.components.v1"] = comp

    # -- record every attribute access we didn't define -- #
    def _note(self, name):
        calls[name] = calls.get(name, 0) + 1

    # layout
    def set_page_config(self, **k): self._note("set_page_config")
    def title(self, *a, **k): self._note("title")
    def caption(self, *a, **k): self._note("caption")
    def header(self, *a, **k): self._note("header")
    def subheader(self, *a, **k): self._note("subheader")
    def markdown(self, *a, **k): self._note("markdown")
    def write(self, *a, **k): self._note("write")
    def code(self, *a, **k): self._note("code")
    def info(self, *a, **k): self._note("info")
    def success(self, *a, **k): self._note("success")

    def error(self, msg, *a, **k):
        self._note("error"); self.errors.append(str(msg))

    def warning(self, msg, *a, **k):
        self._note("warning"); self.warnings_.append(str(msg))

    def expander(self, *a, **k): self._note("expander"); return _Ctx()
    def spinner(self, *a, **k): self._note("spinner"); return _Ctx()
    def columns(self, spec, **k):
        self._note("columns")
        n = spec if isinstance(spec, int) else len(spec)
        return [_Ctx() for _ in range(n)]

    def tabs(self, names, **k):
        self._note("tabs"); return [_Ctx() for _ in names]

    def stop(self):
        self._note("stop"); self.stopped = True

    def rerun(self):
        self._note("rerun")
        raise RerunException()
    def divider(self): self._note("divider")

    # widgets -- return the declared default
    def _remember(self, key, value):
        if key and key in self._overrides:
            value = self._overrides[key]
        elif key and key in self.session_state:
            # real Streamlit: a value already in session_state wins over the
            # widget's declared default. This is what makes presets stick.
            value = self.session_state[key]
        if key:
            # write before claiming: a widget may legitimately seed its own key
            dict.__setitem__(self.session_state, key, value)
            self.session_state._claim(key)
        return value

    def text_input(self, label, value="", key=None, **k):
        self._note("text_input"); return self._remember(key, value)

    def selectbox(self, label, options, index=0, key=None, **k):
        self._note("selectbox")
        options = list(options)
        return self._remember(key, options[index] if options else None)

    def multiselect(self, label, options, default=None, key=None, **k):
        self._note("multiselect")
        return self._remember(key, list(default or []))

    def color_picker(self, label, value="#000000", key=None, **k):
        self._note("color_picker"); return self._remember(key, value)

    def slider(self, label, mn=0, mx=1, value=None, key=None, **k):
        self._note("slider")
        return self._remember(key, value if value is not None else mn)

    def number_input(self, label, mn=None, mx=None, value=None, key=None, **k):
        self._note("number_input")
        return self._remember(key, value if value is not None else (mn or 0))

    def radio(self, label, options, index=0, key=None, **k):
        self._note("radio")
        options = list(options)
        return self._remember(key, options[index] if options else None)

    def select_slider(self, label, options=None, value=None, key=None, **k):
        self._note("select_slider")
        return self._remember(key, value)

    def checkbox(self, label, value=False, key=None, **k):
        self._note("checkbox"); return self._remember(key, value)

    def button(self, label, key=None, **k):
        self._note("button")
        return bool(self._buttons.get(label, False))

    def file_uploader(self, label, type=None, key=None, **k):
        self._note("file_uploader")
        return self._uploads.get(key or label)

    def download_button(self, label, data=None, file_name=None, key=None, **k):
        self._note("download_button")
        self.downloads.append((label, file_name, len(data) if data else 0))
        return False

    def image(self, data, **k):
        self._note("image"); self.images.append(len(data) if data else 0)

    def dataframe(self, df, **k): self._note("dataframe")

    def data_editor(self, df, **k):
        self._note("data_editor"); return df

    def cache_data(self, *a, **k):
        """Pass-through decorator, with or without arguments."""
        if a and callable(a[0]):
            return a[0]

        def deco(fn):
            return fn
        return deco


# --------------------------------------------------------------------------- #
#  Stub the UCSC network layer with real captured payloads
# --------------------------------------------------------------------------- #

TP53_ROWS = [
    {"name": "NM_000546.6", "chrom": "chr17", "strand": "-",
     "txStart": 7668420, "txEnd": 7687490, "cdsStart": 7669608,
     "cdsEnd": 7676594,
     "exonStarts": "7668420,7670608,7673534,7673700,7674180,7674858,7675052,7675993,7676381,7676520,7687376,",
     "exonEnds": "7669690,7670715,7673608,7673837,7674290,7674971,7675236,7676272,7676403,7676622,7687490,",
     "name2": "TP53"},
    {"name": "NM_001126112.3", "chrom": "chr17", "strand": "-",
     "txStart": 7668420, "txEnd": 7687490, "cdsStart": 7669608,
     "cdsEnd": 7676594,
     "exonStarts": "7668420,7670608,7673534,7673700,7674180,7674858,7675052,7675993,7676381,7676520,7687376,",
     "exonEnds": "7669690,7670715,7673608,7673837,7674290,7674971,7675236,7676272,7676403,7676619,7687490,",
     "name2": "TP53"},
]


def fake_get(self, endpoint, params, use_cache=True):
    if endpoint == "search":
        return {"positionMatches": [{"trackName": "mane", "matches": [
            {"posName": "TP53", "position": "chr17:7668421-7687490"}]}]}
    if endpoint == "list/tracks":
        return {params.get("genome", "hg38"): {
            "ncbiRefSeqSelect": {"type": "genePred"},
            "ncbiRefSeqCurated": {"type": "genePred"},
            "gap": {"type": "bed 3"}}}
    if endpoint == "getData/track":
        return {params["track"]: [dict(r) for r in TP53_ROWS]}
    if endpoint == "getData/sequence":
        return {"dna": "acgt" * 25}
    raise AssertionError(f"unexpected endpoint {endpoint}")


# --------------------------------------------------------------------------- #
#  Run
# --------------------------------------------------------------------------- #

def build_variant_upload():
    import pandas as pd
    buf = io.BytesIO()
    pd.DataFrame({
        "HGVSc": ["c.524G>A", "c.743G>A", "c.818G>A", "c.376-2A>G",
                  "c.-28C>T", "c.*38G>A", "p.Arg248Gln"],
        "HGVSp": ["R175H", "R248Q", "R273H", "splice", "5'UTR", "3'UTR", "bad"],
        "Clinical Significance": ["Pathogenic"] * 3 + ["VUS"] * 3 + ["VUS"],
        "No. of cases": [9, 12, 8, 1, 1, 2, 1],
    }).to_csv(buf, index=False)
    buf.seek(0)
    buf.name = "variants.csv"
    return buf


def run_app(uploads=None, overrides=None, buttons=None, state=None):
    mock = MockStreamlit(uploads=uploads, overrides=overrides, buttons=buttons)
    if state:
        for k, v in state.items():
            dict.__setitem__(mock.session_state, k, v)
    sys.modules["streamlit"] = mock
    for m in ("app",):
        sys.modules.pop(m, None)

    import ucsc_gene_cartoon as ugc
    ugc.UCSCClient._get = fake_get

    import app
    try:
        app.main()
    except RerunException:
        mock.rerun_requested = True
    return mock, app


print("=== app with no variant file ===")
mock, appmod = run_app()
check("app runs end to end", not mock.errors and not mock.stopped,
      "; ".join(mock.errors))
check("figure rendered to the page", mock.images and mock.images[0] > 4000,
      f"{mock.images} bytes")
check("download buttons offered", len(mock.downloads) >= 4,
      ", ".join(f"{d[1]}" for d in mock.downloads if d[1]))
fmts = {d[1].rsplit(".", 1)[-1] for d in mock.downloads if d[1]}
check("all four figure formats downloadable",
      {"svg", "pdf", "png", "tif"} <= fmts, str(sorted(fmts)))
check("every download carries real bytes",
      all(size > 1000 for _, name, size in mock.downloads
          if name and name.endswith(("svg", "pdf", "png", "tif"))))
check("style widgets populate session state",
      mock.session_state.get("sty_cds_color") == "#2C6FA6"
      and mock.session_state.get("sty_intron_mode") == "compress")

print("\n=== app with a variant spreadsheet ===")
up = build_variant_upload()
mock2, _ = run_app(uploads={"Variant table": up})
check("app runs with variants", not mock2.errors and not mock2.stopped,
      "; ".join(mock2.errors))
# with variants present the preview becomes the interactive component
check("figure still renders", bool(mock2.embeds) or bool(mock2.images),
      "interactive embed" if mock2.embeds else "static image")
check("unplaced variant is reported",
      any("could not be placed" in w for w in mock2.warnings_)
      or any("p." in e for e in mock2.errors),
      f"warnings={mock2.warnings_}")

print("\n=== sidebar controls actually reach the renderer ===")
# The bug this guards: a control can appear in the sidebar and be silently
# dropped on the way to the figure, so assert the drawn output really changes.
import ucsc_gene_cartoon as _ugc
STEM = _ugc.Style().variant_stem_color.lower().lstrip("#")


def svg_for(overrides):
    m, a = run_app(uploads={"Variant table": build_variant_upload()},
                   overrides=overrides)
    style = a.style_from_state()
    txs, links, track = a.fetch_models("TP53", "hg38", "auto")
    drawn = _ugc.select_transcripts(txs, "TP53", "longest-coding", None)
    vf = __import__("variants").load_variant_table(build_variant_upload())
    placed, _ = _ugc.place_variants(vf.variants, drawn[0])
    return style, _ugc.GeneCartoon(drawn, style, "TP53", "hg38", track,
                                   variants=placed).to_svg().lower()


st_stacked, svg_stacked = svg_for({"sty_variant_style": "stacked"})
st_lolli, svg_lolli = svg_for({"sty_variant_style": "lollipop"})
check("choosing 'stacked' reaches Style", st_stacked.variant_style == "stacked")
check("choosing 'lollipop' reaches Style", st_lolli.variant_style == "lollipop")
check("'stacked' really draws no stems", svg_stacked.count(STEM) == 0,
      f"{svg_stacked.count(STEM)} stem-coloured strokes")
check("'lollipop' really draws stems", svg_lolli.count(STEM) > 0,
      f"{svg_lolli.count(STEM)} stem-coloured strokes")
check("the two settings produce different figures", svg_stacked != svg_lolli)

st_sq, svg_sq = svg_for({"sty_variant_marker": "s"})
check("marker shape reaches Style", st_sq.variant_marker == "s")
check("marker shape changes the drawing", svg_sq != svg_stacked)

st_col, _ = svg_for({"sty_cds_color": "#7A5195", "sty_intron_mode": "linear",
                     "sty_figure_width_in": 12.0})
check("colour / intron mode / width all reach Style",
      (st_col.cds_color, st_col.intron_mode, st_col.figure_width_in)
      == ("#7A5195", "linear", 12.0))

print("\n=== interactive figure ===")
import json as _json
import re as _re

_m, _a = run_app(uploads={"Variant table": build_variant_upload()})
_vf = __import__("variants").load_variant_table(build_variant_upload())
_txs, _links, _track = _a.fetch_models("TP53", "hg38", "auto")
_drawn = _ugc.select_transcripts(_txs, "TP53", "longest-coding", None)
_placed, _ = _ugc.place_variants(_vf.variants, _drawn[0])
_cart = _ugc.GeneCartoon(_drawn, _ugc.Style(), "TP53", "hg38", _track,
                         variants=_placed)
_m.embeds.clear()
_a.show_interactive(_cart)
check("interactive view is emitted as an HTML component", len(_m.embeds) == 1)
_html, _height = _m.embeds[0]
check("no unfilled template placeholders",
      "__SVG__" not in _html and "__DATA__" not in _html)
_ids = set(_re.findall(r'id="(variant-\d+)"', _html))
_data = _json.loads(_re.search(r'const DATA = (\{.*?\});\n', _html, _re.S).group(1))
check("every marker in the SVG has a tooltip record",
      _ids == set(_data) and len(_ids) == len(_placed),
      f"{len(_ids)} markers, {len(_data)} records")
check("SVG is made fluid so it scales to the page",
      'width="100%"' in _html and 'height="auto"' not in _html.split("<svg")[0])
for _ev in ("mouseenter", "mousemove", "mouseleave", "click"):
    check(f"{_ev} handler is wired", f"'{_ev}'" in _html)
check("click-to-pin logic present",
      "function toggle" in _html and "pinned.add" in _html
      and "pinned.delete" in _html)
check("iframe height is sized from the figure", 300 < _height < 4000, str(_height))
check("tooltip payload survives JSON escaping",
      "</script>" not in _json.dumps(_data))

_static, _ = run_app()
check("a figure with no variants falls back to a static image",
      _static.images and not _static.embeds)

print("\n=== stale-engine guard ===")
_real = _ugc.__version__
_ugc.__version__ = "1.0.0"
mock_old, _ = run_app()
check("an out-of-date ucsc_gene_cartoon.py is caught loudly",
      mock_old.stopped and any("out of date" in e for e in mock_old.errors),
      (mock_old.errors[:1] or ["no error raised"])[0][:80])
_ugc.__version__ = _real
mock_ok, _ = run_app()
check("the current version passes the guard",
      not mock_ok.stopped and not mock_ok.errors)

print("\n=== app-level helpers ===")
style = appmod.style_from_state()
check("style_from_state returns a Style", isinstance(style, appmod.Style))
# a real rerun starts with state already populated, so seed it the way
# Streamlit itself would rather than assigning to a live widget key
dict.__setitem__(appmod.st.session_state, "sty_cds_color", "#123456")
check("style_from_state picks up edits",
      appmod.style_from_state().cds_color == "#123456")

print("\n=== density presets (the Compact / Roomy buttons) ===")
m_c, a_c = run_app(buttons={"Compact": True})
check("clicking Compact does not raise", not m_c.errors, "; ".join(m_c.errors))
check("Compact asks for a rerun", m_c.rerun_requested)
pend = dict.get(m_c.session_state, a_c.PENDING_STYLE)
check("Compact queues the preset instead of writing widget keys",
      pend == a_c.COMPACT_PRESET, str(pend)[:60])
check("Compact did not touch any widget key directly",
      not any(k.startswith("sty_") and dict.get(m_c.session_state, k)
              != a_c.ROOMY_PRESET.get(k[4:], dict.get(m_c.session_state, k))
              for k in a_c.COMPACT_PRESET))

# the next run picks the preset up
m_c2, a_c2 = run_app(state={a_c.PENDING_STYLE: dict(a_c.COMPACT_PRESET)})
check("the queued preset is applied on the next run",
      not m_c2.errors and not m_c2.rerun_requested, "; ".join(m_c2.errors))
applied = a_c2.style_from_state()
check("Compact really reaches the Style",
      (applied.variant_head_size, applied.variant_row_spacing,
       applied.row_height_in) == (14.0, 1.05, 0.65),
      f"{applied.variant_head_size}, {applied.variant_row_spacing}, "
      f"{applied.row_height_in}")
check("the queue is cleared once applied",
      a_c2.PENDING_STYLE not in m_c2.session_state)

m_r, a_r = run_app(buttons={"Roomy": True})
check("clicking Roomy does not raise", not m_r.errors and m_r.rerun_requested)
check("Roomy queues the roomy preset",
      dict.get(m_r.session_state, a_r.PENDING_STYLE) == a_r.ROOMY_PRESET)

m_s, a_s = run_app(state={"_pending_style": {"cds_color": "#ABCDEF"}})
check("a loaded style JSON is applied the same way",
      a_s.style_from_state().cds_color == "#ABCDEF")
check("unknown keys in a loaded style are ignored, not crashed on",
      run_app(state={"_pending_style": {"not_a_field": 1,
                                        "cds_color": "#101010"}})[1]
      .style_from_state().cds_color == "#101010")

import pandas as pd
ann = appmod.annotations_from_editor(pd.DataFrame([
    {"start": 7674220, "end": 7674230, "label": "dom", "color": "#E8A33D",
     "style": "box", "row": 0},
    {"start": None, "end": None, "label": "", "color": "", "style": "",
     "row": 0},                                    # blank row -> dropped
    {"start": 7676594, "end": "", "label": "ATG", "color": "", "style": "marker",
     "row": 1},
]))
check("annotation editor drops blank rows", len(ann) == 2, str(len(ann)))
check("annotation editor defaults a missing end to start",
      ann[1]["start"] == ann[1]["end"] == 7676594)
check("annotation editor uses browser (1-based) coords",
      all(a["coords"] == "browser" for a in ann))
check("annotation editor fills a default colour",
      ann[1]["color"] == "#E8A33D")

check("template file is a real xlsx",
      appmod._template_bytes()[:2] == b"PK",
      f"{len(appmod._template_bytes())} bytes")

print("\n=== launcher ===")
import run_gui
check("run_gui finds app.py", run_gui.APP.exists())
check("run_gui dependency check runs", isinstance(run_gui.missing_packages(), list),
      f"missing here: {run_gui.missing_packages()}")

print()
print(f"widget calls exercised: {sum(calls.values())} across "
      f"{len(calls)} Streamlit APIs")
print(f"{'ALL TESTS PASSED' if not failures else str(len(failures)) + ' FAILURE(S): ' + ', '.join(failures)}")
sys.exit(1 if failures else 0)
