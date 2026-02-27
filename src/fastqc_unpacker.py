#!/usr/bin/env python3
# raw_qc_report.py
# Mimic FastQC table/image style: reuse FastQC-like CSS + embed FastQC original Images/*.png

from __future__ import annotations
import argparse
import base64
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# --------- Parsing fastqc_data.txt ---------

@dataclass
class Module:
    name: str
    status: str  # pass/warn/fail
    header: List[str]
    rows: List[List[str]]  # raw table rows
    cols: List[str]        # table columns (no leading '#')


def parse_fastqc_data_txt(txt_path: Path) -> Dict[str, Module]:
    lines = txt_path.read_text(encoding="utf-8", errors="replace").splitlines()
    mods: Dict[str, Module] = {}

    i = 0
    while i < len(lines):
        line = lines[i].rstrip("\n")
        if line.startswith(">>") and not line.startswith(">>END_MODULE"):
            parts = line[2:].split("\t")
            name = parts[0].strip()
            status = parts[1].strip() if len(parts) > 1 else "unknown"
            i += 1

            header: List[str] = []
            cols: List[str] = []
            rows: List[List[str]] = []

            while i < len(lines) and not lines[i].startswith(">>END_MODULE"):
                if lines[i].startswith("#"):
                    header.append(lines[i])
                    cols = lines[i][1:].split("\t")  # last header wins
                else:
                    if lines[i].strip():
                        rows.append(lines[i].split("\t"))
                i += 1

            # normalize row length
            if cols:
                norm = []
                for r in rows:
                    if len(r) < len(cols):
                        r = r + [""] * (len(cols) - len(r))
                    elif len(r) > len(cols):
                        r = r[:len(cols)]
                    norm.append(r)
                rows = norm

            mods[name] = Module(name=name, status=status, header=header, rows=rows, cols=[c.strip() for c in cols])

        i += 1

    return mods


def get_basic_stats(mods: Dict[str, Module]) -> Dict[str, str]:
    m = mods.get("Basic Statistics")
    out: Dict[str, str] = {}
    if not m:
        return out
    # Basic Statistics table is Measure/Value
    # rows: [["Filename","..."], ["Total Sequences","..."], ...]
    for r in m.rows:
        if len(r) >= 2:
            out[r[0].strip()] = r[1].strip()
    return out


def q30_proxy_from_per_seq_qual(mods: Dict[str, Module]) -> Optional[float]:
    """
    FastQC txt doesn't provide base-level Q30 ratio. We compute a read-level proxy:
    fraction of reads with mean quality >= 30 using "Per sequence quality scores".
    """
    m = mods.get("Per sequence quality scores")
    if not m or not m.cols:
        return None
    # expect columns: Quality, Count
    try:
        q_idx = m.cols.index("Quality")
        c_idx = m.cols.index("Count")
    except ValueError:
        return None

    total = 0.0
    q30 = 0.0
    for r in m.rows:
        try:
            q = float(r[q_idx])
            c = float(r[c_idx])
        except Exception:
            continue
        total += c
        if q >= 30:
            q30 += c

    if total <= 0:
        return None
    return q30 / total


def median_iqr_from_length_dist(mods: Dict[str, Module]) -> Tuple[Optional[float], Optional[Tuple[float, float]]]:
    """
    Strict quantiles from FastQC "Sequence Length Distribution" binned table.

    FastQC provides bins like:
      Length   Count
      35-250   12345
      250      67890

    We compute weighted quantiles using bin LOWER/UPPER bounds.
    For a bin "a-b", we assume reads are uniformly distributed within [a, b].
    Then quantile is interpolated inside the bin, guaranteeing results within the global range.
    """
    m = mods.get("Sequence Length Distribution")
    if not m or not m.cols:
        return None, None
    if "Length" not in m.cols or "Count" not in m.cols:
        return None, None

    L = m.cols.index("Length")
    C = m.cols.index("Count")

    # Parse bins as (low, high, count)
    bins: List[Tuple[float, float, float]] = []
    total = 0.0

    for r in m.rows:
        try:
            length_str = str(r[L]).strip()
            count = float(r[C])
        except Exception:
            continue
        if count <= 0:
            continue

        if "-" in length_str:
            a, b = length_str.split("-", 1)
            low = float(a)
            high = float(b)
        else:
            low = float(length_str)
            high = float(length_str)

        # ensure low<=high
        if high < low:
            low, high = high, low

        bins.append((low, high, count))
        total += count

    if total <= 0 or not bins:
        return None, None

    # Sort by low then high
    bins.sort(key=lambda t: (t[0], t[1]))

    def quantile(p: float) -> float:
        """Return p-quantile with interpolation inside the bin."""
        target = total * p
        acc = 0.0
        for low, high, count in bins:
            prev = acc
            acc += count
            if acc >= target:
                # target is inside this bin
                if count == 0:
                    return low
                frac = (target - prev) / count  # 0..1 within bin mass
                if low == high:
                    return low
                # Uniform interpolation within [low, high]
                return low + frac * (high - low)
        # fallback (numerical edge)
        return bins[-1][1]

    q1 = quantile(0.25)
    q2 = quantile(0.50)
    q3 = quantile(0.75)
    return q2, (q1, q3)

# --------- HTML helpers (FastQC-like look) ---------

FASTQC_LIKE_CSS = """
<style type="text/css">
/* FastQC-like layout with fixed header height */

:root {
  --header-height: 90px;   /* 统一控制顶部高度 */
  --footer-height: 40px;
}

@media screen {

  div.header {
    background: #EEE;
    padding: 10px 20px;
    font-size: 28px;
    font-weight: bold;
    position: fixed;
    top: 0;
    left: 0;
    width: 100%;
    height: var(--header-height);
    box-sizing: border-box;
    z-index: 2;
  }

  div.summary {
    width: 18em;
    position: fixed;
    top: var(--header-height);
    left: 0;
    margin: 1em;
    z-index: 1;
  }

  div.main {
    position: absolute;
    top: var(--header-height);
    bottom: var(--footer-height);
    left: 18em;
    right: 0;
    border-left: 1px solid #CCC;
    padding: 0 0 0 1em;
    background: #FFF;
    overflow: auto;
  }

  div.footer {
    background: #EEE;
    padding: 0.5em;
    height: var(--footer-height);
    font-weight: bold;
    position: fixed;
    bottom: 0;
    left: 0;
    width: 100%;
    z-index: 2;
    overflow: hidden;
  }

  img.indented { margin-left: 3em; }
}

body { font-family: sans-serif; margin: 0; padding: 0; background: #FFF; color: #000; }
h2 { color: #800000; margin-bottom: 0; }
div.module { padding: 1.5em 0; }
table { margin-left: 3em; border-collapse: collapse; }
th { text-align: center; background-color: #000080; color: #FFF; padding: 0.4em; }
td { font-family: monospace; text-align: left; background-color: #EEE; color: #000; padding: 0.4em; }
a { color: #000080; }
a:hover { color: #800000; }

.badge { display:inline-block; padding: 2px 8px; border-radius: 10px; font-size: 12px; margin-right: 6px; }
.badge-pass { background:#1b5e20; color:#fff; }
.badge-warn { background:#f9a825; color:#000; }
.badge-fail { background:#b71c1c; color:#fff; }

.small { font-size: 12px; color: #444; margin-left: 3em; }

</style>
"""

def html_escape(s: str) -> str:
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
              .replace('"', "&quot;").replace("'", "&#39;"))


def status_badge(status: str) -> str:
    s = status.lower()
    if s.startswith("pass"):
        return '<span class="badge badge-pass">PASS</span>'
    if s.startswith("warn"):
        return '<span class="badge badge-warn">WARNING</span>'
    if s.startswith("fail"):
        return '<span class="badge badge-fail">FAIL</span>'
    return f'<span class="badge">{html_escape(status)}</span>'


def df_like_table(title: str, rows: List[Tuple[str, str]]) -> str:
    # FastQC basic stats table vibe
    out = [f'<h2>{html_escape(title)}</h2>']
    out.append('<table><tr><th>Measure</th><th>Value</th></tr>')
    for k, v in rows:
        out.append(f"<tr><td>{html_escape(k)}</td><td>{html_escape(v)}</td></tr>")
    out.append("</table>")
    return "\n".join(out)


def module_status_table(mods: Dict[str, Module], want: List[str]) -> str:
    out = ['<h2>Module Summary</h2>']
    out.append('<table><tr><th>Module</th><th>Status</th></tr>')
    for name in want:
        m = mods.get(name)
        if not m:
            out.append(f"<tr><td>{html_escape(name)}</td><td>{status_badge('missing')}</td></tr>")
        else:
            out.append(f"<tr><td>{html_escape(name)}</td><td>{status_badge(m.status)}</td></tr>")
    out.append("</table>")
    return "\n".join(out)


def embed_png_as_data_uri(png_path: Path) -> Optional[str]:
    if not png_path.exists():
        return None
    b64 = base64.b64encode(png_path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{b64}"


def module_plot_block(anchor: str, title: str, status: str, png_path: Path) -> str:
    uri = embed_png_as_data_uri(png_path)
    out = [f'<div class="module" id="{html_escape(anchor)}">']
    out.append(f"<h2>{status_badge(status)} {html_escape(title)}</h2>")
    if uri is None:
        out.append(f'<p class="small">Missing image: {html_escape(str(png_path.name))}</p>')
    else:
        out.append(f'<img class="indented" src="{uri}" alt="{html_escape(title)}"/>')
    out.append("</div>")
    return "\n".join(out)

def module_table_block(anchor: str, title: str, status: str, m: Optional[Module], empty_msg: str = "No data") -> str:
    out = [f'<div class="module" id="{html_escape(anchor)}">']
    out.append(f"<h2>{status_badge(status)} {html_escape(title)}</h2>")

    if (m is None) or (not m.cols) or (not m.rows):
        out.append(f'<p class="small">{html_escape(empty_msg)}</p>')
        out.append("</div>")
        return "\n".join(out)

    out.append("<table>")
    out.append("<tr>" + "".join(f"<th>{html_escape(c)}</th>" for c in m.cols) + "</tr>")
    for r in m.rows:
        out.append("<tr>" + "".join(f"<td>{html_escape(x)}</td>" for x in r) + "</tr>")
    out.append("</table>")
    out.append("</div>")
    return "\n".join(out)

# --------- Main ---------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastqc_txt", required=True, help="Path to FastQC fastqc_data.txt")
    ap.add_argument("--outdir", required=True, help="Output directory for HTML")
    ap.add_argument("--images_dir", default=None,
                    help="Path to FastQC Images directory (default: sibling 'Images' next to fastqc_data.txt)")
    args = ap.parse_args()

    txt_path = Path(args.fastqc_txt)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    images_dir = Path(args.images_dir) if args.images_dir else (txt_path.parent / "Images")

    mods = parse_fastqc_data_txt(txt_path)
    basic = get_basic_stats(mods)

    # Required metrics
    total_reads = basic.get("Total Sequences", "")
    total_bases = basic.get("Total Bases", "")
    length_range = basic.get("Sequence length", basic.get("Sequence Length", ""))
    poor = basic.get("Sequences flagged as poor quality", "")

    median_len, iqr = median_iqr_from_length_dist(mods)
    q30p = q30_proxy_from_per_seq_qual(mods)

    length_extra = []
    if median_len is not None:
        if iqr is not None:
            length_extra.append(f"median≈{median_len:.1f}, IQR≈({iqr[0]:.1f}, {iqr[1]:.1f})")
        else:
            length_extra.append(f"median≈{median_len:.1f}")

    q30_str = "" if q30p is None else f"{q30p:.4f}"

    # --- Build HTML ---
    title = f"Raw QC Dashboard - {basic.get('Filename','FastQC')}"
    header_html = f"""
<div class="header">
  <div id="header_title">Raw QC Dashboard</div>
  <div id="header_filename" style="float:right; font-size:50%; margin-right:2em; text-align:right;">
    {html_escape(basic.get("Filename",""))}<br/>
  </div>
</div>
"""

    # Left summary nav (FastQC-like)
    want_modules = [
        ("M1", "Per base sequence quality", "per_base_quality.png"),
        ("M2", "Per tile sequence quality", "per_tile_quality.png"),
        ("M3", "Per sequence quality scores", "per_sequence_quality.png"),
        ("M4", "Per base N content", "per_base_n_content.png"),
        ("M5", "Sequence Length Distribution", "sequence_length_distribution.png"),
        ("M6", "Adapter Content", "adapter_content.png"),
        ("M7", "Sequence Duplication Levels", "duplication_levels.png"),
    ]

    want_table_modules = [
        ("M8", "Overrepresented sequences"),
    ]

    summary_items = []
    summary_items.append('<li><a href="#M0">Basic Statistics</a></li>')
    summary_items.append('<li><a href="#M0b">Module Summary</a></li>')

    for anchor, name in want_table_modules:
        summary_items.append(f'<li><a href="#{html_escape(anchor)}">{html_escape(name)}</a></li>')
    for anchor, name, _png in want_modules:
        summary_items.append(f'<li><a href="#{html_escape(anchor)}">{html_escape(name)}</a></li>')

    summary_html = f"""
<div class="summary">
  <h2>Summary</h2>
  <ul>
    {''.join(summary_items)}
  </ul>
</div>
"""

    # Main content
    bs_rows = [
        ("Total Sequences (reads)", total_reads),
        ("Total Bases", total_bases),
        ("Sequence length (range)", length_range),
        ("Q30 ratio (proxy: reads mean Q≥30)", q30_str),
        ("Sequences flagged as poor quality", poor),
    ]
    basic_block = f'<div class="module" id="M0">{df_like_table("Basic Statistics", bs_rows)}</div>'

    status_block = f'<div class="module" id="M0b">{module_status_table(mods, [m[1] for m in want_modules] + [m[1] for m in want_table_modules])}</div>'
    note_block = """
<p class="small">
Note: FastQC fastqc_data.txt typically does not provide base-level Q30 bases ratio; we report a read-level proxy from "Per sequence quality scores".
</p>
"""

    plots = []
    # Image modules
    for anchor, mod_name, png in want_modules:
        status = mods.get(mod_name).status if mods.get(mod_name) else "missing"
        plots.append(module_plot_block(anchor, mod_name, status, images_dir / png))

    # Table modules (e.g., Overrepresented sequences)
    for anchor, mod_name in want_table_modules:
        m = mods.get(mod_name)
        status = m.status if m else "missing"
        plots.append(module_table_block(anchor, mod_name, status, m, empty_msg="No overrepresented sequences detected."))


    footer_html = '<div class="footer">Generated by raw_qc_report.py (FastQC-style)</div>'

    html = "\n".join([
        "<!DOCTYPE html><html><head>",
        f"<title>{html_escape(title)}</title>",
        FASTQC_LIKE_CSS,
        "</head><body>",
        header_html,
        summary_html,
        '<div class="main">',
        basic_block,
        note_block,
        status_block,
        "\n".join(plots),
        "</div>",
        footer_html,
        "</body></html>"
    ])

    out_html = outdir / "raw_qc_dashboard.html"
    out_html.write_text(html, encoding="utf-8")
    print(f"[OK] Wrote: {out_html}")


if __name__ == "__main__":
    main()
