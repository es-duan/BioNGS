#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
entry_qc.py

Fancy QC runner that repeatedly calls fastqc_unpacker.py without changing report structure.

Supports:
  1) Single sample:
      python src/entry_qc.py --r1_path ... --r2_path ... --outdir ... --title ...
  2) Multi-sample via manifest (CSV/TSV with columns: label,r1,r2):
      python src/entry_qc.py --manifest manifest.csv --outdir ... --title ...
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path
from typing import List, Tuple, Optional


Pair = Tuple[str, Path, Path]  # (label, r1, r2)


def die(msg: str, code: int = 2) -> None:
    print(f"[entry_qc] ERROR: {msg}", file=sys.stderr)
    raise SystemExit(code)


def run_cmd(cmd: List[str]) -> None:
    print("[entry_qc] RUN:", " ".join(cmd))
    p = subprocess.run(cmd, text=True, capture_output=True)
    if p.returncode != 0:
        if p.stdout.strip():
            print("[entry_qc] STDOUT:\n" + p.stdout)
        if p.stderr.strip():
            print("[entry_qc] STDERR:\n" + p.stderr, file=sys.stderr)
        raise SystemExit(p.returncode)


def ensure_fastqc_available() -> None:
    if shutil.which("fastqc") is None:
        die("fastqc not found in PATH. Install FastQC and ensure `fastqc` is callable.")


def strip_fastq_ext(filename: str) -> str:
    """
    FastQC output is based on input filename, removing .gz and .fastq/.fq.
    """
    name = filename
    if name.endswith(".gz"):
        name = name[:-3]
    for ext in (".fastq", ".fq"):
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    return name


def expected_fastqc_zip_name(fq: Path) -> str:
    # <stem>_fastqc.zip
    return f"{strip_fastq_ext(fq.name)}_fastqc.zip"


def unzip_fastqc(zip_path: Path, extract_to: Path) -> Path:
    """
    Unzips zip to extract_to and returns extracted folder path (zip_path.stem).
    """
    if not zip_path.exists():
        die(f"FastQC zip not found: {zip_path}")

    extract_to.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "r") as z:
        z.extractall(extract_to)

    extracted_dir = extract_to / zip_path.stem
    if not extracted_dir.exists():
        die(f"Expected extracted dir not found after unzip: {extracted_dir}")
    return extracted_dir


def load_manifest(path: Path) -> List[Pair]:
    """
    Manifest CSV/TSV columns: label,r1,r2
    """
    if not path.exists():
        die(f"Manifest not found: {path}")

    pairs: List[Pair] = []
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        sample = f.read(4096)
        f.seek(0)
        dialect = csv.excel_tab if ("\t" in sample and "," not in sample) else csv.excel
        reader = csv.DictReader(f, dialect=dialect)

        if reader.fieldnames is None:
            die("Manifest has no header row.")
        fields = {h.strip() for h in reader.fieldnames}
        need = {"label", "r1", "r2"}
        if not need.issubset(fields):
            die(f"Manifest must contain columns: label,r1,r2. Got: {reader.fieldnames}")

        for row in reader:
            label = (row.get("label") or "").strip()
            r1 = Path((row.get("r1") or "").strip()).expanduser().resolve()
            r2 = Path((row.get("r2") or "").strip()).expanduser().resolve()
            if not label:
                die("Manifest row missing label.")
            if not r1.exists():
                die(f"R1 not found (label={label}): {r1}")
            if not r2.exists():
                die(f"R2 not found (label={label}): {r2}")
            pairs.append((label, r1, r2))

    if not pairs:
        die("Manifest has no data rows.")
    return pairs


def process_fastq(
    fq: Path,
    fastqc_out_dir: Path,
    unpacker_py: Path,
    report_out_dir: Path,
) -> Path:
    """
    Run fastqc on fq -> unzip -> call unpacker -> return path to raw_qc_dashboard.html
    """
    # 1) run FastQC
    fastqc_out_dir.mkdir(parents=True, exist_ok=True)
    run_cmd(["fastqc", str(fq), "-o", str(fastqc_out_dir)])

    # 2) unzip
    zip_path = fastqc_out_dir / expected_fastqc_zip_name(fq)
    extracted = unzip_fastqc(zip_path, fastqc_out_dir)

    fastqc_txt = extracted / "fastqc_data.txt"
    images_dir = extracted / "Images"
    if not fastqc_txt.exists():
        die(f"fastqc_data.txt not found: {fastqc_txt}")
    if not images_dir.exists():
        die(f"FastQC Images dir not found: {images_dir}")

    # 3) run unpacker (report structure unchanged)
    report_out_dir.mkdir(parents=True, exist_ok=True)
    run_cmd(
        [
            sys.executable,
            str(unpacker_py),
            "--fastqc_txt",
            str(fastqc_txt),
            "--outdir",
            str(report_out_dir),
            "--images_dir",
            str(images_dir),
        ]
    )

    html = report_out_dir / "raw_qc_dashboard.html"
    if not html.exists():
        die(f"Expected report not found: {html}")
    return html


def main() -> None:
    ap = argparse.ArgumentParser(description="Run fancy QC reports (FastQC + unpacker) for raw or per-pop fastqs.")
    ap.add_argument("--outdir", required=True, help="Output directory for fancy QC reports")
    ap.add_argument("--title", default="Fancy QC", help="Title for index.html")

    ap.add_argument("--manifest", default=None, help="CSV/TSV with columns: label,r1,r2 (multi-sample mode)")
    ap.add_argument("--r1_path", default=None, help="Full path to R1 fastq (single-sample mode)")
    ap.add_argument("--r2_path", default=None, help="Full path to R2 fastq (single-sample mode)")

    args = ap.parse_args()

    ensure_fastqc_available()

    # unpacker is expected next to this script: src/fastqc_unpacker.py
    this_dir = Path(__file__).resolve().parent
    repo_root = this_dir.parent.parent  # .../BioNGS
    unpacker_py = repo_root / "src" / "fastqc_unpacker.py"
    if not unpacker_py.exists():
        die(f"fastqc_unpacker.py not found next to entry_qc.py: {unpacker_py}")

    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # store fastqc outputs here (kept for audit/debug)
    fastqc_out_base = outdir / "_fastqc_out"
    fastqc_out_base.mkdir(parents=True, exist_ok=True)

    # decide pairs
    pairs: List[Pair]
    single_mode = False

    if args.manifest:
        pairs = load_manifest(Path(args.manifest).expanduser().resolve())
    else:
        if not (args.r1_path and args.r2_path):
            die("Provide either --manifest OR both --r1_path and --r2_path.")
        r1 = Path(args.r1_path).expanduser().resolve()
        r2 = Path(args.r2_path).expanduser().resolve()
        if not r1.exists():
            die(f"R1 not found: {r1}")
        if not r2.exists():
            die(f"R2 not found: {r2}")
        pairs = [("SAMPLE", r1, r2)]
        single_mode = True

    # generate reports
    index_items: List[str] = []

    for label, r1, r2 in pairs:
        if single_mode:
            # keep legacy structure: outdir/R1 and outdir/R2
            r1_report_dir = outdir / "R1"
            r2_report_dir = outdir / "R2"
            rel_r1 = "R1/raw_qc_dashboard.html"
            rel_r2 = "R2/raw_qc_dashboard.html"
            show = "Sample"
            fastqc_out_dir = fastqc_out_base
        else:
            # multi-sample structure: outdir/<label>/R1 and outdir/<label>/R2
            r1_report_dir = outdir / label / "R1"
            r2_report_dir = outdir / label / "R2"
            rel_r1 = f"{label}/R1/raw_qc_dashboard.html"
            rel_r2 = f"{label}/R2/raw_qc_dashboard.html"
            show = label
            fastqc_out_dir = fastqc_out_base / label

        process_fastq(r1, fastqc_out_dir, unpacker_py, r1_report_dir)
        process_fastq(r2, fastqc_out_dir, unpacker_py, r2_report_dir)

        index_items.append(
            f'<li><b>{show}</b>: <a href="{rel_r1}">R1</a> | <a href="{rel_r2}">R2</a></li>'
        )

    # write index.html
    index_html = outdir / "index.html"
    index_html.write_text(
        f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>{args.title}</title>
</head>
<body style="font-family:sans-serif; margin:24px;">
  <h1>{args.title}</h1>
  <ul>
    {chr(10).join(index_items)}
  </ul>
</body>
</html>
""",
        encoding="utf-8",
    )

    print(f"[entry_qc] OK: {index_html}")


if __name__ == "__main__":
    main()