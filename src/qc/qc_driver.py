# src/qc/qc_driver.py
from __future__ import annotations
import argparse
import subprocess
import sys
from pathlib import Path
from typing import List

from src.qc.io.paths import repo_root, results_dir
from src.qc.io.manifest import write_manifest
from src.qc.stages.stage_raw import build_raw_pairs
from src.qc.stages.stage_after_demux import build_demux_pairs


def run_entry_qc(entry_qc: Path, manifest: Path, outdir: Path, title: str) -> None:
    cmd = [
        sys.executable,
        str(entry_qc),
        "--manifest",
        str(manifest),
        "--outdir",
        str(outdir),
        "--title",
        title,
    ]
    print("\n[qc_driver] RUN:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    ap = argparse.ArgumentParser(
        description="Generate QC reports (RAW + After Script2 Demux). Future stages can be added without changing report renderer."
    )
    ap.add_argument("experiment", help="Experiment name, e.g. example")
    ap.add_argument("--gw_name", default=None, help="GW_name for raw fastq naming, e.g. P22R1 (used if raw paths not provided)")
    ap.add_argument("--raw_r1", default=None, help="Explicit raw R1 fastq path (optional)")
    ap.add_argument("--raw_r2", default=None, help="Explicit raw R2 fastq path (optional)")
    ap.add_argument("--raw_subdir", default="example_fastq", help="Subdir under input_data/{exp}/ for raw fastq (default: example_fastq)")
    ap.add_argument("--stages", default="raw,after_demux", help="Comma-separated stages (default: raw,after_demux)")
    args = ap.parse_args()

    exp = args.experiment
    stage_list: List[str] = [s.strip() for s in args.stages.split(",") if s.strip()]

    root = repo_root()
    entry_qc = root / "src" / "qc" / "entry_qc.py"
    if not entry_qc.exists():
        raise SystemExit(f"entry_qc.py not found: {entry_qc} (expected at src/entry_qc.py)")

    qc_base = results_dir(exp) / "qc"
    qc_base.mkdir(parents=True, exist_ok=True)

    # ---- Stage: RAW ----
    if "raw" in stage_list:
        outdir = qc_base / "00_raw"
        manifest = outdir / "manifest.csv"
        pairs = build_raw_pairs(exp, args.gw_name, args.raw_r1, args.raw_r2, raw_subdir=args.raw_subdir)
        write_manifest(pairs, manifest)
        run_entry_qc(entry_qc, manifest, outdir, f"Fancy QC (RAW) - {exp}")
        print("[qc_driver] RAW report:", outdir / "index.html")

    # ---- Stage: After Demux (Script2) ----
    if "after_demux" in stage_list:
        outdir = qc_base / "02_after_demux"
        manifest = outdir / "manifest.csv"
        pairs = build_demux_pairs(exp)
        write_manifest(pairs, manifest)
        run_entry_qc(entry_qc, manifest, outdir, f"Fancy QC (After Script2 Demux) - {exp}")
        print("[qc_driver] After-demux report:", outdir / "index.html")

    print("\n[qc_driver] DONE. Reports live under:")
    print(qc_base)


if __name__ == "__main__":
    main()