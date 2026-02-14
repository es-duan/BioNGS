#!/usr/bin/env python3
# entry_qc.py
# Entry in the SAME folder as fastqc_unpacker.py
# Input: R1/R2 FASTQ in input_data/example/example_fastq
# Steps: FastQC -> unzip -> call fastqc_unpacker.py -> produce HTML

from __future__ import annotations
import argparse
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path


def run(cmd: list[str]) -> None:
    print("[RUN]", " ".join(cmd))
    p = subprocess.run(cmd, text=True, capture_output=True)
    if p.returncode != 0:
        print("[STDOUT]\n", p.stdout)
        print("[STDERR]\n", p.stderr, file=sys.stderr)
        raise SystemExit(p.returncode)
    if p.stdout.strip():
        print("[STDOUT]\n", p.stdout)


def ensure_fastqc() -> None:
    if shutil.which("fastqc") is None:
        raise SystemExit("ERROR: fastqc not found in PATH. Try: fastqc --version")


def strip_fastq_ext(name: str) -> str:
    for ext in [".fastq.gz", ".fastq", ".fq.gz", ".fq"]:
        if name.endswith(ext):
            return name[:-len(ext)]
    return Path(name).stem


def fastqc_zip_name(fq: Path) -> str:
    # FastQC zip is <stem>_fastqc.zip
    return f"{strip_fastq_ext(fq.name)}_fastqc.zip"


def unzip_fastqc(zip_path: Path, outdir: Path) -> Path:
    # Returns extracted folder path <zip_stem> (e.g., P22R1_R1_001_fastqc/)
    if not zip_path.exists():
        raise SystemExit(f"ERROR: FastQC zip not found: {zip_path}")

    with zipfile.ZipFile(zip_path, "r") as z:
        z.extractall(outdir)

    extracted = outdir / zip_path.stem
    if not extracted.exists():
        raise SystemExit(f"ERROR: unzip done but folder not found: {extracted}")
    return extracted


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--input_dir",
        default="/home/uniscan/BioNGS/input_data/example/example_fastq",
        help="Directory containing R1/R2 fastq files",
    )
    ap.add_argument(
        "--outdir",
        default="/home/uniscan/BioNGS/draft/results/raw_qc_example",
        help="Where to write FastQC outputs + custom HTML reports",
    )
    ap.add_argument(
        "--r1",
        default="P22R1_R1_001.fastq",
        help="R1 fastq filename inside input_dir",
    )
    ap.add_argument(
        "--r2",
        default="P22R1_R2_001.fastq",
        help="R2 fastq filename inside input_dir",
    )
    args = ap.parse_args()

    ensure_fastqc()

    here = Path(__file__).resolve().parent
    imagereader = (here / "fastqc_unpacker.py").resolve()
    if not imagereader.exists():
        raise SystemExit(f"ERROR: fastqc_unpacker.py not found next to entry script: {imagereader}")

    input_dir = Path(args.input_dir).resolve()
    if not input_dir.exists():
        raise SystemExit(f"ERROR: input_dir not found: {input_dir}")

    r1 = (input_dir / args.r1).resolve()
    r2 = (input_dir / args.r2).resolve()
    if not r1.exists():
        raise SystemExit(f"ERROR: R1 not found: {r1}")
    if not r2.exists():
        raise SystemExit(f"ERROR: R2 not found: {r2}")

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Put fastqc outputs in one place
    fastqc_out = outdir / "fastqc_out"
    fastqc_out.mkdir(parents=True, exist_ok=True)

    def process_one(label: str, fq: Path) -> Path:
        # 1) FastQC
        run(["fastqc", str(fq), "-o", str(fastqc_out)])

        # 2) unzip
        zip_path = fastqc_out / fastqc_zip_name(fq)
        extracted_dir = unzip_fastqc(zip_path, fastqc_out)

        fastqc_txt = extracted_dir / "fastqc_data.txt"
        if not fastqc_txt.exists():
            raise SystemExit(f"ERROR: fastqc_data.txt not found: {fastqc_txt}")

        # 3) call imagereader to create dashboard
        report_out = outdir / label
        report_out.mkdir(parents=True, exist_ok=True)

        run([
            sys.executable, str(imagereader),
            "--fastqc_txt", str(fastqc_txt),
            "--outdir", str(report_out),
        ])

        html = report_out / "raw_qc_dashboard.html"
        if not html.exists():
            raise SystemExit(f"ERROR: expected HTML not found: {html}")
        print(f"[OK] {label} HTML: {html}")
        return html

    r1_html = process_one("R1", r1)
    r2_html = process_one("R2", r2)

    # Write a simple combined index.html
    index = outdir / "index.html"
    index.write_text(
        f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Raw QC Example</title></head>
<body style="font-family:sans-serif; margin:24px;">
<h1>Raw QC Example</h1>
<ul>
  <li><a href="R1/{r1_html.name}">R1 report</a></li>
  <li><a href="R2/{r2_html.name}">R2 report</a></li>
</ul>
</body></html>
""",
        encoding="utf-8",
    )
    print(f"[OK] index: {index}")
    print("[DONE]")


if __name__ == "__main__":
    main()
