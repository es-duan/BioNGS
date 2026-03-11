"""Step 0a: run FastQC on raw reads and keep HTML reports only."""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path


def _strip_fastq_ext(filename: str) -> str:
    """Match FastQC naming by removing .gz then .fastq/.fq extension."""
    name = filename
    if name.endswith(".gz"):
        name = name[:-3]
    for ext in (".fastq", ".fq"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def find_input_fastqs(experiment_name: str, repo_root: Path) -> list[Path]:
    """Find all FASTQ inputs in input_data/{experiment}/{experiment}_fastq."""
    fastq_dir = repo_root / "input_data" / experiment_name / f"{experiment_name}_fastq"
    if not fastq_dir.exists() or not fastq_dir.is_dir():
        raise FileNotFoundError(
            f"FASTQ directory not found: {fastq_dir}. Run Step 0 to confirm upload structure."
        )

    candidates = sorted(
        [
            path
            for path in fastq_dir.iterdir()
            if path.is_file()
            and (
                path.name.endswith(".fastq")
                or path.name.endswith(".fq")
                or path.name.endswith(".fastq.gz")
                or path.name.endswith(".fq.gz")
            )
        ],
        key=lambda path: path.name,
    )

    if not candidates:
        raise FileNotFoundError(f"No FASTQ files found in: {fastq_dir}")

    return candidates


def run_fastqc(fastq_files: list[Path], outdir: Path) -> None:
    """Execute FastQC on all input files."""
    if shutil.which("fastqc") is None:
        raise RuntimeError("fastqc is not available in PATH. Install FastQC before running Step 0a.")

    outdir.mkdir(parents=True, exist_ok=True)
    cmd = ["fastqc", "--quiet", "--outdir", str(outdir)] + [str(path) for path in fastq_files]
    subprocess.run(cmd, check=True)


def _expected_html_names(fastq_files: list[Path]) -> set[str]:
    """Return expected FastQC html filenames for provided FASTQ files."""
    return {f"{_strip_fastq_ext(path.name)}_fastqc.html" for path in fastq_files}


def cleanup_non_html_fastqc_outputs(outdir: Path, expected_html: set[str]) -> None:
    """Remove zip/extracted outputs and any unexpected files from the QC folder."""
    for path in outdir.iterdir():
        if path.is_dir() and path.name.endswith("_fastqc"):
            shutil.rmtree(path)
            continue

        if path.is_file() and path.name.endswith("_fastqc.zip"):
            path.unlink()
            continue

        if path.is_file() and path.suffix.lower() == ".html" and path.name in expected_html:
            continue

        if path.is_file() and path.name == ".DS_Store":
            continue

        if path.is_file():
            path.unlink()


def run_script_0_5(experiment_name: str, repo_root: Path) -> list[Path]:
    """Run FastQC and return generated HTML report paths."""
    fastq_files = find_input_fastqs(experiment_name, repo_root)
    outdir = repo_root / "results" / experiment_name / "qc"
    run_fastqc(fastq_files, outdir)

    expected_html = _expected_html_names(fastq_files)
    cleanup_non_html_fastqc_outputs(outdir, expected_html)

    html_paths = sorted([outdir / name for name in expected_html if (outdir / name).exists()])
    if not html_paths:
        raise RuntimeError("FastQC completed but no HTML reports were found.")

    return html_paths


def main() -> None:
    """CLI entry point for Step 0a."""
    parser = argparse.ArgumentParser(
        description="Run FastQC for raw FASTQ files and keep only HTML outputs.",
    )
    parser.add_argument("experiment_name", help="Experiment folder name under input_data/ (e.g., example)")
    args = parser.parse_args()

    html_paths = run_script_0_5(args.experiment_name, Path.cwd())

    print(f"Generated {len(html_paths)} FastQC HTML report(s):")
    for html_path in html_paths:
        print(f"- {html_path}")


if __name__ == "__main__":
    main()
