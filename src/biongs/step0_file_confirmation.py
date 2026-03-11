"""Step 0: confirm input files and report runnable pipeline stage."""

from __future__ import annotations

import argparse
import difflib
import shutil
from dataclasses import dataclass
from pathlib import Path


@dataclass
class PresenceSummary:
    """Presence/absence status for required Step 0 inputs."""

    experiment_dir_exists: bool
    fastq_dir_exists: bool
    has_fastq_pairs: bool
    multiplexing_csv_exists: bool
    umi_csv_exists: bool
    fastqc_installed: bool


def _find_fastq_pairs(fastq_dir: Path) -> list[tuple[Path, Path]]:
    """Return matched R1/R2 fastq pairs from a fastq directory."""
    if not fastq_dir.exists() or not fastq_dir.is_dir():
        return []

    r1_files = sorted(
        [
            path
            for path in fastq_dir.iterdir()
            if path.is_file()
            and path.suffix.lower() in {".fastq", ".fq"}
            and "_R1" in path.name
        ]
    )

    pairs: list[tuple[Path, Path]] = []
    for r1_file in r1_files:
        paired_name = r1_file.name.replace("_R1", "_R2", 1)
        paired_path = r1_file.with_name(paired_name)
        if paired_path.exists() and paired_path.is_file():
            pairs.append((r1_file, paired_path))

    return pairs


def collect_presence(experiment_name: str, repo_root: Path) -> tuple[PresenceSummary, list[tuple[Path, Path]]]:
    """Collect required-file presence information for an experiment."""
    experiment_dir = repo_root / "input_data" / experiment_name
    fastq_dir = experiment_dir / f"{experiment_name}_fastq"
    multiplexing_csv = experiment_dir / f"{experiment_name}_multiplexing_info.csv"
    umi_csv = experiment_dir / f"{experiment_name}_UMI_primers.csv"

    pairs = _find_fastq_pairs(fastq_dir)
    summary = PresenceSummary(
        experiment_dir_exists=experiment_dir.exists() and experiment_dir.is_dir(),
        fastq_dir_exists=fastq_dir.exists() and fastq_dir.is_dir(),
        has_fastq_pairs=len(pairs) > 0,
        multiplexing_csv_exists=multiplexing_csv.exists() and multiplexing_csv.is_file(),
        umi_csv_exists=umi_csv.exists() and umi_csv.is_file(),
        fastqc_installed=shutil.which("fastqc") is not None,
    )
    return summary, pairs


def determine_furthest_script(summary: PresenceSummary) -> str:
    """Determine furthest pipeline script that can run from uploaded inputs."""
    if not summary.experiment_dir_exists:
        return "No pipeline scripts can run (experiment folder missing)."

    if not summary.has_fastq_pairs:
        return "Only Step 0 can run (missing at least one R1/R2 fastq pair)."

    if summary.has_fastq_pairs and not summary.multiplexing_csv_exists:
        return "Step 0a can run (FastQC on raw reads)."

    if summary.has_fastq_pairs and summary.multiplexing_csv_exists and not summary.umi_csv_exists:
        return "Step 3 can run (through indexing quality checks)."

    return "Step 6 can run (all required upload inputs are present)."


def find_possible_typos(experiment_name: str, repo_root: Path) -> list[tuple[str, str]]:
    """Find similarly named files/folders that may be typos of required names."""
    experiment_dir = repo_root / "input_data" / experiment_name
    if not experiment_dir.exists() or not experiment_dir.is_dir():
        return []

    expected_names = [
        f"{experiment_name}_fastq",
        f"{experiment_name}_multiplexing_info.csv",
        f"{experiment_name}_UMI_primers.csv",
    ]

    suggestions: list[tuple[str, str]] = []
    for item in sorted(experiment_dir.iterdir(), key=lambda path: path.name.lower()):
        if item.name in expected_names:
            continue

        if item.is_file() and item.suffix.lower() not in {".csv", ".fastq", ".fq"}:
            continue

        close = difflib.get_close_matches(item.name, expected_names, n=1, cutoff=0.6)
        if close:
            suggestions.append((item.name, close[0]))

    return suggestions


def build_report(
    experiment_name: str,
    summary: PresenceSummary,
    fastq_pairs: list[tuple[Path, Path]],
    typo_suggestions: list[tuple[str, str]],
) -> str:
    """Build text report for terminal and file output."""
    required_status = {
        f"input_data/{experiment_name}/": "present" if summary.experiment_dir_exists else "missing",
        f"input_data/{experiment_name}/{experiment_name}_fastq/": "present" if summary.fastq_dir_exists else "missing",
        "At least one R1/R2 fastq pair": "present" if summary.has_fastq_pairs else "missing",
        f"input_data/{experiment_name}/{experiment_name}_multiplexing_info.csv": "present"
        if summary.multiplexing_csv_exists
        else "missing",
        f"input_data/{experiment_name}/{experiment_name}_UMI_primers.csv": "present"
        if summary.umi_csv_exists
        else "missing",
        "FastQC (required for Step 0a)": "installed" if summary.fastqc_installed else "NOT INSTALLED",
    }

    lines = [
        "=" * 72,
        f"Step 0 Report: File Confirmation for experiment '{experiment_name}'",
        "=" * 72,
        "",
        "Required inputs:",
    ]

    for label, status in required_status.items():
        lines.append(f"- {label}: {status}")

    lines.extend(["", f"Detected fastq pairs: {len(fastq_pairs)}"])
    for r1_file, r2_file in fastq_pairs:
        lines.append(f"  - {r1_file.name} <-> {r2_file.name}")

    lines.extend(["", "Furthest runnable stage:", f"- {determine_furthest_script(summary)}", ""])

    if not summary.fastqc_installed:
        lines.append(
            "Warning (Step 0a): FastQC is not installed and Step 0a cannot run. "
            "See https://github.com/s-andrews/FastQC/blob/master/INSTALL.md for instructions, "
            "then re-run Step 0 to confirm."
        )
        lines.append("")

    if typo_suggestions:
        lines.append("Possible naming typos:")
        for observed, expected in typo_suggestions:
            lines.append(f"- '{observed}' looks similar to '{expected}'")
    else:
        lines.append("Possible naming typos:")
        lines.append("- none detected")

    lines.append("")
    return "\n".join(lines)


def write_report(experiment_name: str, report: str, repo_root: Path) -> Path:
    """Write Step 0 report text file to the experiment logs folder."""
    report_dir = repo_root / "results" / experiment_name / "logs"
    report_dir.mkdir(parents=True, exist_ok=True)
    report_path = report_dir / "file_confirmation_report.txt"
    report_path.write_text(report, encoding="utf-8")
    return report_path


def main() -> None:
    """CLI entry point for Step 0."""
    parser = argparse.ArgumentParser(
        description="Check required upload files and report the furthest runnable pipeline stage.",
    )
    parser.add_argument("experiment_name", help="Experiment folder name under input_data/ (e.g., example)")
    args = parser.parse_args()

    repo_root = Path.cwd()
    summary, fastq_pairs = collect_presence(args.experiment_name, repo_root)
    typo_suggestions = find_possible_typos(args.experiment_name, repo_root)
    report = build_report(args.experiment_name, summary, fastq_pairs, typo_suggestions)

    print(report)
    report_path = write_report(args.experiment_name, report, repo_root)
    print(f"Report saved to: {report_path}")


if __name__ == "__main__":
    main()
