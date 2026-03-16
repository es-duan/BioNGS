from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Dict, Any


class AlignmentQCError(RuntimeError):
    pass


def _run_cmd(cmd: list[str]) -> str:
    """Run command and return stdout"""
    print("Running:", " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)

    if res.returncode != 0:
        print(res.stderr)
        raise AlignmentQCError("Command failed")

    return res.stdout


def _load_stage8_metrics(results_exp_dir: Path) -> Dict[str, Any]:
    metrics_dir = results_exp_dir / "metrics"
    files = sorted(metrics_dir.glob("metrics_stage8_alignment_*.json"))

    if not files:
        raise FileNotFoundError("Stage8 metrics file not found")

    return json.loads(files[-1].read_text())


def _infer_version_from_metrics(metrics: Dict[str, Any]) -> str:
    """
    Prefer an explicit version field from metrics.
    Fall back to inferring version from BAM path:
      results/<exp>/alignment/<run_tag>/<version>/<pop>/<pop>.bam
    """
    version = metrics.get("version") or metrics.get("trimmed_version")
    if version:
        return str(version)

    per_population = metrics.get("per_population", {})
    if not per_population:
        raise KeyError("No per_population field found in Stage8 metrics")

    first_pop_info = next(iter(per_population.values()))
    bam_str = first_pop_info.get("bam")
    if not bam_str:
        raise KeyError("No bam field found in Stage8 metrics per_population entry")

    bam_path = Path(bam_str)

    # 期望路径中包含 alignment/<run_tag>/<version>/<pop>/<file>.bam
    parts = bam_path.parts
    if "alignment" in parts:
        idx = parts.index("alignment")
        if idx + 2 < len(parts):
            return parts[idx + 2]

    raise KeyError(
        "Could not determine version from Stage8 metrics. "
        "Expected metrics['version'] or a BAM path containing alignment/<run_tag>/<version>/..."
    )


def _parse_flagstat(text: str) -> Dict[str, Any]:
    total = None
    mapped = None
    paired = None

    for line in text.splitlines():
        if "in total" in line:
            total = int(line.split()[0])

        if "mapped (" in line:
            mapped = int(line.split()[0])

        if "properly paired" in line:
            paired = int(line.split()[0])

    mapping_rate = mapped / total if total else 0

    return {
        "total_reads": total,
        "mapped_reads": mapped,
        "properly_paired": paired,
        "mapping_rate": mapping_rate,
    }


def _parse_coverage(text: str) -> Dict[str, Any]:
    lines = text.splitlines()

    if len(lines) < 2:
        return {}

    fields = lines[1].split()

    return {
        "coverage_percent": float(fields[5]),
        "mean_depth": float(fields[6]),
    }


def _parse_depth(text: str) -> Dict[str, Any]:
    depths = []

    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 3:
            depths.append(int(parts[2]))

    if not depths:
        return {}

    depths.sort()
    n = len(depths)

    mean_depth = sum(depths) / n
    median_depth = depths[n // 2]
    max_depth = max(depths)

    return {
        "mean_depth": mean_depth,
        "median_depth": median_depth,
        "max_depth": max_depth,
    }


def stage8_5_alignment_qc(
    exp: str,
    run_tag: str,
    results_exp_dir: Path,
) -> Dict[str, Any]:

    metrics = _load_stage8_metrics(results_exp_dir)
    version = _infer_version_from_metrics(metrics)

    qc_root = results_exp_dir / "alignment_qc" / run_tag / version
    qc_root.mkdir(parents=True, exist_ok=True)

    payload = {
        "stage": "stage8_5_alignment_qc",
        "exp": exp,
        "run_tag": run_tag,
        "version": version,
        "per_population": {},
    }

    for pop, info in metrics["per_population"].items():
        bam = Path(info["bam"])

        pop_dir = qc_root / pop
        pop_dir.mkdir(parents=True, exist_ok=True)

        print("\nAlignment QC for population", pop)
        print("Input BAM:", bam)

        flagstat_txt = pop_dir / "flagstat.txt"
        coverage_txt = pop_dir / "coverage.txt"
        depth_txt = pop_dir / "depth.txt"

        flagstat = _run_cmd(["samtools", "flagstat", str(bam)])
        flagstat_txt.write_text(flagstat, encoding="utf-8")

        coverage = _run_cmd(["samtools", "coverage", str(bam)])
        coverage_txt.write_text(coverage, encoding="utf-8")

        depth = _run_cmd(["samtools", "depth", str(bam)])
        depth_txt.write_text(depth, encoding="utf-8")

        flagstat_stats = _parse_flagstat(flagstat)
        coverage_stats = _parse_coverage(coverage)
        depth_stats = _parse_depth(depth)

        summary = {
            **flagstat_stats,
            **coverage_stats,
            **depth_stats,
        }

        qc_status = "PASS"
        if summary.get("mapping_rate", 0) < 0.8:
            qc_status = "WARNING"
        if summary.get("mean_depth", 0) < 30:
            qc_status = "WARNING"

        summary["qc_status"] = qc_status
        payload["per_population"][pop] = {
            **summary,
            "bam": str(bam),
            "flagstat_txt": str(flagstat_txt),
            "coverage_txt": str(coverage_txt),
            "depth_txt": str(depth_txt),
        }

        print("\nQC summary")
        for k, v in summary.items():
            print(f"{k:20s} {v}")

    out_json = results_exp_dir / "metrics" / f"metrics_stage8_5_alignment_qc_{version}.json"
    out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print("\nAlignment QC results written to:")
    print(qc_root)
    print("Metrics written to:")
    print(out_json)

    return payload