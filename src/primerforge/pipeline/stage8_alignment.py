from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Dict, Any


class AlignmentError(RuntimeError):
    pass


def _run_cmd(cmd: list[str]) -> None:
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise AlignmentError("Command failed")


def _load_stage7_metrics(results_exp_dir: Path) -> Dict[str, Any]:
    metrics_dir = results_exp_dir / "metrics"
    files = sorted(metrics_dir.glob("metrics_stage7_family_strategy_*.json"))

    if not files:
        raise FileNotFoundError("Stage7 metrics file not found")

    return json.loads(files[-1].read_text())


def _bowtie2_align(
    r1: Path,
    r2: Path,
    reference_index: Path,
    out_bam: Path,
):
    sam = out_bam.with_suffix(".sam")

    cmd = [
        "bowtie2",
        "-x", str(reference_index),
        "-1", str(r1),
        "-2", str(r2),
        "-S", str(sam),
    ]

    _run_cmd(cmd)

    sort_cmd = [
        "samtools",
        "sort",
        "-o",
        str(out_bam),
        str(sam),
    ]

    _run_cmd(sort_cmd)

    index_cmd = [
        "samtools",
        "index",
        str(out_bam),
    ]

    _run_cmd(index_cmd)

    sam.unlink()


def stage8_alignment(
    exp: str,
    run_tag: str,
    results_exp_dir: Path,
    reference_index: Path,
) -> Dict[str, Any]:

    metrics = _load_stage7_metrics(results_exp_dir)

    version = metrics["trimmed_version"]

    out_root = results_exp_dir / "alignment" / run_tag / version
    out_root.mkdir(parents=True, exist_ok=True)

    payload = {
        "stage": "stage8_alignment",
        "per_population": {},
    }

    for pop, info in metrics["per_population"].items():

        strategy_info = info["strategy_output"]

        r1 = Path(strategy_info["out_r1"])
        r2 = Path(strategy_info["out_r2"])

        pop_dir = out_root / pop
        pop_dir.mkdir(parents=True, exist_ok=True)

        bam = pop_dir / f"{pop}.bam"

        print(f"\nAligning population {pop}")
        print("Input R1:", r1)
        print("Input R2:", r2)

        _bowtie2_align(r1, r2, reference_index, bam)

        payload["per_population"][pop] = {
            "bam": str(bam),
            "bai": str(bam) + ".bai",
        }

    out_json = results_exp_dir / "metrics" / f"metrics_stage8_alignment_{version}.json"
    out_json.write_text(json.dumps(payload, indent=2))

    print("\nStage8 alignment finished")
    print("Metrics written to:", out_json)

    return payload