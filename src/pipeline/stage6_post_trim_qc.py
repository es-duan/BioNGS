from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple
import subprocess

from src.qc.overview_post_trim import (
    compute_post_trim_overview,
    aggregate_post_trim_metrics,
    write_post_trim_overview,
)


def _run_cmd(cmd: List[str]) -> None:
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed ({p.returncode}): {' '.join(cmd)}\n{p.stdout}")


def _find_population_pairs(trimmed_dir: Path) -> List[Tuple[str, Path, Path]]:
    """
    trimmed_dir expected:
      trimmed_dir/<pop>/*R1*.fastq(.gz)
      trimmed_dir/<pop>/*R2*.fastq(.gz)
    Returns list of (pop, r1, r2).
    """
    pairs: List[Tuple[str, Path, Path]] = []

    for pop_dir in sorted([p for p in trimmed_dir.iterdir() if p.is_dir()]):
        pop = pop_dir.name
        r1s = sorted(list(pop_dir.glob("*R1*.fastq")) + list(pop_dir.glob("*R1*.fastq.gz")))
        r2s = sorted(list(pop_dir.glob("*R2*.fastq")) + list(pop_dir.glob("*R2*.fastq.gz")))
        if not r1s or not r2s:
            continue
        pairs.append((pop, r1s[0], r2s[0]))

    if not pairs:
        raise FileNotFoundError(f"No trimmed FASTQ pairs found under: {trimmed_dir}")

    return pairs


def stage6_post_trim_qc(
    exp: str,
    run_tag: str,
    trimmed_dir: Path,
    results_exp_dir: Path,
) -> Dict:
    """
    Step 6) Post-preprocess QC
    - qc_details: FastQC per population (照旧)
    - qc_overview: per-pop summary txt + read length frequency plot (NO bar chart)
    - qc_overview overall: aggregate across all populations (R1 overall + R2 overall)

    Outputs:
      results/<exp>/qc_details/01_post_trim/<run_tag>/<pop>/
      results/<exp>/qc_overview/01_post_trim/<run_tag>/<pop>/R1|R2/
      results/<exp>/qc_overview/01_post_trim/<run_tag>/overall/R1|R2/
    """
    qc_details_root = results_exp_dir / "qc_details" / "01_post_trim" / run_tag
    qc_overview_root = results_exp_dir / "qc_overview" / "01_post_trim" / run_tag
    overall_root = qc_overview_root / "overall"

    qc_details_root.mkdir(parents=True, exist_ok=True)
    qc_overview_root.mkdir(parents=True, exist_ok=True)
    overall_root.mkdir(parents=True, exist_ok=True)

    pop_pairs = _find_population_pairs(trimmed_dir)

    # 收集 per-pop metrics，用于 overall 聚合
    all_r1_metrics: List[dict] = []
    all_r2_metrics: List[dict] = []

    # 1) per-pop overview
    overview_written: Dict[str, Dict[str, str]] = {}
    for pop, r1, r2 in pop_pairs:
        out_pop = qc_overview_root / pop
        out_pop.mkdir(parents=True, exist_ok=True)

        m1 = compute_post_trim_overview(r1)
        m2 = compute_post_trim_overview(r2)

        all_r1_metrics.append(m1)
        all_r2_metrics.append(m2)

        write_post_trim_overview(m1, out_pop / "R1", title=f"{exp} | {run_tag} | {pop} | Post-trim QC | R1")
        write_post_trim_overview(m2, out_pop / "R2", title=f"{exp} | {run_tag} | {pop} | Post-trim QC | R2")

        overview_written[pop] = {
            "r1_dir": str(out_pop / "R1"),
            "r2_dir": str(out_pop / "R2"),
        }

    # 2) overall overview (aggregate)
    overall_r1 = aggregate_post_trim_metrics(all_r1_metrics)
    overall_r2 = aggregate_post_trim_metrics(all_r2_metrics)

    write_post_trim_overview(
        overall_r1,
        overall_root / "R1",
        title=f"{exp} | {run_tag} | OVERALL | Post-trim QC | R1"
    )
    write_post_trim_overview(
        overall_r2,
        overall_root / "R2",
        title=f"{exp} | {run_tag} | OVERALL | Post-trim QC | R2"
    )

    # 3) qc_details: FastQC per pop
    for pop, r1, r2 in pop_pairs:
        outd = qc_details_root / pop
        outd.mkdir(parents=True, exist_ok=True)
        _run_cmd(["fastqc", "-o", str(outd), str(r1), str(r2)])

    print("\n[Stage 6] Post-trim QC generated.")
    print(f"[Stage 6] qc_details:  {qc_details_root}")
    print(f"[Stage 6] qc_overview: {qc_overview_root}")
    print(f"[Stage 6] overall:     {overall_root}")

    return {
        "status": "success",
        "trimmed_dir": str(trimmed_dir),
        "qc_details_dir": str(qc_details_root),
        "qc_overview_dir": str(qc_overview_root),
        "overall_overview_dir": str(overall_root),
        "per_population_overview": overview_written,
        "n_populations": len(pop_pairs),
        "overall_q30_rate_r1": overall_r1["q30_rate"],
        "overall_q30_rate_r2": overall_r2["q30_rate"],
    }