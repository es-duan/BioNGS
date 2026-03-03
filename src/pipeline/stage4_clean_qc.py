from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple
import subprocess

from src.qc.overview_raw import compute_raw_overview, write_overview_report


def _run_cmd(cmd: List[str]) -> None:
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed ({p.returncode}): {' '.join(cmd)}\n{p.stdout}")


def _find_pop_fastq_pairs(umi_extracted_dir: Path) -> List[Tuple[str, Path, Path]]:
    """
    umi_extracted_dir = results/<exp>/umi_extracted/<run_tag>/
    期望：
      umi_extracted_dir/<pop>/<pop>_R1.fastq(.gz)
      umi_extracted_dir/<pop>/<pop>_R2.fastq(.gz)
      umi_extracted_dir/unmatched/unmatched_R1.fastq(.gz) (可选)
    """
    pairs: List[Tuple[str, Path, Path]] = []

    for pop_dir in sorted([p for p in umi_extracted_dir.iterdir() if p.is_dir()]):
        pop = pop_dir.name
        r1s = sorted(list(pop_dir.glob("*_R1.fastq")) + list(pop_dir.glob("*_R1.fastq.gz")))
        r2s = sorted(list(pop_dir.glob("*_R2.fastq")) + list(pop_dir.glob("*_R2.fastq.gz")))
        if not r1s or not r2s:
            continue
        pairs.append((pop, r1s[0], r2s[0]))

    if not pairs:
        raise FileNotFoundError(f"No FASTQ pairs found under: {umi_extracted_dir}")

    return pairs


def stage4_clean_qc(
    exp: str,
    run_tag: str,
    results_exp_dir: Path,
) -> Dict:
    """
    Stage 4: QC on "clean sequence" after stripping index+UMI+primers.
    - qc_details: fastqc per-pop
    - qc_overview: old-style overview (txt + length histogram) per-pop for R1/R2
    """
    umi_extracted_dir = results_exp_dir / "umi_extracted" / run_tag
    qc_details_root = results_exp_dir / "qc_details" / "02_clean_after_umi" / run_tag
    qc_overview_root = results_exp_dir / "qc_overview" / "02_clean_after_umi" / run_tag

    qc_details_root.mkdir(parents=True, exist_ok=True)
    qc_overview_root.mkdir(parents=True, exist_ok=True)

    pop_pairs = _find_pop_fastq_pairs(umi_extracted_dir)

    per_pop: Dict[str, Dict] = {}

    print("\n[Stage 4] Clean QC (after UMI/primer stripping) ...")

    for pop, r1, r2 in pop_pairs:
        # ---- overview (old style) ----
        out_over = qc_overview_root / pop
        (out_over / "R1").mkdir(parents=True, exist_ok=True)
        (out_over / "R2").mkdir(parents=True, exist_ok=True)

        m1 = compute_raw_overview(r1, short_len=150)  # 这里不用于 stop，只是复用长度统计/Q30统计能力
        m2 = compute_raw_overview(r2, short_len=150)

        write_overview_report(m1, out_over / "R1", title=f"{pop} - Clean R1 length distribution")
        write_overview_report(m2, out_over / "R2", title=f"{pop} - Clean R2 length distribution")

        # ---- details (fastqc) ----
        out_det = qc_details_root / pop
        out_det.mkdir(parents=True, exist_ok=True)
        _run_cmd(["fastqc", "-o", str(out_det), str(r1), str(r2)])

        per_pop[pop] = {
            "r1": str(r1),
            "r2": str(r2),
            "qc_overview_dir": str(out_over),
            "qc_details_dir": str(out_det),
            "r1_total_reads": m1["total_reads"],
            "r2_total_reads": m2["total_reads"],
            "r1_q30_rate": m1["q30_rate"],
            "r2_q30_rate": m2["q30_rate"],
        }

    print("[Stage 4] QC outputs written to:")
    print(f"  - {qc_overview_root}")
    print(f"  - {qc_details_root}")

    return {
        "status": "success",
        "umi_extracted_dir": str(umi_extracted_dir),
        "qc_overview_root": str(qc_overview_root),
        "qc_details_root": str(qc_details_root),
        "per_population": per_pop,
    }