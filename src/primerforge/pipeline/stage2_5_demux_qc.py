from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import subprocess

from primerforge.qc.overview_raw import compute_raw_overview, write_overview_report


def _run_cmd(cmd: List[str]) -> None:
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed ({p.returncode}): {' '.join(cmd)}\n{p.stdout}")


def _find_population_pairs(demux_dir: Path) -> List[Tuple[str, Path, Path]]:
    """
    demux_dir = results/<exp>/demultiplexing/<run_tag>/
    期望结构：
      demux_dir/P*/..._R1.fastq(.gz)
      demux_dir/P*/..._R2.fastq(.gz)
      demux_dir/unmatched_R1.fastq(.gz), unmatched_R2.fastq(.gz)
    返回: [(pop_name, r1_path, r2_path), ...] 包括 unmatched
    """
    pairs: List[Tuple[str, Path, Path]] = []

    # populations
    for pop_dir in sorted([p for p in demux_dir.iterdir() if p.is_dir()]):
        pop = pop_dir.name
        r1s = sorted(list(pop_dir.glob("*_R1.fastq")) + list(pop_dir.glob("*_R1.fastq.gz")))
        r2s = sorted(list(pop_dir.glob("*_R2.fastq")) + list(pop_dir.glob("*_R2.fastq.gz")))
        if not r1s or not r2s:
            continue
        pairs.append((pop, r1s[0], r2s[0]))

    # unmatched
    um1 = list(demux_dir.glob("unmatched_R1.fastq")) + list(demux_dir.glob("unmatched_R1.fastq.gz"))
    um2 = list(demux_dir.glob("unmatched_R2.fastq")) + list(demux_dir.glob("unmatched_R2.fastq.gz"))
    if um1 and um2:
        pairs.append(("unmatched", um1[0], um2[0]))

    if not pairs:
        raise FileNotFoundError(f"No demux FASTQ pairs found under: {demux_dir}")

    return pairs


def stage2_5_demux_qc(
    *,
    exp: str,
    run_tag: str,
    # ✅ 新接口：run_pipeline 传这个
    results_exp_dir: Optional[Path] = None,
    # ✅ 兼容旧接口：你也可以直接传这三个（用于单独测试）
    demux_dir: Optional[Path] = None,
    qc_details_root: Optional[Path] = None,
    qc_overview_root: Optional[Path] = None,
) -> Dict[str, Any]:
    """
    Stage 2.5：
      - qc_details: 对每个群跑 fastqc（照旧）
      - qc_overview: 对每个群用老代码 compute_raw_overview + write_overview_report（png+txt）
      - 额外：生成一个 demux_summary.txt + demux_distribution.png（总体分布，不报告 short）
    """

    # -------------------------
    # ✅ 参数归一化：优先用 results_exp_dir 推导路径
    # -------------------------
    if results_exp_dir is not None:
        results_exp_dir = Path(results_exp_dir)
        demux_dir = results_exp_dir / "demultiplexing" / run_tag
        qc_details_root = results_exp_dir / "qc_details"
        qc_overview_root = results_exp_dir / "qc_overview"
    else:
        # 走旧接口时，这三个必须齐
        if demux_dir is None or qc_details_root is None or qc_overview_root is None:
            raise ValueError(
                "stage2_5_demux_qc requires either results_exp_dir OR "
                "(demux_dir + qc_details_root + qc_overview_root)."
            )

    assert demux_dir is not None
    assert qc_details_root is not None
    assert qc_overview_root is not None

    # -------------------------
    # 输出目录
    # -------------------------
    qc_details_dir = qc_details_root / run_tag
    qc_overview_dir = qc_overview_root / run_tag
    per_pop_overview = qc_overview_dir / "per_population"

    qc_details_dir.mkdir(parents=True, exist_ok=True)
    per_pop_overview.mkdir(parents=True, exist_ok=True)

    pop_pairs = _find_population_pairs(demux_dir)

    # -------------------------
    # 统计 demux 分布（pairs 口径：用 R1 的 total_reads 作为 pairs count）
    # -------------------------
    pop_pair_counts: Dict[str, int] = {}
    for pop, r1, _r2 in pop_pairs:
        m = compute_raw_overview(r1, short_len=150)
        pop_pair_counts[pop] = m["total_reads"]

    total_pairs = sum(pop_pair_counts.values())
    unmatched_pairs = pop_pair_counts.get("unmatched", 0)
    unmatched_rate = unmatched_pairs / total_pairs if total_pairs else 0.0

    # -------------------------
    # 1) per-pop overview（老代码风）
    # -------------------------
    for pop, r1, r2 in pop_pairs:
        out_pop = per_pop_overview / pop
        out_pop.mkdir(parents=True, exist_ok=True)

        m1 = compute_raw_overview(r1, short_len=150)
        m2 = compute_raw_overview(r2, short_len=150)

        write_overview_report(m1, out_pop / "R1", title=f"{pop} - R1 Read Length Distribution")
        write_overview_report(m2, out_pop / "R2", title=f"{pop} - R2 Read Length Distribution")

    # -------------------------
    # 2) qc_details：fastqc per-pop（照旧）
    # -------------------------
    for pop, r1, r2 in pop_pairs:
        outd = qc_details_dir / pop
        outd.mkdir(parents=True, exist_ok=True)
        _run_cmd(["fastqc", "-o", str(outd), str(r1), str(r2)])

    # -------------------------
    # 3) demux_summary.txt + demux_distribution.png
    # -------------------------
    qc_overview_dir.mkdir(parents=True, exist_ok=True)
    summary_txt = qc_overview_dir / "demux_summary.txt"
    with summary_txt.open("w", encoding="utf-8") as f:
        f.write(f"Experiment: {exp}\n")
        f.write(f"Run tag: {run_tag}\n")
        f.write("=" * 60 + "\n")
        for pop in sorted(pop_pair_counts.keys()):
            f.write(f"{pop}\t{pop_pair_counts[pop]}\n")
        f.write("=" * 60 + "\n")
        f.write(f"total_pairs\t{total_pairs}\n")
        f.write(f"unmatched_pairs\t{unmatched_pairs}\n")
        f.write(f"unmatched_rate\t{unmatched_rate:.3%}\n")

    import matplotlib.pyplot as plt

    labels = [p for p in sorted(pop_pair_counts.keys())]
    values = [pop_pair_counts[p] for p in labels]

    plt.figure(figsize=(8, 5))
    plt.bar(labels, values, width=0.75, edgecolor="black", linewidth=0.4)
    plt.title("Demultiplex Read Pair Distribution", fontsize=12, pad=10)
    plt.ylabel("Read pairs (count)")
    plt.xlabel("Population")
    plt.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.6)
    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(qc_overview_dir / "demux_distribution.png", dpi=150)
    plt.close()

    # -------------------------
    # 终端输出 unmatched rate（按你的阈值）
    # WARNING/STRONG WARNING：只 print（不停机、不询问）
    # ABNORMAL：才询问 Stop now? (y/N)
    # -------------------------
    print("\n[Stage 2.5] Demux QC overview generated.")
    print(f"[Stage 2.5] Overall unmatched_rate = {unmatched_rate:.3%}")

    if unmatched_rate > 0.50:
        print("ABNORMAL: unmatched_rate > 50%. Recommended: stop and check index/anchor/orientation.")
        ans = input("Stop now? (y/N) ").strip().lower()
        if ans == "y":
            # ✅ 不用 SystemExit（run_pipeline 抓不到），用统一返回值
            return {
                "status": "aborted",
                "qc_details_dir": str(qc_details_dir),
                "qc_overview_dir": str(qc_overview_dir),
                "unmatched_rate": unmatched_rate,
                "reason": "ABNORMAL unmatched_rate > 50% and user chose to stop",
            }
    elif unmatched_rate > 0.35:
        print("STRONG WARNING: unmatched_rate > 35%. Please inspect demux setup.")
    elif unmatched_rate > 0.20:
        print("WARNING: unmatched_rate > 20%. Please inspect index/anchor/orientation.")

    return {
        "status": "success",
        "qc_details_dir": str(qc_details_dir),
        "qc_overview_dir": str(qc_overview_dir),
        "unmatched_rate": unmatched_rate,
    }