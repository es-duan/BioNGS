from __future__ import annotations

import pickle
from pathlib import Path
from statistics import median
from typing import Dict, List, Tuple


def _load_umi_dict(pkl_path: Path) -> Dict[str, int]:
    with open(pkl_path, "rb") as f:
        return pickle.load(f)


def _summarize_family_sizes(umi_counts: Dict[str, int]) -> Dict:
    sizes = list(umi_counts.values())
    if not sizes:
        return {
            "unique_families": 0,
            "mean_reads_per_family": 0.0,
            "median_family_size": 0,
            "max_family_size": 0,
            "n_size_1": 0, "n_size_2": 0, "n_size_3": 0, "n_size_gt3": 0,
            "pct_size_1": 0.0, "pct_size_2": 0.0, "pct_size_3": 0.0, "pct_size_gt3": 0.0,
            "n_ge2": 0, "pct_ge2": 0.0,
            "n_ge3": 0, "pct_ge3": 0.0,
        }

    total_fam = len(sizes)
    mean_sz = sum(sizes) / total_fam
    med_sz = int(median(sizes))
    max_sz = max(sizes)

    n1 = sum(1 for x in sizes if x == 1)
    n2 = sum(1 for x in sizes if x == 2)
    n3 = sum(1 for x in sizes if x == 3)
    ngt3 = sum(1 for x in sizes if x > 3)

    n_ge2 = sum(1 for x in sizes if x >= 2)
    n_ge3 = sum(1 for x in sizes if x >= 3)

    return {
        "unique_families": total_fam,
        "mean_reads_per_family": mean_sz,
        "median_family_size": med_sz,
        "max_family_size": max_sz,
        "n_size_1": n1, "n_size_2": n2, "n_size_3": n3, "n_size_gt3": ngt3,
        "pct_size_1": n1 / total_fam, "pct_size_2": n2 / total_fam, "pct_size_3": n3 / total_fam, "pct_size_gt3": ngt3 / total_fam,
        "n_ge2": n_ge2, "pct_ge2": n_ge2 / total_fam,
        "n_ge3": n_ge3, "pct_ge3": n_ge3 / total_fam,
    }


def stage3_5_umi_quality(
    exp: str,
    run_tag: str,
    results_exp_dir: Path,
    stage3_result: Dict,
) -> Dict:
    """
    读取 Stage3 每群的 *_UMI_dict.pkl，输出：
      results/<exp>/UMI_quality/<run_tag>/<pop>/UMI_quality_summary.txt
    并在终端打印完全相同内容。
    """
    out_root = results_exp_dir / "UMI_quality" / run_tag
    out_root.mkdir(parents=True, exist_ok=True)

    per_pop = stage3_result.get("per_population", {})
    if not per_pop:
        raise ValueError("stage3_result.per_population missing or empty.")

    all_summaries: Dict[str, Dict] = {}

    for pop, info in per_pop.items():
        pkl_path = Path(info["umi_dict_path"])
        umi_counts = _load_umi_dict(pkl_path)

        total_pairs = int(info.get("total_pairs", 0))
        extracted_pairs = int(info.get("extracted_pairs", 0))
        unmatched_pairs = int(info.get("unmatched_pairs", 0))

        fam_summary = _summarize_family_sizes(umi_counts)

        lines: List[str] = []
        lines.append(f"Experiment: {exp}")
        lines.append(f"Run tag: {run_tag}")
        lines.append(f"Population: {pop}")
        lines.append("-" * 60)
        lines.append(f"Total read pairs: {total_pairs}")
        lines.append(f"Pairs with UMI extracted: {extracted_pairs}")
        lines.append(f"Unmatched (failed extraction) pairs: {unmatched_pairs}")
        lines.append(f"Unique UMI families: {fam_summary['unique_families']}")
        lines.append(f"Mean reads per family: {fam_summary['mean_reads_per_family']:.3f}")
        lines.append(f"Median family size: {fam_summary['median_family_size']}")
        lines.append(f"Max family size: {fam_summary['max_family_size']}")
        lines.append("-" * 60)
        lines.append("Family size distribution (count / percent of families):")
        lines.append(f"  size=1 : {fam_summary['n_size_1']}  ({fam_summary['pct_size_1']:.3%})")
        lines.append(f"  size=2 : {fam_summary['n_size_2']}  ({fam_summary['pct_size_2']:.3%})")
        lines.append(f"  size=3 : {fam_summary['n_size_3']}  ({fam_summary['pct_size_3']:.3%})")
        lines.append(f"  size>3 : {fam_summary['n_size_gt3']}  ({fam_summary['pct_size_gt3']:.3%})")
        lines.append("-" * 60)
        lines.append(f"Families size>=2: {fam_summary['n_ge2']} ({fam_summary['pct_ge2']:.3%})")
        lines.append(f"Families size>=3: {fam_summary['n_ge3']} ({fam_summary['pct_ge3']:.3%})")
        lines.append("")

        text = "\n".join(lines)

        pop_out = out_root / pop
        pop_out.mkdir(parents=True, exist_ok=True)
        out_txt = pop_out / "UMI_quality_summary.txt"
        out_txt.write_text(text, encoding="utf-8")

        # 终端打印完全相同内容
        print("\n" + text)

        all_summaries[pop] = {
            "out_txt": str(out_txt),
            "family_summary": fam_summary,
            "counts": {
                "total_pairs": total_pairs,
                "extracted_pairs": extracted_pairs,
                "unmatched_pairs": unmatched_pairs,
            },
        }

    return {
        "status": "success",
        "umi_quality_root": str(out_root),
        "per_population": all_summaries,
    }