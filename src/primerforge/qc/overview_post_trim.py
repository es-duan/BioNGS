from __future__ import annotations

import gzip
from pathlib import Path
from collections import Counter
import matplotlib.pyplot as plt


def open_fastq(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def compute_post_trim_overview(fastq: Path) -> dict:
    """
    Post-trim overview (no short-rate judgement here).
    Computes:
      - total_reads
      - total_bases
      - q30_rate (base-level)
      - mean_len
      - median_len
      - length_counter
    """
    length_counter = Counter()
    total_reads = 0
    total_bases = 0
    q30_bases = 0

    with open_fastq(fastq) as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # plus
            qual = f.readline().strip()

            L = len(seq)
            length_counter[L] += 1
            total_reads += 1
            total_bases += L

            for ch in qual:
                if ord(ch) - 33 >= 30:
                    q30_bases += 1

    q30_rate = (q30_bases / total_bases) if total_bases else 0.0
    mean_len = (total_bases / total_reads) if total_reads else 0.0

    # median length from counter
    median_len = 0
    if total_reads:
        target = (total_reads + 1) // 2
        running = 0
        for L in sorted(length_counter.keys()):
            running += length_counter[L]
            if running >= target:
                median_len = L
                break

    return {
        "total_reads": total_reads,
        "total_bases": total_bases,
        "q30_bases": q30_bases,  # <- 为 overall 聚合保留
        "q30_rate": q30_rate,
        "mean_len": mean_len,
        "median_len": median_len,
        "length_counter": length_counter,
    }


def aggregate_post_trim_metrics(metrics_list: list[dict]) -> dict:
    """
    Aggregate multiple compute_post_trim_overview() outputs into one.
    Used for overall summary across populations.
    """
    total_reads = sum(m.get("total_reads", 0) for m in metrics_list)
    total_bases = sum(m.get("total_bases", 0) for m in metrics_list)
    q30_bases = sum(m.get("q30_bases", 0) for m in metrics_list)

    length_counter = Counter()
    for m in metrics_list:
        length_counter.update(m.get("length_counter", {}))

    q30_rate = (q30_bases / total_bases) if total_bases else 0.0
    mean_len = (total_bases / total_reads) if total_reads else 0.0

    median_len = 0
    if total_reads:
        target = (total_reads + 1) // 2
        running = 0
        for L in sorted(length_counter.keys()):
            running += length_counter[L]
            if running >= target:
                median_len = L
                break

    return {
        "total_reads": total_reads,
        "total_bases": total_bases,
        "q30_bases": q30_bases,
        "q30_rate": q30_rate,
        "mean_len": mean_len,
        "median_len": median_len,
        "length_counter": length_counter,
    }


def write_post_trim_overview(metrics: dict, outdir: Path, title: str) -> None:
    """
    Writes:
      - post_trim_overview_summary.txt
      - read_length_frequency.png  (line/step plot, NOT bar chart)
    """
    outdir.mkdir(parents=True, exist_ok=True)

    summary_txt = outdir / "post_trim_overview_summary.txt"
    with summary_txt.open("w", encoding="utf-8") as f:
        f.write(f"{title}\n")
        f.write("=" * 60 + "\n")
        f.write(f"Total reads: {metrics['total_reads']}\n")
        f.write(f"Total bases: {metrics['total_bases']}\n")
        f.write(f"Q30 rate (base-level): {metrics['q30_rate']:.3%}\n")
        f.write(f"Mean read length: {metrics['mean_len']:.2f}\n")
        f.write(f"Median read length: {metrics['median_len']}\n")

    lengths = sorted(metrics["length_counter"].keys())
    counts = [metrics["length_counter"][L] for L in lengths]

    plt.figure(figsize=(7, 4))
    plt.plot(lengths, counts, drawstyle="steps-mid", linewidth=1.0)
    plt.title(title, fontsize=12, pad=10)
    plt.xlabel("Read length (bp)")
    plt.ylabel("Read count")
    plt.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.6)

    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(outdir / "read_length_frequency.png", dpi=150)
    plt.close()