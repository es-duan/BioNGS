from __future__ import annotations
import gzip
from pathlib import Path
from collections import Counter
import matplotlib.pyplot as plt


def open_fastq(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def compute_raw_overview(fastq: Path, short_len: int = 150):
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
            _plus = f.readline()
            qual = f.readline().strip()

            L = len(seq)
            length_counter[L] += 1
            total_reads += 1
            total_bases += L

            for ch in qual:
                if ord(ch) - 33 >= 30:
                    q30_bases += 1

    short_reads = sum(c for l, c in length_counter.items() if l < short_len)
    short_rate = short_reads / total_reads if total_reads else 0
    q30_rate = q30_bases / total_bases if total_bases else 0

    return {
        "total_reads": total_reads,
        "short_rate": short_rate,
        "q30_rate": q30_rate,
        "length_counter": length_counter
    }


def write_overview_report(metrics: dict, outdir: Path, title: str = "Read Length Distribution"):
    outdir.mkdir(parents=True, exist_ok=True)

    summary_txt = outdir / "raw_overview_summary.txt"
    with open(summary_txt, "w", encoding="utf-8") as f:
        f.write(f"Total reads: {metrics['total_reads']}\n")
        f.write(f"Short rate (<150bp): {metrics['short_rate']:.3%}\n")
        f.write(f"Q30 rate (base-level): {metrics['q30_rate']:.3%}\n")

    lengths = list(metrics["length_counter"].keys())
    counts = list(metrics["length_counter"].values())

    plt.figure(figsize=(7, 4))
    plt.bar(lengths, counts, width=2, edgecolor="black", linewidth=0.3)
    plt.title(title, fontsize=12, pad=10)
    plt.xlabel("Read length (bp)")
    plt.ylabel("Read count")

    # 横线更好读
    plt.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.6)

    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(outdir / "length_distribution.png", dpi=150)
    plt.close()