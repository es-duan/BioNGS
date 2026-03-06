from __future__ import annotations

import gzip
from pathlib import Path
from typing import Tuple, Dict, Any, Optional
import subprocess


def open_text_maybe_gz(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8", errors="replace")


def count_fastq_reads(path: Path) -> int:
    """
    通过 wc -l 计数（gz 用 zcat），read数 = 行数/4
    注意：这是一次性统计，可能较慢，但只在抽样比例计算时用。
    """
    if str(path).endswith(".gz"):
        cmd = f"zcat {str(path)!s} | wc -l"
    else:
        cmd = f"wc -l < {str(path)!s}"
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Failed counting FASTQ lines: {path}\n{p.stdout}")
    lines = int(p.stdout.strip())
    return lines // 4


def write_fastq_sample_pairs(
    r1: Path,
    r2: Path,
    out_r1: Path,
    out_r2: Path,
    n_pairs: int,
) -> Dict[str, Any]:
    """
    顺序抽样前 n_pairs 对（不随机，工程上更可复现）
    """
    out_r1.parent.mkdir(parents=True, exist_ok=True)
    out_r2.parent.mkdir(parents=True, exist_ok=True)

    wrote = 0
    with open_text_maybe_gz(r1) as f1, open_text_maybe_gz(r2) as f2, open(out_r1, "w") as o1, open(out_r2, "w") as o2:
        while wrote < n_pairs:
            h1 = f1.readline()
            h2 = f2.readline()
            if not h1 or not h2:
                break
            s1 = f1.readline(); s2 = f2.readline()
            p1 = f1.readline(); p2 = f2.readline()
            q1 = f1.readline(); q2 = f2.readline()
            if not q1 or not q2:
                break

            o1.write(h1); o1.write(s1); o1.write(p1); o1.write(q1)
            o2.write(h2); o2.write(s2); o2.write(p2); o2.write(q2)
            wrote += 1

    return {"sample_pairs": wrote, "requested_pairs": n_pairs, "out_r1": str(out_r1), "out_r2": str(out_r2)}


def preview_pair_short_rate_from_trimmed(
    r1: Path,
    r2: Path,
    min_len: int,
) -> Dict[str, Any]:
    """
    pair-level short_rate: 若任一 mate 长度 < min_len，该pair为 short
    输出：total_pairs, short_pairs, kept_pairs, short_rate, kept_rate
    """
    total = 0
    short = 0

    with open_text_maybe_gz(r1) as f1, open_text_maybe_gz(r2) as f2:
        while True:
            h1 = f1.readline()
            h2 = f2.readline()
            if not h1 or not h2:
                break
            s1 = f1.readline().strip()
            s2 = f2.readline().strip()
            f1.readline(); f2.readline()
            q1 = f1.readline()
            q2 = f2.readline()
            if not q1 or not q2:
                break

            total += 1
            if len(s1) < min_len or len(s2) < min_len:
                short += 1

    kept = total - short
    return {
        "total_pairs": total,
        "short_pairs": short,
        "kept_pairs": kept,
        "short_rate": (short / total) if total else 0.0,
        "kept_rate": (kept / total) if total else 0.0,
        "estimated_dropped_reads": short,  # pair口径：预计丢弃 pairs
    }


def filter_pairs_by_min_len(
    r1_in: Path,
    r2_in: Path,
    r1_out_gz: Path,
    r2_out_gz: Path,
    min_len: int,
) -> Dict[str, Any]:
    """
    真正执行 length filtering：trim 后任一 mate < min_len 的 pair 直接丢
    输出为 .fastq.gz
    """
    r1_out_gz.parent.mkdir(parents=True, exist_ok=True)
    r2_out_gz.parent.mkdir(parents=True, exist_ok=True)

    total = 0
    kept = 0
    dropped = 0

    with open_text_maybe_gz(r1_in) as f1, open_text_maybe_gz(r2_in) as f2, gzip.open(r1_out_gz, "wt") as o1, gzip.open(r2_out_gz, "wt") as o2:
        while True:
            h1 = f1.readline()
            h2 = f2.readline()
            if not h1 or not h2:
                break
            s1 = f1.readline()
            s2 = f2.readline()
            p1 = f1.readline()
            p2 = f2.readline()
            q1 = f1.readline()
            q2 = f2.readline()
            if not q1 or not q2:
                break

            total += 1
            L1 = len(s1.strip())
            L2 = len(s2.strip())

            if L1 < min_len or L2 < min_len:
                dropped += 1
                continue

            kept += 1
            o1.write(h1); o1.write(s1); o1.write(p1); o1.write(q1)
            o2.write(h2); o2.write(s2); o2.write(p2); o2.write(q2)

    return {
        "total_pairs": total,
        "kept_pairs": kept,
        "dropped_pairs": dropped,
        "kept_rate": (kept / total) if total else 0.0,
        "short_rate": (dropped / total) if total else 0.0,
    }