from __future__ import annotations

import gzip
import pickle
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

import pandas as pd


def _open_text(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def _open_out(path: Path):
    # 输出统一 gzip，便于后续跑
    if str(path).endswith(".gz"):
        return gzip.open(path, "wt")
    return gzip.open(str(path) + ".gz", "wt")


def _iter_fastq(path: Path) -> Iterator[Tuple[str, str, str, str]]:
    with _open_text(path, "rt") as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline().rstrip("\n")
            p = f.readline()
            q = f.readline().rstrip("\n")
            if not p or not q:
                raise ValueError(f"Malformed FASTQ (incomplete record): {path}")
            yield h.rstrip("\n"), s, p.rstrip("\n"), q


def _write_fastq_record(out, h: str, s: str, p: str, q: str) -> None:
    out.write(h + "\n")
    out.write(s + "\n")
    out.write(p + "\n")
    out.write(q + "\n")


@dataclass
class UmiPattern:
    prefix: str
    suffix: str
    umi_len: int

    @staticmethod
    def from_primer_with_N(primer: str) -> "UmiPattern":
        m = re.search(r"(N+)", primer)
        if not m:
            raise ValueError("Primer does not contain N-run (UMI).")
        umi_len = len(m.group(1))
        prefix = primer[: m.start(1)]
        suffix = primer[m.end(1) :]
        if not prefix or not suffix:
            raise ValueError("Primer must have both prefix and suffix around N-run.")
        return UmiPattern(prefix=prefix, suffix=suffix, umi_len=umi_len)


def _extract_umi_and_trim(seq: str, qual: str, pat: UmiPattern, search_window: int = 60) -> Optional[Tuple[str, str, str]]:
    """
    在 read 前段 search_window 内寻找 prefix + UMI + suffix（严格匹配 prefix/suffix）
    成功：返回 (umi, trimmed_seq, trimmed_qual)
    失败：None
    """
    window = seq[:search_window]
    start = window.find(pat.prefix)
    if start < 0:
        return None

    umi_start = start + len(pat.prefix)
    umi_end = umi_start + pat.umi_len
    suf_start = umi_end
    suf_end = suf_start + len(pat.suffix)

    if len(seq) < suf_end:
        return None

    if seq[suf_start:suf_end] != pat.suffix:
        return None

    umi = seq[umi_start:umi_end]
    # 去掉 primer(含 UMI) 部分：保留 suffix 后的序列
    trimmed_seq = seq[suf_end:]
    trimmed_qual = qual[suf_end:]
    return umi, trimmed_seq, trimmed_qual


def _find_demux_pairs(demux_dir: Path) -> List[Tuple[str, Path, Path]]:
    """
    demux_dir = results/<exp>/demultiplexing/<run_tag>/
      P*/xxx_R1.fastq(.gz), xxx_R2.fastq(.gz)
      unmatched_R1.fastq(.gz), unmatched_R2.fastq(.gz)
    """
    pairs: List[Tuple[str, Path, Path]] = []

    for pop_dir in sorted([p for p in demux_dir.iterdir() if p.is_dir()]):
        pop = pop_dir.name
        r1s = sorted(list(pop_dir.glob("*_R1.fastq")) + list(pop_dir.glob("*_R1.fastq.gz")))
        r2s = sorted(list(pop_dir.glob("*_R2.fastq")) + list(pop_dir.glob("*_R2.fastq.gz")))
        if r1s and r2s:
            pairs.append((pop, r1s[0], r2s[0]))

    # demux unmatched
    um1 = list(demux_dir.glob("unmatched_R1.fastq")) + list(demux_dir.glob("unmatched_R1.fastq.gz"))
    um2 = list(demux_dir.glob("unmatched_R2.fastq")) + list(demux_dir.glob("unmatched_R2.fastq.gz"))
    if um1 and um2:
        pairs.append(("unmatched", um1[0], um2[0]))

    if not pairs:
        raise FileNotFoundError(f"No demux FASTQ pairs found under: {demux_dir}")

    return pairs


def stage3_umi_extract(
    exp: str,
    run_tag: str,
    results_exp_dir: Path,
    umi_primers_csv: Path,
) -> Dict:
    """
    Stage 3:
      输入：results/<exp>/demultiplexing/<run_tag>/ 下的每个 P* FASTQ
      输出：
        results/<exp>/umi_extracted/<run_tag>/P*/P*_R1.fastq.gz
        results/<exp>/umi_extracted/<run_tag>/P*/P*_R2.fastq.gz
        results/<exp>/umi_extracted/<run_tag>/P*/P*_UMI_dict.pkl
        results/<exp>/umi_extracted/<run_tag>/P*/P*_unmatched_UMI_R1.fastq.gz
        results/<exp>/umi_extracted/<run_tag>/P*/P*_unmatched_UMI_R2.fastq.gz
      统计 unmatched_rate（pair 口径）并按阈值提示/询问停机（仅 ABNORMAL 才询问）
    """
    demux_dir = results_exp_dir / "demultiplexing" / run_tag
    out_root = results_exp_dir / "umi_extracted" / run_tag
    out_root.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(umi_primers_csv)
    if not {"f", "r"}.issubset(set(df.columns)):
        raise ValueError(f"UMI primers CSV must contain columns: f,r. Found: {list(df.columns)}")
    f_primer = str(df.loc[0, "f"]).strip()
    r_primer = str(df.loc[0, "r"]).strip()

    pat_r1 = UmiPattern.from_primer_with_N(f_primer)
    pat_r2 = UmiPattern.from_primer_with_N(r_primer)

    pop_pairs = _find_demux_pairs(demux_dir)

    stage_total_pairs = 0
    stage_unmatched_pairs = 0
    per_pop_stats: Dict[str, Dict] = {}

    for pop, r1_path, r2_path in pop_pairs:
        pop_dir = out_root / pop
        pop_dir.mkdir(parents=True, exist_ok=True)

        out_r1 = pop_dir / f"{pop}_R1.fastq.gz"
        out_r2 = pop_dir / f"{pop}_R2.fastq.gz"
        out_um1 = pop_dir / f"{pop}_unmatched_UMI_R1.fastq.gz"
        out_um2 = pop_dir / f"{pop}_unmatched_UMI_R2.fastq.gz"
        out_pkl = pop_dir / f"{pop}_UMI_dict.pkl"

        umi_counts: Dict[str, int] = {}
        total_pairs = 0
        ok_pairs = 0
        bad_pairs = 0

        with _open_out(out_r1) as w1, _open_out(out_r2) as w2, _open_out(out_um1) as wu1, _open_out(out_um2) as wu2:
            it1 = _iter_fastq(r1_path)
            it2 = _iter_fastq(r2_path)

            for rec1, rec2 in zip(it1, it2):
                (h1, s1, p1, q1) = rec1
                (h2, s2, p2, q2) = rec2

                total_pairs += 1

                ex1 = _extract_umi_and_trim(s1, q1, pat_r1)
                ex2 = _extract_umi_and_trim(s2, q2, pat_r2)

                if (ex1 is None) or (ex2 is None):
                    bad_pairs += 1
                    _write_fastq_record(wu1, h1, s1, p1, q1)
                    _write_fastq_record(wu2, h2, s2, p2, q2)
                    continue

                umi1, ts1, tq1 = ex1
                umi2, ts2, tq2 = ex2
                umi = f"{umi1}-{umi2}"

                # 记录 UMI family
                umi_counts[umi] = umi_counts.get(umi, 0) + 1

                # 写入 header tag：在 header 后加一个空格字段
                # e.g. "@readid ... UMI=AAAAA-BBBBB"
                new_h1 = f"{h1} UMI={umi}"
                new_h2 = f"{h2} UMI={umi}"

                ok_pairs += 1
                _write_fastq_record(w1, new_h1, ts1, p1, tq1)
                _write_fastq_record(w2, new_h2, ts2, p2, tq2)

        # 保存 pkl
        with open(out_pkl, "wb") as pf:
            pickle.dump(umi_counts, pf)

        unmatched_rate = (bad_pairs / total_pairs) if total_pairs else 0.0
        per_pop_stats[pop] = {
            "total_pairs": total_pairs,
            "extracted_pairs": ok_pairs,
            "unmatched_pairs": bad_pairs,
            "unmatched_rate": unmatched_rate,
            "umi_dict_path": str(out_pkl),
            "out_dir": str(pop_dir),
        }

        stage_total_pairs += total_pairs
        stage_unmatched_pairs += bad_pairs

    stage_unmatched_rate = (stage_unmatched_pairs / stage_total_pairs) if stage_total_pairs else 0.0

    print("\n[Stage 3] UMI extraction completed.")
    print(f"[Stage 3] Overall unmatched_rate = {stage_unmatched_rate:.3%}")

    if stage_unmatched_rate > 0.50:
        print("ABNORMAL: unmatched_rate > 50%. Recommended: stop and check anchor/primer/orientation/mismatch tolerance.")
        ans = input("Stop now? (y/N) ").strip().lower()
        if ans == "y":
            return {"status": "aborted", "reason": "ABNORMAL unmatched_rate > 50% and user chose to stop"}
    elif stage_unmatched_rate > 0.35:
        print("STRONG WARNING: unmatched_rate > 35%. Please inspect anchor/primer/orientation/mismatch tolerance.")
    elif stage_unmatched_rate > 0.20:
        print("WARNING: unmatched_rate > 20%. Please inspect anchor/primer/orientation/mismatch tolerance.")

    return {
        "status": "success",
        "umi_extracted_root": str(out_root),
        "overall_unmatched_rate": stage_unmatched_rate,
        "per_population": per_pop_stats,
        "umi_patterns": {
            "r1_prefix": pat_r1.prefix,
            "r1_suffix": pat_r1.suffix,
            "r1_umi_len": pat_r1.umi_len,
            "r2_prefix": pat_r2.prefix,
            "r2_suffix": pat_r2.suffix,
            "r2_umi_len": pat_r2.umi_len,
        },
    }