from __future__ import annotations

import gzip
import pickle
import subprocess
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

import pandas as pd


# -----------------------------
# IO helpers
# -----------------------------
def _open_in(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def _open_out(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "wt")
    return open(path, "w", encoding="utf-8")


def _iter_fastq(path: Path) -> Iterator[Tuple[str, str, str, str]]:
    with _open_in(path) as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline()
            p = f.readline()
            q = f.readline()
            if (not s) or (not p) or (not q):
                # 读到半截：直接抛错（这是“坏记录”，按你要求应停机）
                raise RuntimeError(f"FASTQ truncated/bad record: {path} header={h.strip()}")
            yield h.strip(), s.strip(), p.strip(), q.strip()


def _write_fastq_record(w, h: str, s: str, p: str, q: str) -> None:
    w.write(h + "\n")
    w.write(s + "\n")
    w.write(p + "\n")
    w.write(q + "\n")


def _run_cmd(cmd: List[str]) -> None:
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed ({p.returncode}): {' '.join(cmd)}\n{p.stdout}")


# -----------------------------
# UMI pattern
# -----------------------------
@dataclass(frozen=True)
class UmiPattern:
    prefix: str     # left primer/anchor part BEFORE Ns
    suffix: str     # right primer/anchor part AFTER Ns
    umi_len: int    # number of N's

    @staticmethod
    def from_primer_with_N(primer: str) -> "UmiPattern":
        """
        primer example: ACGTNNNNNNNNNNGGTACA
        """
        primer = primer.strip().upper()
        if "N" not in primer:
            raise ValueError(f"Primer must contain Ns (UMI region). Got: {primer}")

        first_n = primer.find("N")
        last_n = primer.rfind("N")
        prefix = primer[:first_n]
        suffix = primer[last_n + 1 :]
        umi_len = (last_n - first_n + 1)
        if umi_len <= 0:
            raise ValueError(f"Invalid N block in primer: {primer}")
        return UmiPattern(prefix=prefix, suffix=suffix, umi_len=umi_len)


def _try_extract_umi_and_trim(seq: str, qual: str, pat: UmiPattern) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
    """
    Returns: (umi, trimmed_seq, trimmed_qual, reason)
      - if success: reason=None
      - if failed: umi/trimmed_seq/trimmed_qual=None and reason is a string
    Reasons:
      - qual_len_mismatch
      - too_short
      - missing_prefix_left
      - missing_anchor_right
      - invalid_umi_chars
      - other_parse_fail
    """
    seq = seq.strip().upper()
    qual = qual.strip()

    if len(seq) != len(qual):
        return None, None, None, "qual_len_mismatch"

    pre_start = 0
    pre_end = len(pat.prefix)

    umi_start = pre_end
    umi_end = umi_start + pat.umi_len

    suf_start = umi_end
    suf_end = suf_start + len(pat.suffix)

    if len(seq) < suf_end:
        return None, None, None, "too_short"

    if seq[pre_start:pre_end] != pat.prefix:
        return None, None, None, "missing_prefix_left"

    if seq[suf_start:suf_end] != pat.suffix:
        return None, None, None, "missing_anchor_right"

    umi = seq[umi_start:umi_end]
    # 允许 N 吗？通常不应出现（测序质量差或读到不明碱基）
    if any(ch not in "ACGT" for ch in umi):
        return None, None, None, "invalid_umi_chars"

    # trim: 去掉 prefix + UMI + suffix，保留后面主体序列
    trimmed_seq = seq[suf_end:]
    trimmed_qual = qual[suf_end:]
    return umi, trimmed_seq, trimmed_qual, None


def _find_demux_pairs(demux_dir: Path) -> List[Tuple[str, Path, Path]]:
    pairs: List[Tuple[str, Path, Path]] = []

    for pop_dir in sorted([p for p in demux_dir.iterdir() if p.is_dir()]):
        pop = pop_dir.name
        r1s = sorted(list(pop_dir.glob("*_R1.fastq")) + list(pop_dir.glob("*_R1.fastq.gz")))
        r2s = sorted(list(pop_dir.glob("*_R2.fastq")) + list(pop_dir.glob("*_R2.fastq.gz")))
        if r1s and r2s:
            pairs.append((pop, r1s[0], r2s[0]))

    um1 = list(demux_dir.glob("unmatched_R1.fastq")) + list(demux_dir.glob("unmatched_R1.fastq.gz"))
    um2 = list(demux_dir.glob("unmatched_R2.fastq")) + list(demux_dir.glob("unmatched_R2.fastq.gz"))
    if um1 and um2:
        pairs.append(("unmatched", um1[0], um2[0]))

    if not pairs:
        raise FileNotFoundError(f"No demux FASTQ pairs found under: {demux_dir}")
    return pairs


def _write_reason_report(
    out_txt: Path,
    exp: str,
    run_tag: str,
    pop: str,
    total_pairs: int,
    extracted_pairs: int,
    unmatched_pairs: int,
    r1_reasons: Counter,
    r2_reasons: Counter,
    pair_reasons: Counter,
) -> None:
    def _dump_counter(title: str, c: Counter) -> List[str]:
        lines = [title]
        if sum(c.values()) == 0:
            lines.append("  (none)")
            return lines
        denom = sum(c.values())
        for k, v in c.most_common():
            lines.append(f"  {k:22s} {v:10d}  ({(v/denom):6.2%})")
        return lines

    out_txt.parent.mkdir(parents=True, exist_ok=True)
    with out_txt.open("w", encoding="utf-8") as f:
        f.write(f"Experiment: {exp}\n")
        f.write(f"Run tag:    {run_tag}\n")
        f.write(f"Population: {pop}\n")
        f.write("=" * 72 + "\n")
        f.write(f"Total read pairs:        {total_pairs}\n")
        f.write(f"Pairs with UMI extracted:{extracted_pairs}\n")
        f.write(f"Unmatched/failed pairs:  {unmatched_pairs}\n")
        f.write(f"Unmatched rate:          {(unmatched_pairs/total_pairs if total_pairs else 0.0):.3%}\n")
        f.write("=" * 72 + "\n\n")

        for line in _dump_counter("[Pair-level fail reasons] (first failure reason)", pair_reasons):
            f.write(line + "\n")
        f.write("\n")
        for line in _dump_counter("[R1 fail reasons]", r1_reasons):
            f.write(line + "\n")
        f.write("\n")
        for line in _dump_counter("[R2 fail reasons]", r2_reasons):
            f.write(line + "\n")


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
        results/<exp>/umi_extracted/<run_tag>/P*/umi_extraction_summary.txt   <-- NEW
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
    stage_pair_reason = Counter()
    stage_r1_reason = Counter()
    stage_r2_reason = Counter()

    per_pop_stats: Dict[str, Dict] = {}

    for pop, r1_path, r2_path in pop_pairs:
        pop_dir = out_root / pop
        pop_dir.mkdir(parents=True, exist_ok=True)

        out_r1 = pop_dir / f"{pop}_R1.fastq.gz"
        out_r2 = pop_dir / f"{pop}_R2.fastq.gz"
        out_um1 = pop_dir / f"{pop}_unmatched_UMI_R1.fastq.gz"
        out_um2 = pop_dir / f"{pop}_unmatched_UMI_R2.fastq.gz"
        out_pkl = pop_dir / f"{pop}_UMI_dict.pkl"
        out_reason_txt = pop_dir / "umi_extraction_summary.txt"

        umi_counts: Dict[str, int] = {}
        total_pairs = 0
        ok_pairs = 0
        bad_pairs = 0

        r1_reason = Counter()
        r2_reason = Counter()
        pair_reason = Counter()

        with _open_out(out_r1) as w1, _open_out(out_r2) as w2, _open_out(out_um1) as wu1, _open_out(out_um2) as wu2:
            it1 = _iter_fastq(r1_path)
            it2 = _iter_fastq(r2_path)

            for rec1, rec2 in zip(it1, it2):
                (h1, s1, p1, q1) = rec1
                (h2, s2, p2, q2) = rec2
                total_pairs += 1

                umi1, ts1, tq1, reason1 = _try_extract_umi_and_trim(s1, q1, pat_r1)
                umi2, ts2, tq2, reason2 = _try_extract_umi_and_trim(s2, q2, pat_r2)

                if reason1 is not None:
                    r1_reason[reason1] += 1
                if reason2 is not None:
                    r2_reason[reason2] += 1

                if (reason1 is not None) or (reason2 is not None):
                    bad_pairs += 1
                    # pair-level: pick “first failure reason” for readability
                    first = reason1 if reason1 is not None else reason2
                    pair_reason[first] += 1
                    _write_fastq_record(wu1, h1, s1, p1, q1)
                    _write_fastq_record(wu2, h2, s2, p2, q2)
                    continue

                umi = f"{umi1}-{umi2}"
                umi_counts[umi] = umi_counts.get(umi, 0) + 1

                new_h1 = f"{h1} UMI={umi}"
                new_h2 = f"{h2} UMI={umi}"

                ok_pairs += 1
                _write_fastq_record(w1, new_h1, ts1, p1, tq1)
                _write_fastq_record(w2, new_h2, ts2, p2, tq2)

        with open(out_pkl, "wb") as pf:
            pickle.dump(umi_counts, pf)

        _write_reason_report(
            out_reason_txt,
            exp=exp,
            run_tag=run_tag,
            pop=pop,
            total_pairs=total_pairs,
            extracted_pairs=ok_pairs,
            unmatched_pairs=bad_pairs,
            r1_reasons=r1_reason,
            r2_reasons=r2_reason,
            pair_reasons=pair_reason,
        )

        unmatched_rate = (bad_pairs / total_pairs) if total_pairs else 0.0
        per_pop_stats[pop] = {
            "total_pairs": total_pairs,
            "extracted_pairs": ok_pairs,
            "unmatched_pairs": bad_pairs,
            "unmatched_rate": unmatched_rate,
            "umi_dict_path": str(out_pkl),
            "out_dir": str(pop_dir),
            "reason_report": str(out_reason_txt),
            "fail_reason_pair": dict(pair_reason),
            "fail_reason_r1": dict(r1_reason),
            "fail_reason_r2": dict(r2_reason),
        }

        stage_total_pairs += total_pairs
        stage_unmatched_pairs += bad_pairs
        stage_pair_reason.update(pair_reason)
        stage_r1_reason.update(r1_reason)
        stage_r2_reason.update(r2_reason)

    stage_unmatched_rate = (stage_unmatched_pairs / stage_total_pairs) if stage_total_pairs else 0.0

    print("\n[Stage 3] UMI extraction completed.")
    print(f"[Stage 3] Overall unmatched_rate = {stage_unmatched_rate:.3%}")

    # 终端给“透明化”摘要（Top 5）
    if sum(stage_pair_reason.values()) > 0:
        print("[Stage 3] Top pair-level fail reasons:")
        denom = sum(stage_pair_reason.values())
        for k, v in stage_pair_reason.most_common(5):
            print(f"  - {k:22s} {v:10d}  ({(v/denom):.2%})")

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
        "overall_fail_reason_pair": dict(stage_pair_reason),
        "overall_fail_reason_r1": dict(stage_r1_reason),
        "overall_fail_reason_r2": dict(stage_r2_reason),
        "umi_patterns": {
            "r1_prefix": pat_r1.prefix,
            "r1_suffix": pat_r1.suffix,
            "r1_umi_len": pat_r1.umi_len,
            "r2_prefix": pat_r2.prefix,
            "r2_suffix": pat_r2.suffix,
            "r2_umi_len": pat_r2.umi_len,
        },
    }