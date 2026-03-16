from __future__ import annotations

import gzip
import pickle
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple
import subprocess
import shutil
import sys

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

        import re

        m = re.search(r"N{10}", primer)
        if m is None:
            raise ValueError(f"Primer must contain exactly 10 consecutive Ns (UMI). Got: {primer}")

        start = m.start()
        end = m.end()

        prefix = primer[:start]
        suffix = primer[end:]
        umi_len = 10

        # 检查是否有额外 N（不允许）
        outside = primer[:start] + primer[end:]
        if "N" in outside:
            raise ValueError(f"Primer contains extra Ns outside the 10N UMI block: {primer}")

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

    # pre_start = 0
    # pre_end = len(pat.prefix)
    #
    # umi_start = pre_end
    # umi_end = umi_start + pat.umi_len
    #
    # suf_start = umi_end
    # suf_end = suf_start + len(pat.suffix)
    #
    # if len(seq) < suf_end:
    #     return None, None, None, "too_short"
    #
    # if seq[pre_start:pre_end] != pat.prefix:
    #     return None, None, None, "missing_prefix_left"
    #
    # if seq[suf_start:suf_end] != pat.suffix:
    #     return None, None, None, "missing_anchor_right"

    pos = seq.find(pat.prefix)
    if pos < 0:
        return None, None, None, "missing_prefix_anywhere"

    umi_start = pos + len(pat.prefix)
    umi_end = umi_start + pat.umi_len

    suf_start = umi_end
    suf_end = suf_start + len(pat.suffix)

    if len(seq) < suf_end:
        return None, None, None, "too_short"

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

def _try_extract_umi_and_trim_relaxed(seq: str, qual: str, pat: UmiPattern) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
    """
    RELAXED extraction:
      - allow pat.prefix to appear anywhere in the read (seq.find)
      - still requires exact match for pat.suffix at the expected position (after UMI)
      - UMI length is fixed by pat.umi_len

    Equivalent to: find(prefix) -> take next umi_len as UMI -> check suffix -> trim after suffix.
    """
    seq = seq.strip().upper()
    qual = qual.strip()

    if len(seq) != len(qual):
        return None, None, None, "qual_len_mismatch"

    pos = seq.find(pat.prefix)
    if pos < 0:
        return None, None, None, "missing_prefix_anywhere"

    umi_start = pos + len(pat.prefix)
    umi_end = umi_start + pat.umi_len
    suf_start = umi_end
    suf_end = suf_start + len(pat.suffix)

    if len(seq) < suf_end:
        return None, None, None, "too_short"

    if seq[suf_start:suf_end] != pat.suffix:
        # DEBUG: inspect failed anchor cases
        #print("\n[DEBUG missing_anchor_right]")
        #print("seq head:", seq[:80])
        #print("expected suffix:", pat.suffix)
        #print("observed around expected:", seq[umi_end - 5:umi_end + 20])

        return None, None, None, "missing_anchor_right"

    umi = seq[umi_start:umi_end]
    if any(ch not in "ACGT" for ch in umi):
        return None, None, None, "invalid_umi_chars"

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
        results/<exp>/umi_extracted/<run_tag>/P*/umi_extraction_summary.txt
      统计 unmatched_rate（pair 口径）并按阈值提示/询问停机（仅 ABNORMAL 才询问）

    New behavior:
      - 先跑严格模式
      - 严格模式 unmatched_rate 过高则询问是否用“宽松模式”重跑 Stage 3
      - 选择重跑：删除 Stage 3 输出目录并重建，覆盖旧数据
      - 最终给出重跑后的报告（如果没重跑则给严格模式报告）
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

    def _run_once(extractor_fn, mode_name: str) -> Dict:
        stage_total_pairs = 0
        stage_unmatched_pairs = 0
        stage_pair_reason = Counter()
        stage_r1_reason = Counter()
        stage_r2_reason = Counter()
        SHORT_MINLEN = 150  # reporting-only: used to compute "effective unmatched rate"

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
            # reporting-only counters (do NOT change writing behavior)
            short_pairs = 0  # pairs where either mate < SHORT_MINLEN
            usable_pairs = 0  # pairs where both mates >= SHORT_MINLEN
            usable_bad_pairs = 0  # failed extraction among usable_pairs
            usable_pair_reason = Counter()

            r1_reason = Counter()
            r2_reason = Counter()
            pair_reason = Counter()

            with _open_out(out_r1) as w1, _open_out(out_r2) as w2, _open_out(out_um1) as wu1, _open_out(out_um2) as wu2:
                it1 = _iter_fastq(r1_path)
                it2 = _iter_fastq(r2_path)

                for rec1, rec2 in zip(it1, it2):

                    def _rid(h: str) -> str:
                        x = h.split()[0]
                        if x.endswith("/1") or x.endswith("/2"):
                            x = x[:-2]
                        return x

                    (h1, s1, p1, q1) = rec1
                    (h2, s2, p2, q2) = rec2
                    total_pairs += 1

                    # sanity: R1/R2 pairing
                    if _rid(h1) != _rid(h2):
                        raise RuntimeError(f"R1/R2 read-id mismatch: {h1} vs {h2}")

                    # reporting-only: short vs usable (does NOT skip processing)
                    is_short = (len(s1) < SHORT_MINLEN) or (len(s2) < SHORT_MINLEN)
                    if is_short:
                        short_pairs += 1
                    else:
                        usable_pairs += 1

                    umi1, ts1, tq1, reason1 = extractor_fn(s1, q1, pat_r1)
                    umi2, ts2, tq2, reason2 = extractor_fn(s2, q2, pat_r2)

                    if reason1 is not None:
                        r1_reason[reason1] += 1
                    if reason2 is not None:
                        r2_reason[reason2] += 1

                    if (reason1 is not None) or (reason2 is not None):
                        bad_pairs += 1
                        first = reason1 if reason1 is not None else reason2
                        pair_reason[first] += 1
                        # reporting-only: failures among usable pairs
                        if not is_short:
                            usable_bad_pairs += 1
                            usable_pair_reason[first] += 1
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

                # 结束后检查是否有一边多记录
                try:
                    next(it1)
                    raise RuntimeError(f"R1 has extra records beyond R2: {r1_path}")
                except StopIteration:
                    pass

                try:
                    next(it2)
                    raise RuntimeError(f"R2 has extra records beyond R1: {r2_path}")
                except StopIteration:
                    pass

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
                # reporting-only
                "short_pairs_lt150": short_pairs,
                "short_rate_lt150": (short_pairs / total_pairs) if total_pairs else 0.0,
                "usable_pairs_ge150": usable_pairs,
                "usable_unmatched_pairs": usable_bad_pairs,
                "usable_unmatched_rate": (usable_bad_pairs / usable_pairs) if usable_pairs else 0.0,
                "usable_fail_reason_pair": dict(usable_pair_reason),
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

        # -----------------------------
        # Friendly terminal reporting
        # -----------------------------
        def _pfx(mode: str) -> str:
            # Hide ":strict" to reduce noise; keep suffix for other modes (e.g., relaxed)
            return "[Stage 3]" if mode == "strict" else f"[Stage 3:{mode}]"

        print(f"\n{_pfx(mode_name)} UMI extraction completed.")
        print(f"{_pfx(mode_name)} Raw unmatched_rate (all pairs) = {stage_unmatched_rate:.3%}")

        # Effective rate (>=SHORT_MINLEN only) aggregated from per_pop_stats
        stage_short = 0
        stage_usable = 0
        stage_usable_bad = 0
        stage_usable_reason = Counter()

        for d in (per_pop_stats or {}).values():
            if not isinstance(d, dict):
                continue
            stage_short += int(d.get("short_pairs_lt150", 0))
            stage_usable += int(d.get("usable_pairs_ge150", 0))
            stage_usable_bad += int(d.get("usable_unmatched_pairs", 0))
            stage_usable_reason.update(d.get("usable_fail_reason_pair") or {})

        short_rate = (stage_short / stage_total_pairs) if stage_total_pairs else 0.0
        usable_rate = (stage_usable_bad / stage_usable) if stage_usable else 0.0

        print(f"{_pfx(mode_name)} Short pairs (<{SHORT_MINLEN} on either mate) = {stage_short} ({short_rate:.3%})")
        print(f"{_pfx(mode_name)} Effective unmatched_rate (>= {SHORT_MINLEN} only) = {usable_rate:.3%}")

        # Keep raw fail reasons for diagnosis
        if sum(stage_pair_reason.values()) > 0:
            print(f"{_pfx(mode_name)} Top pair-level fail reasons (all pairs):")
            denom = sum(stage_pair_reason.values())
            for k, v in stage_pair_reason.most_common(5):
                print(f"  - {k:22s} {v:10d}  ({(v / denom):.2%})")

        # Usable-only fail reasons (what users should focus on)
        if sum(stage_usable_reason.values()) > 0:
            print(f"{_pfx(mode_name)} Top pair-level fail reasons (>= {SHORT_MINLEN} only):")
            denom = sum(stage_usable_reason.values())
            for k, v in stage_usable_reason.most_common(5):
                print(f"  - {k:22s} {v:10d}  ({(v / denom):.2%})")

        # ---- RETURN RESULT (CRITICAL) ----
        return {
            "status": "success",
            "mode": mode_name,
            "umi_extracted_root": str(out_root),
            "overall_unmatched_rate": stage_unmatched_rate,
            "per_population": per_pop_stats,
            "overall_fail_reason_pair": dict(stage_pair_reason),
            "overall_fail_reason_r1": dict(stage_r1_reason),
            "overall_fail_reason_r2": dict(stage_r2_reason),
        }

    # ---------- run strict once ----------
    strict_result = _run_once(_try_extract_umi_and_trim, "strict")
    # 先保证最小可运行：直接返回 strict_result
    return strict_result