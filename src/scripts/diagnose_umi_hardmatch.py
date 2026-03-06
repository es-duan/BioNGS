#!/usr/bin/env python3
# diagnose_umi_hardmatch.py
# Purpose: Diagnose why strict (prefix at pos0 + exact suffix) UMI extraction fails.
# Works on paired FASTQ (plain or .gz). No dependency besides Python stdlib.

from __future__ import annotations
import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Iterator, Tuple, Optional, Dict, List


def open_maybe_gz(path: str, mode: str = "rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def fastq_iter(path: str) -> Iterator[Tuple[str, str, str]]:
    """
    Yield (name, seq, qual) from FASTQ.
    """
    with open_maybe_gz(path, "rt") as f:
        while True:
            name = f.readline().rstrip("\n")
            if not name:
                break
            seq = f.readline().rstrip("\n")
            plus = f.readline().rstrip("\n")
            qual = f.readline().rstrip("\n")
            if not plus.startswith("+"):
                raise ValueError(f"FASTQ parse error in {path}: missing '+' line.")
            yield (name, seq, qual)


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        raise ValueError("Hamming requires equal length strings.")
    return sum(1 for x, y in zip(a, b) if x != y)


@dataclass
class PrimerPattern:
    prefix: str
    umi_len: int
    suffix: str

    @staticmethod
    def from_primer_with_N(primer: str) -> "PrimerPattern":
        """
        Parse primer like: PREFIX + NNNNN... + SUFFIX, where N-run defines UMI length.
        Uses the LONGEST contiguous N-run as UMI (robust if someone accidentally has short N elsewhere).
        """
        p = primer.strip()
        if not p:
            raise ValueError("Empty primer string.")
        # find longest contiguous run of N
        best_start, best_len = -1, 0
        i = 0
        while i < len(p):
            if p[i].upper() == "N":
                j = i
                while j < len(p) and p[j].upper() == "N":
                    j += 1
                run_len = j - i
                if run_len > best_len:
                    best_start, best_len = i, run_len
                i = j
            else:
                i += 1
        if best_start < 0 or best_len <= 0:
            raise ValueError(f"Primer has no contiguous N-run: {primer}")
        prefix = p[:best_start]
        suffix = p[best_start + best_len :]
        return PrimerPattern(prefix=prefix, umi_len=best_len, suffix=suffix)


def load_primers_csv(path: str) -> Tuple[str, str]:
    """
    Expect columns: f, r (forward primer template, reverse primer template),
    each contains N-run for UMI.
    """
    with open(path, "rt", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    if not rows:
        raise ValueError("UMI primers csv is empty.")
    row = rows[0]
    if "f" not in row or "r" not in row:
        raise ValueError("UMI primers csv must have columns: f, r")
    return row["f"].strip(), row["r"].strip()


@dataclass
class MatchResult:
    ok: bool
    reason: str
    prefix_pos: Optional[int] = None
    shift_used: Optional[int] = None
    suffix_mism: Optional[int] = None
    umi: Optional[str] = None


def strict_match(seq: str, pat: PrimerPattern, mismatch: int = 0, shift: int = 0) -> MatchResult:
    """
    Strict stage3-like match with optional mismatch & shift:
    require prefix at pos=0+shift, suffix at pos=shift+len(prefix)+umi_len.
    """
    min_len = shift + len(pat.prefix) + pat.umi_len + len(pat.suffix)
    if len(seq) < min_len:
        return MatchResult(False, "too_short")

    p0 = shift
    if seq[p0 : p0 + len(pat.prefix)] != pat.prefix:
        return MatchResult(False, "missing_prefix_left")

    umi_start = p0 + len(pat.prefix)
    umi = seq[umi_start : umi_start + pat.umi_len]
    if any(c not in "ACGTacgt" for c in umi):
        return MatchResult(False, "invalid_umi_chars")

    suf_start = umi_start + pat.umi_len
    suf_obs = seq[suf_start : suf_start + len(pat.suffix)]
    suf_m = hamming(suf_obs.upper(), pat.suffix.upper())
    if suf_m > mismatch:
        return MatchResult(False, "missing_anchor_right", prefix_pos=p0, shift_used=shift, suffix_mism=suf_m, umi=umi)

    return MatchResult(True, "ok", prefix_pos=p0, shift_used=shift, suffix_mism=suf_m, umi=umi)


def window_find_best(seq: str, pat: PrimerPattern, max_window: int, max_shift: int, max_mismatch: int) -> MatchResult:
    """
    "Gengsong-like" constrained search:
    - search prefix occurrences within first max_window bases
    - allow shift up to max_shift (prefix can start at 0..max_shift)
    - choose candidate with lowest suffix mism (<=max_mismatch); break ties by smallest prefix_pos
    """
    best: Optional[MatchResult] = None

    # Candidate prefix positions: 0..max_shift, plus any find() hits in the window.
    cand_positions = set(range(0, max_shift + 1))
    # find exact prefix in early window
    if pat.prefix:
        start = 0
        while True:
            pos = seq.find(pat.prefix, start, max_window)
            if pos == -1:
                break
            cand_positions.add(pos)
            start = pos + 1

    for p0 in sorted(cand_positions):
        # require enough length
        min_len = p0 + len(pat.prefix) + pat.umi_len + len(pat.suffix)
        if len(seq) < min_len:
            continue
        if seq[p0 : p0 + len(pat.prefix)] != pat.prefix:
            continue

        umi_start = p0 + len(pat.prefix)
        umi = seq[umi_start : umi_start + pat.umi_len]
        if any(c not in "ACGTacgt" for c in umi):
            continue

        suf_start = umi_start + pat.umi_len
        suf_obs = seq[suf_start : suf_start + len(pat.suffix)]
        suf_m = hamming(suf_obs.upper(), pat.suffix.upper())
        if suf_m > max_mismatch:
            continue

        mr = MatchResult(True, "ok", prefix_pos=p0, shift_used=p0, suffix_mism=suf_m, umi=umi)
        if best is None:
            best = mr
        else:
            # choose lower mism, then smaller pos
            if (mr.suffix_mism, mr.prefix_pos) < (best.suffix_mism, best.prefix_pos):
                best = mr

    if best is not None:
        return best

    # If nothing found, decide most informative reason: too_short vs missing_prefix vs missing_anchor
    # We'll approximate: if seq length too short for p0=0, call too_short; else missing_prefix_left.
    if len(seq) < (len(pat.prefix) + pat.umi_len + len(pat.suffix)):
        return MatchResult(False, "too_short")
    return MatchResult(False, "no_hit_in_window")


def diagnose_pair(
    r1_seq: str,
    r2_seq: str,
    f_pat: PrimerPattern,
    r_pat: PrimerPattern,
    mismatch: int,
    shift: int,
    window: int,
    search_mode: str,
) -> Tuple[MatchResult, MatchResult]:
    if search_mode == "strict":
        m1 = strict_match(r1_seq, f_pat, mismatch=mismatch, shift=shift)
        m2 = strict_match(r2_seq, r_pat, mismatch=mismatch, shift=shift)
        return m1, m2
    else:
        # constrained search
        m1 = window_find_best(r1_seq, f_pat, max_window=window, max_shift=shift, max_mismatch=mismatch)
        m2 = window_find_best(r2_seq, r_pat, max_window=window, max_shift=shift, max_mismatch=mismatch)
        return m1, m2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r1", required=True, help="R1 FASTQ (.fastq or .fastq.gz)")
    ap.add_argument("--r2", required=True, help="R2 FASTQ (.fastq or .fastq.gz)")
    ap.add_argument("--primers_csv", required=True, help="CSV with columns f,r (N-run denotes UMI)")
    ap.add_argument("--n", type=int, default=20000, help="Number of read pairs to sample from start")
    ap.add_argument("--mismatch", type=int, default=0, help="Allowed hamming mismatches on suffix/anchor")
    ap.add_argument("--shift", type=int, default=0, help="Shift for strict mode; max_shift for search mode")
    ap.add_argument("--window", type=int, default=60, help="Search window size for prefix in search mode")
    ap.add_argument("--mode", choices=["strict", "search"], default="strict", help="strict vs constrained-search")
    ap.add_argument("--try_r2_revcomp", action="store_true", help="Also test R2 primer as reverse-complement; compare hit rates")
    ap.add_argument("--out_json", default="", help="Write summary JSON to this path (optional)")
    args = ap.parse_args()

    f_primer, r_primer = load_primers_csv(args.primers_csv)
    f_pat = PrimerPattern.from_primer_with_N(f_primer)
    r_pat = PrimerPattern.from_primer_with_N(r_primer)
    r_pat_rc = PrimerPattern.from_primer_with_N(revcomp(r_primer))

    print("== Parsed primer patterns ==")
    print(f"F: prefix_len={len(f_pat.prefix)} umi_len={f_pat.umi_len} suffix_len={len(f_pat.suffix)}")
    print(f"R: prefix_len={len(r_pat.prefix)} umi_len={r_pat.umi_len} suffix_len={len(r_pat.suffix)}")
    if args.try_r2_revcomp:
        print(f"R(rc): prefix_len={len(r_pat_rc.prefix)} umi_len={r_pat_rc.umi_len} suffix_len={len(r_pat_rc.suffix)}")
    print()

    counters = Counter()
    prefix_pos_hist_r1 = Counter()
    prefix_pos_hist_r2 = Counter()

    # optional A/B: R2 primer vs R2 primer revcomp
    ab_counts = Counter()

    it1 = fastq_iter(args.r1)
    it2 = fastq_iter(args.r2)

    total = 0
    for (n1, s1, q1), (n2, s2, q2) in zip(it1, it2):
        total += 1
        if total > args.n:
            break

        # baseline: using provided r primer
        m1, m2 = diagnose_pair(
            s1, s2, f_pat, r_pat,
            mismatch=args.mismatch, shift=args.shift, window=args.window, search_mode=args.mode
        )
        pair_ok = m1.ok and m2.ok
        if pair_ok:
            counters["pair_ok"] += 1
        else:
            counters["pair_fail"] += 1
            # pair-level reason: prioritize prefix missing / too_short / anchor missing / other
            # mimic your stage3 "pair-level" notion roughly:
            if (m1.reason == "too_short") or (m2.reason == "too_short"):
                counters["pair_too_short"] += 1
            elif (m1.reason in ("missing_prefix_left", "no_hit_in_window")) or (m2.reason in ("missing_prefix_left", "no_hit_in_window")):
                counters["pair_missing_prefix_left"] += 1
            elif (m1.reason == "missing_anchor_right") or (m2.reason == "missing_anchor_right"):
                counters["pair_missing_anchor_right"] += 1
            elif (m1.reason == "invalid_umi_chars") or (m2.reason == "invalid_umi_chars"):
                counters["pair_invalid_umi_chars"] += 1
            else:
                counters[f"pair_other_{m1.reason}_{m2.reason}"] += 1

        if m1.prefix_pos is not None:
            prefix_pos_hist_r1[m1.prefix_pos] += 1
        if m2.prefix_pos is not None:
            prefix_pos_hist_r2[m2.prefix_pos] += 1

        if args.try_r2_revcomp:
            # compare which R2 primer makes more sense for this read
            _, m2_rc = diagnose_pair(
                s1, s2, f_pat, r_pat_rc,
                mismatch=args.mismatch, shift=args.shift, window=args.window, search_mode=args.mode
            )
            key = ("R2_ok" if m2.ok else "R2_fail", "R2rc_ok" if m2_rc.ok else "R2rc_fail")
            ab_counts["/".join(key)] += 1

    def pct(x: int) -> str:
        if total == 0:
            return "0.00%"
        return f"{100.0 * x / total:.2f}%"

    print("== Summary ==")
    print(f"Total pairs checked: {total}")
    print(f"Pair OK:   {counters['pair_ok']} ({pct(counters['pair_ok'])})")
    print(f"Pair FAIL: {counters['pair_fail']} ({pct(counters['pair_fail'])})")
    print()
    print("Top pair-level fail buckets:")
    for k in ["pair_missing_anchor_right", "pair_missing_prefix_left", "pair_too_short", "pair_invalid_umi_chars"]:
        print(f"  - {k}: {counters[k]} ({pct(counters[k])})")
    print()

    def print_hist(title: str, hist: Counter, topk: int = 12):
        if not hist:
            print(f"{title}: (no prefix positions recorded)")
            return
        print(title)
        for pos, c in hist.most_common(topk):
            print(f"  pos={pos:<3d}  count={c:<7d}  ({100.0*c/max(1,sum(hist.values())):.2f}% of matched)")
        print()

    print_hist("R1 prefix position histogram (among reads where we recorded a prefix_pos):", prefix_pos_hist_r1)
    print_hist("R2 prefix position histogram (among reads where we recorded a prefix_pos):", prefix_pos_hist_r2)

    if args.try_r2_revcomp:
        print("== R2 primer vs R2 primer reverse-complement A/B ==")
        for k, v in ab_counts.most_common():
            print(f"  {k}: {v} ({pct(v)})")
        print()

    out = {
        "total_pairs": total,
        "counters": dict(counters),
        "r1_prefix_pos_hist": dict(prefix_pos_hist_r1),
        "r2_prefix_pos_hist": dict(prefix_pos_hist_r2),
        "r2_ab": dict(ab_counts),
        "params": {
            "mode": args.mode,
            "mismatch": args.mismatch,
            "shift": args.shift,
            "window": args.window,
            "n": args.n,
        },
        "primer_patterns": {
            "f": {"prefix_len": len(f_pat.prefix), "umi_len": f_pat.umi_len, "suffix_len": len(f_pat.suffix)},
            "r": {"prefix_len": len(r_pat.prefix), "umi_len": r_pat.umi_len, "suffix_len": len(r_pat.suffix)},
        }
    }
    if args.try_r2_revcomp:
        out["primer_patterns"]["r_revcomp"] = {"prefix_len": len(r_pat_rc.prefix), "umi_len": r_pat_rc.umi_len, "suffix_len": len(r_pat_rc.suffix)}

    if args.out_json:
        with open(args.out_json, "wt") as f:
            json.dump(out, f, indent=2)
        print(f"Wrote JSON summary to: {args.out_json}")


if __name__ == "__main__":
    main()