#!/usr/bin/env python3
import argparse, gzip, csv, re
from collections import Counter

def open_txt(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "rt")

def iter_fastq(p):
    with open_txt(p) as f:
        while True:
            h = f.readline()
            if not h: return
            s = f.readline(); plus = f.readline(); q = f.readline()
            if not q: return
            yield h.strip(), s.strip(), plus.strip(), q.strip()

def rid(h):
    x = h.split()[0]
    if x.endswith("/1") or x.endswith("/2"):
        x = x[:-2]
    return x

def parse_10N(primer, umi_len=10):
    primer = primer.strip().upper()
    m = list(re.finditer(rf"N{{{umi_len}}}", primer))
    if len(m) != 1:
        raise ValueError(f"primer must contain exactly one N{{{umi_len}}} block: {primer}")
    st, ed = m[0].start(), m[0].end()
    if "N" in (primer[:st] + primer[ed:]):
        raise ValueError(f"primer has extra N outside UMI block: {primer}")
    before, after = primer[:st], primer[ed:]
    if not before or not after:
        raise ValueError(f"primer before/after empty: {primer}")
    return before, after

def load_primers(csv_path):
    rows = list(csv.DictReader(open(csv_path, "r")))
    f = rows[0]["f"]; r = rows[0]["r"]
    bf, af = parse_10N(f, 10)
    br, ar = parse_10N(r, 10)
    return (bf, af), (br, ar)

def extract(seq, qual, before, after, umi_len=10):
    seq = seq.upper()
    pos = seq.find(before)
    if pos < 0:
        return "missing_prefix"
    umi_start = pos + len(before)
    umi_end = umi_start + umi_len
    suf_start = umi_end
    suf_end = suf_start + len(after)
    if len(seq) < suf_end or len(qual) < suf_end:
        return "too_short"
    umi = seq[umi_start:umi_end]
    if not re.fullmatch(r"[ACGT]{10}", umi):
        return "invalid_umi"
    if seq[suf_start:suf_end] != after:
        return "missing_anchor"
    return "ok"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--primers", required=True)
    ap.add_argument("--n", type=int, default=200000)
    ap.add_argument("--minlen", type=int, default=150)
    args = ap.parse_args()

    (b1,a1), (b2,a2) = load_primers(args.primers)

    c_pair = Counter()
    c_r1 = Counter()
    c_r2 = Counter()

    it1 = iter_fastq(args.r1)
    it2 = iter_fastq(args.r2)

    total = 0
    short_pairs = 0
    usable = 0
    ok_pairs = 0

    for rec1, rec2 in zip(it1, it2):
        total += 1
        if total > args.n: break
        h1,s1,_,q1 = rec1
        h2,s2,_,q2 = rec2

        if rid(h1) != rid(h2):
            c_pair["read_id_mismatch"] += 1
            break

        # Stage2-style short filter (pair-level)
        if len(s1) < args.minlen or len(s2) < args.minlen:
            short_pairs += 1
            continue

        usable += 1
        r1 = extract(s1,q1,b1,a1)
        r2 = extract(s2,q2,b2,a2)
        c_r1[r1] += 1
        c_r2[r2] += 1

        if r1 == "ok" and r2 == "ok":
            ok_pairs += 1
            c_pair["ok"] += 1
        else:
            # pair-level reason (prioritize anchor/prefix)
            if "missing_anchor" in (r1,r2):
                c_pair["missing_anchor"] += 1
            elif "missing_prefix" in (r1,r2):
                c_pair["missing_prefix"] += 1
            elif "too_short" in (r1,r2):
                c_pair["too_short"] += 1
            elif "invalid_umi" in (r1,r2):
                c_pair["invalid_umi"] += 1
            else:
                c_pair["other_fail"] += 1

    print("=== UNIFIED DEBUG METRICS ===")
    print(f"Total pairs scanned: {total}")
    print(f"Short pairs (<{args.minlen} on either mate): {short_pairs} ({short_pairs/max(1,total):.2%})")
    print(f"Usable (non-short) pairs: {usable} ({usable/max(1,total):.2%})")
    if usable:
        print(f"OK pairs among usable: {ok_pairs} ({ok_pairs/usable:.2%})")
    print("\nPair-level fails (usable only):")
    for k,v in c_pair.most_common():
        if k=="ok": continue
        print(f"  {k:15s} {v:10d} ({v/max(1,usable):.2%})")

    print("\nR1 reasons (usable only):")
    for k,v in c_r1.most_common():
        print(f"  {k:15s} {v:10d} ({v/max(1,usable):.2%})")

    print("\nR2 reasons (usable only):")
    for k,v in c_r2.most_common():
        print(f"  {k:15s} {v:10d} ({v/max(1,usable):.2%})")

if __name__ == "__main__":
    main()