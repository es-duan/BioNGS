# src/demux/demux_inline_index.py
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple, Optional

from primerforge.utils.fastq_io import (
    FastqFormatError,
    FastqRecord,
    normalize_read_id,
    open_text_maybe_gz,
    read_fastq,
    write_fastq_record,
)


@dataclass
class DemuxResult:
    total_pairs: int
    unmatched_pairs: int
    per_population_pairs: Dict[str, int]
    index_lengths: Dict[str, int]


def load_multiplex_csv(multiplex_csv: Path) -> Dict[Tuple[str, str], str]:
    """
    读取 multiplexing csv，返回 (R1_index, R2_index) -> population
    要求列：Population, R1_index, R2_index（大小写不敏感）
    """
    with multiplex_csv.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise RuntimeError("multiplex csv has no header row")

        # normalize headers
        header_map = {name.lower().strip(): name for name in reader.fieldnames}

        # population 必须有
        if "population" not in header_map:
            raise RuntimeError(
                f"multiplex csv must contain column: Population. Found: {reader.fieldnames}"
            )
        pop_col = header_map["population"]

        # 双 index 必须有
        if "r1_index" not in header_map or "r2_index" not in header_map:
            raise RuntimeError(
                f"multiplex csv must contain columns: R1_index and R2_index. Found: {reader.fieldnames}"
            )
        r1_col = header_map["r1_index"]
        r2_col = header_map["r2_index"]

        mapping: Dict[Tuple[str, str], str] = {}
        for row in reader:
            r1 = (row.get(r1_col) or "").strip().upper()
            r2 = (row.get(r2_col) or "").strip().upper()
            pop = (row.get(pop_col) or "").strip()
            if not r1 or not r2 or not pop:
                continue

            key = (r1, r2)
            if key in mapping:
                raise RuntimeError(f"Duplicate (R1_index,R2_index) in multiplex csv: {key}")
            mapping[key] = pop

    if not mapping:
        raise RuntimeError("multiplex csv produced empty (R1_index,R2_index)->population mapping")

    print(f"[Stage 2] Loaded {len(mapping)} (R1_index,R2_index)->population mappings.")
    return mapping


def strip_inline_index(rec: FastqRecord, n: int) -> FastqRecord:
    if n <= 0:
        return rec
    if len(rec.seq) < n:
        # 这属于“坏记录”（不足以读 index），按你的规范应报错停机
        rid = normalize_read_id(rec.header)
        raise FastqFormatError(
            f"Read shorter than required inline index length: len(seq)={len(rec.seq)} < {n}",
            read_id=rid,
            record_no=None,
        )
    return FastqRecord(
        header=rec.header,
        seq=rec.seq[n:],
        plus=rec.plus,
        qual=rec.qual[n:],
    )


def demux_paired_inline_index(
    *,
    r1: Path,
    r2: Path,
    multiplex_csv: Path,
    outdir: Path,
    compress: bool = True,
    strip_matched: bool = False
) -> DemuxResult:
    """
    Stage2 核心：
    - 逐对读取 R1/R2
    - 用 R1 前缀 inline index 精确匹配（不做错配容忍）
    - matched -> population 文件夹
    - unmatched -> unmatched_R1/R2
    - demux 只负责分群
    - 是否 strip index 由 strip_matched 控制（默认不 strip）
    - unmatched 永远原样输出（你现在就是这样做的）
    - 遇到坏记录：立即报错停机（FastqFormatError）
    """
    idx2pop = load_multiplex_csv(multiplex_csv)

    # 分别统计 R1 和 R2 的 index 长度（因为你文库是 R1/R2 都有 index）
    r1_lens = [len(r1_idx) for (r1_idx, _r2_idx) in idx2pop.keys()]
    r2_lens = [len(r2_idx) for (_r1_idx, r2_idx) in idx2pop.keys()]

    max_r1_len = max(r1_lens)
    max_r2_len = max(r2_lens)

    # 保留原字段，但改成双端长度信息（供 metrics 用）
    index_lengths = {
        "r1_unique_lengths": sorted(set(r1_lens)),
        "r2_unique_lengths": sorted(set(r2_lens)),
        "r1_max_len": max_r1_len,
        "r2_max_len": max_r2_len,
    }

    outdir.mkdir(parents=True, exist_ok=True)

    # output handles
    suffix = ".fastq.gz" if compress else ".fastq"
    handles: Dict[str, Tuple] = {}

    def get_handles(pop: str):
        if pop in handles:
            return handles[pop]
        pop_dir = outdir / pop
        pop_dir.mkdir(parents=True, exist_ok=True)
        f1 = pop_dir / f"{pop}_R1{suffix}"
        f2 = pop_dir / f"{pop}_R2{suffix}"
        h1 = open_text_maybe_gz(f1, "wt")
        h2 = open_text_maybe_gz(f2, "wt")
        handles[pop] = (h1, h2)
        return h1, h2

    unmatched1 = open_text_maybe_gz(outdir / f"unmatched_R1{suffix}", "wt")
    unmatched2 = open_text_maybe_gz(outdir / f"unmatched_R2{suffix}", "wt")

    per_pop_pairs: Dict[str, int] = {}
    total_pairs = 0
    unmatched_pairs = 0

    it1 = read_fastq(r1)
    it2 = read_fastq(r2)

    record_no = 0
    try:
        for rec1, rec2 in zip(it1, it2):
            record_no += 1
            total_pairs += 1

            id1 = normalize_read_id(rec1.header)
            id2 = normalize_read_id(rec2.header)
            if id1 != id2:
                raise FastqFormatError(
                    f"R1/R2 not paired: id mismatch (R1={id1}, R2={id2})",
                    read_id=id1,
                    record_no=record_no,
                )

            # dual inline index：要求 R1 和 R2 都匹配同一对 (R1_index, R2_index)
            r1_head = rec1.seq[:max_r1_len].upper()
            r2_head = rec2.seq[:max_r2_len].upper()

            matched_key: Optional[Tuple[str, str]] = None
            matched_pop: Optional[str] = None

            for (r1_idx, r2_idx), pop in idx2pop.items():
                if r1_head.startswith(r1_idx) and r2_head.startswith(r2_idx):
                    matched_key = (r1_idx, r2_idx)
                    matched_pop = pop
                    break

            if matched_pop is None:
                unmatched_pairs += 1

                # unmatched：原样输出（不要 strip，保留证据，方便你排查到底是 index 错还是测序错）
                write_fastq_record(unmatched1, rec1)
                write_fastq_record(unmatched2, rec2)
                continue

            # matched：默认不裁剪；如需裁剪 inline index，则打开 strip_matched
            if strip_matched:
                r1_idx, r2_idx = matched_key  # matched_key 在这里一定不是 None
                rec1_out = strip_inline_index(rec1, len(r1_idx))
                rec2_out = strip_inline_index(rec2, len(r2_idx))
            else:
                rec1_out = rec1
                rec2_out = rec2

            h1, h2 = get_handles(matched_pop)
            write_fastq_record(h1, rec1_out)
            write_fastq_record(h2, rec2_out)

            per_pop_pairs[matched_pop] = per_pop_pairs.get(matched_pop, 0) + 1

        # 检查 zip 是否提前结束（记录数不一致）
        # 如果一个迭代器还有剩余，zip 会截断；我们要把它当“坏记录/坏输入”报错停机
        try:
            next(it1)
            # R1 还有
            raise FastqFormatError("R1 has extra records beyond R2 (record count mismatch).", None, record_no + 1)
        except StopIteration:
            pass
        try:
            next(it2)
            # R2 还有
            raise FastqFormatError("R2 has extra records beyond R1 (record count mismatch).", None, record_no + 1)
        except StopIteration:
            pass

    finally:
        # close all handles
        for h1, h2 in handles.values():
            h1.close()
            h2.close()
        unmatched1.close()
        unmatched2.close()

    return DemuxResult(
        total_pairs=total_pairs,
        unmatched_pairs=unmatched_pairs,
        per_population_pairs=dict(sorted(per_pop_pairs.items())),
        index_lengths=index_lengths,
    )