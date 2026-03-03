# src/demux/demux_inline_index.py
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple

from src.utils.fastq_io import (
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


def load_multiplex_csv(multiplex_csv: Path) -> Dict[str, str]:
    """
    读取 multiplexing csv，返回 index -> population
    支持：
      - index
      - R1_index
    population 必须存在
    不区分大小写
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

        # index 列自动识别
        if "index" in header_map:
            idx_col = header_map["index"]
        elif "r1_index" in header_map:
            idx_col = header_map["r1_index"]
        else:
            raise RuntimeError(
                f"multiplex csv must contain column: index or R1_index. Found: {reader.fieldnames}"
            )

        mapping: Dict[str, str] = {}
        for row in reader:
            idx = (row.get(idx_col) or "").strip()
            pop = (row.get(pop_col) or "").strip()
            if not idx or not pop:
                continue
            if idx in mapping:
                raise RuntimeError(f"Duplicate index in multiplex csv: {idx}")
            mapping[idx] = pop

    if not mapping:
        raise RuntimeError("multiplex csv produced empty index->population mapping")

    print(f"[Stage 2] Loaded {len(mapping)} index->population mappings.")
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
) -> DemuxResult:
    """
    Stage2 核心：
    - 逐对读取 R1/R2
    - 用 R1 前缀 inline index 精确匹配（不做错配容忍）
    - matched -> population 文件夹
    - unmatched -> unmatched_R1/R2
    - demux 阶段就 strip index（matched 一定 strip；unmatched 也 strip 用统一 index_len）
    - 遇到坏记录：立即报错停机（FastqFormatError）
    """
    idx2pop = load_multiplex_csv(multiplex_csv)
    index_lengths = {idx: len(idx) for idx in idx2pop.keys()}
    # 常见情况：全等长
    unique_lens = sorted(set(index_lengths.values()))
    index_len_for_unmatched = max(unique_lens)

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

            # inline index 从 R1 头部取（精确匹配）
            # 用“所有 index 的最大长度”先取前缀，再精确比对（避免不同长度 index 时切片错位）
            prefix = rec1.seq[:index_len_for_unmatched]

            matched_idx = None
            matched_pop = None
            # 精确：尝试所有 idx（规模一般很小）
            for idx, pop in idx2pop.items():
                if rec1.seq.startswith(idx):
                    matched_idx = idx
                    matched_pop = pop
                    break

            if matched_pop is None:
                unmatched_pairs += 1

                # unmatched：仍 strip（用统一 index_len_for_unmatched）
                rec1s = strip_inline_index(rec1, index_len_for_unmatched)
                rec2s = strip_inline_index(rec2, index_len_for_unmatched)

                write_fastq_record(unmatched1, rec1s)
                write_fastq_record(unmatched2, rec2s)
                continue

            # matched：按 matched_idx 的长度 strip（更严格）
            strip_len = len(matched_idx)
            rec1s = strip_inline_index(rec1, strip_len)
            rec2s = strip_inline_index(rec2, strip_len)

            h1, h2 = get_handles(matched_pop)
            write_fastq_record(h1, rec1s)
            write_fastq_record(h2, rec2s)

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