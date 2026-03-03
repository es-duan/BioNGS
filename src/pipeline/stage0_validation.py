# src/pipeline/stage0_validation.py
from __future__ import annotations

import csv
import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional

POP_SAFE_RE = re.compile(r"^[A-Za-z0-9._-]+$")

REQUIRED_MULTIPLEX_COLS = {"Population", "R1_index"}  # R2_index 可选（看你数据是否有）

def _open_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")

def validate_gzip_readable(path: Path) -> None:
    # 允许非 gz；若是 gz 必须可解
    if path.suffix != ".gz":
        return
    try:
        with gzip.open(path, "rb") as f:
            f.read(1)
    except Exception as e:
        raise ValueError(f"gzip cannot be decompressed: {path} ({e})") from e

def count_fastq_records(path: Path) -> int:
    """
    记录数 = FASTQ 行数 / 4
    若行数不是4的倍数，直接报错（属于坏输入）
    """
    n_lines = 0
    with _open_maybe_gzip(path) as f:
        for _ in f:
            n_lines += 1
    if n_lines % 4 != 0:
        raise ValueError(f"FASTQ line count not divisible by 4: {path} (lines={n_lines})")
    return n_lines // 4

def validate_fastq_pair(r1: Path, r2: Path) -> Dict[str, Any]:
    if not r1.exists():
        raise FileNotFoundError(f"R1 not found: {r1}")
    if not r2.exists():
        raise FileNotFoundError(f"R2 not found: {r2}")

    validate_gzip_readable(r1)
    validate_gzip_readable(r2)

    n1 = count_fastq_records(r1)
    n2 = count_fastq_records(r2)
    if n1 != n2:
        raise ValueError(f"R1/R2 record counts differ: R1={n1}, R2={n2}")

    return {"r1_records": n1, "r2_records": n2}

def read_multiplex_csv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"multiplex csv not found: {path}")
    rows: List[Dict[str, str]] = []
    # 兼容带 BOM 的 csv
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"multiplex csv has no header: {path}")
        for r in reader:
            rows.append({k: (v or "").strip() for k, v in r.items()})
    return rows

def validate_multiplex(rows: List[Dict[str, str]]) -> Dict[str, Any]:
    if len(rows) == 0:
        raise ValueError("multiplex csv is empty (no data rows)")

    cols = set(rows[0].keys())
    missing = REQUIRED_MULTIPLEX_COLS - cols
    if missing:
        raise ValueError(f"multiplex csv missing required columns: {sorted(missing)}")

    # index 不重复（以 R1_index 为主；若你后续需要 R2_index 也可加）
    seen_r1: Dict[str, int] = {}
    pop_names: List[str] = []

    for i, r in enumerate(rows, start=1):
        pop = r.get("Population", "")
        idx1 = r.get("R1_index", "")

        if not pop:
            raise ValueError(f"multiplex csv row {i}: Population is empty")
        if not idx1:
            raise ValueError(f"multiplex csv row {i}: R1_index is empty")

        if not POP_SAFE_RE.match(pop):
            raise ValueError(f"illegal Population name '{pop}' (row {i}); allowed: {POP_SAFE_RE.pattern}")

        if idx1 in seen_r1:
            j = seen_r1[idx1]
            raise ValueError(f"duplicate R1_index '{idx1}' in rows {j} and {i}")
        seen_r1[idx1] = i

        pop_names.append(pop)

    return {
        "n_rows": len(rows),
        "unique_populations": sorted(set(pop_names)),
        "n_populations": len(set(pop_names)),
    }

def stage0_validate(
    *,
    r1: Path,
    r2: Path,
    multiplex_csv: Path,
) -> Dict[str, Any]:
    pair_stats = validate_fastq_pair(r1, r2)
    rows = read_multiplex_csv(multiplex_csv)
    mux_stats = validate_multiplex(rows)

    return {
        "ok": True,
        "fastq": {"r1": str(r1), "r2": str(r2), **pair_stats},
        "multiplex_csv": {"path": str(multiplex_csv), **mux_stats},
    }