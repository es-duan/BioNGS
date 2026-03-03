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
    """
    Open a text file that may be gzipped.

    If the file suffix is ".gz", it will be opened using gzip in text mode.
    Otherwise, it will be opened as a regular UTF-8 text file.

    Args:
        path (Path): Path to the input file.

    Returns:
        TextIO: A readable text file handle.
    """
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")

def validate_gzip_readable(path: Path) -> None:
    """
        Validate that a gzip-compressed file can be decompressed.

        If the file does not have a ".gz" suffix, this function does nothing.
        If it is a gzip file but cannot be decompressed, a ValueError is raised.

        Args:
            path (Path): Path to the input file.

        Raises:
            ValueError: If the gzip file cannot be decompressed.
        """
    if path.suffix != ".gz":
        return
    try:
        with gzip.open(path, "rb") as f:
            f.read(1)
    except Exception as e:
        raise ValueError(f"gzip cannot be decompressed: {path} ({e})") from e

def count_fastq_records(path: Path) -> int:
    """
    Count the number of records in a FASTQ file.

    A FASTQ record consists of 4 lines. The total number of lines must
    therefore be divisible by 4. If not, the file is considered invalid.

    Args:
        path (Path): Path to the FASTQ file (plain text or gzipped).

    Returns:
        int: Number of FASTQ records.

    Raises:
        ValueError: If the total number of lines is not divisible by 4.
    """
    n_lines = 0
    with _open_maybe_gzip(path) as f:
        for _ in f:
            n_lines += 1
    if n_lines % 4 != 0:
        raise ValueError(f"FASTQ line count not divisible by 4: {path} (lines={n_lines})")
    return n_lines // 4

def validate_fastq_pair(r1: Path, r2: Path) -> Dict[str, Any]:
    """
    Validate a paired-end FASTQ input (R1 and R2).

    This function checks:
    - Both files exist.
    - Gzip files (if applicable) are readable.
    - R1 and R2 contain the same number of FASTQ records.

    Args:
        r1 (Path): Path to the R1 FASTQ file.
        r2 (Path): Path to the R2 FASTQ file.

    Returns:
        Dict[str, Any]: A dictionary containing the number of records
                        in R1 and R2.

    Raises:
        FileNotFoundError: If either file does not exist.
        ValueError: If gzip validation fails or record counts differ.
    """
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
    """
    Read a multiplex configuration CSV file.

    The CSV file must contain a header row. UTF-8 files with BOM
    are supported. All values are stripped of surrounding whitespace.

    Args:
        path (Path): Path to the multiplex CSV file.

    Returns:
        List[Dict[str, str]]: A list of rows as dictionaries.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the CSV file has no header.
    """
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

REQUIRED_MULTIPLEX_COLS = {"Population", "R1_index", "R2_index"}

def validate_multiplex(rows: List[Dict[str, str]]) -> Dict[str, Any]:
    """
    Validate multiplex configuration rows for dual-index (R1 + R2).

    Rules:
      - Must have at least one data row.
      - Must contain Population, R1_index, R2_index.
      - Population name must match POP_SAFE_RE.
      - (R1_index, R2_index) pair must be unique.
    """
    if len(rows) == 0:
        raise ValueError("multiplex csv is empty (no data rows)")

    cols = set(rows[0].keys())
    missing = REQUIRED_MULTIPLEX_COLS - cols
    if missing:
        raise ValueError(f"multiplex csv missing required columns: {sorted(missing)}")

    seen_pairs: Dict[Tuple[str, str], int] = {}
    pop_names: List[str] = []

    for i, r in enumerate(rows, start=1):
        pop = r.get("Population", "")
        r1 = r.get("R1_index", "")
        r2 = r.get("R2_index", "")

        if not pop:
            raise ValueError(f"multiplex csv row {i}: Population is empty")
        if not r1:
            raise ValueError(f"multiplex csv row {i}: R1_index is empty")
        if not r2:
            raise ValueError(f"multiplex csv row {i}: R2_index is empty")

        if not POP_SAFE_RE.match(pop):
            raise ValueError(
                f"illegal Population name '{pop}' (row {i}); allowed: {POP_SAFE_RE.pattern}"
            )

        key = (r1, r2)
        if key in seen_pairs:
            j = seen_pairs[key]
            raise ValueError(
                f"duplicate index pair {key} in rows {j} and {i}"
            )
        seen_pairs[key] = i

        pop_names.append(pop)

    return {
        "n_rows": len(rows),
        "unique_populations": sorted(set(pop_names)),
        "n_populations": len(set(pop_names)),
        "index_key_mode": "R1+R2",
    }

def stage0_validate(
    *,
    r1: Path,
    r2: Path,
    multiplex_csv: Path,
) -> Dict[str, Any]:
    """
       Perform Stage 0 input validation for the pipeline.

       This stage validates:
       - Paired-end FASTQ integrity and consistency.
       - Multiplex configuration file structure and logic.

       Args:
           r1 (Path): Path to R1 FASTQ.
           r2 (Path): Path to R2 FASTQ.
           multiplex_csv (Path): Path to multiplex configuration CSV.

       Returns:
           Dict[str, Any]: A structured validation summary containing:
               - FASTQ statistics
               - Multiplex configuration statistics
               - Overall status flag ("ok")

       Raises:
           FileNotFoundError: If required files are missing.
           ValueError: If validation rules are violated.
       """
    pair_stats = validate_fastq_pair(r1, r2)
    rows = read_multiplex_csv(multiplex_csv)
    mux_stats = validate_multiplex(rows)

    return {
        "ok": True,
        "fastq": {"r1": str(r1), "r2": str(r2), **pair_stats},
        "multiplex_csv": {"path": str(multiplex_csv), **mux_stats},
    }

import json
from datetime import datetime
from typing import Dict, Any

def print_stage0_summary(stage0_result: Dict[str, Any], *, as_json: bool = False) -> None:
    """
    Print Stage 0 validation details to terminal.

    Args:
        stage0_result: dict returned by stage0_validate()
        as_json: if True, also print the raw JSON summary (pretty-printed)
    """
    fastq = stage0_result.get("fastq", {})
    mux = stage0_result.get("multiplex_csv", {})

    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    print("========== Stage 0 Summary ==========")
    print(f"Time         : {ts}")
    print(f"Status       : {'OK' if stage0_result.get('ok') else 'FAILED'}")
    print("-------------------------------------")
    print("FASTQ Pair")
    print(f"  R1          : {fastq.get('r1', '')}")
    print(f"  R2          : {fastq.get('r2', '')}")
    print(f"  R1 records  : {fastq.get('r1_records', 'NA')}")
    print(f"  R2 records  : {fastq.get('r2_records', 'NA')}")
    print("-------------------------------------")
    print("Multiplex CSV")
    print(f"  Path        : {mux.get('path', '')}")
    print(f"  Rows        : {mux.get('n_rows', 'NA')}")
    print(f"  Populations : {mux.get('n_populations', 'NA')}")

    pops = mux.get("unique_populations", [])
    if isinstance(pops, list):
        show = ", ".join(pops[:12])
        suffix = "" if len(pops) <= 12 else f" ... (+{len(pops)-12} more)"
        print(f"  Names       : {show}{suffix}")

    print(f"  Index mode  : {mux.get('index_key_mode', 'NA')}")
    print("=====================================\n")

    if as_json:
        print(json.dumps(stage0_result, indent=2, ensure_ascii=False))
        print()