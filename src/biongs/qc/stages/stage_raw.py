# src/qc/stages/stage_raw.py
from __future__ import annotations
from pathlib import Path
from typing import List, Tuple, Optional

from biongs.qc.io.paths import input_dir
from biongs.qc.io.manifest import Pair

def build_raw_pairs(
    experiment: str,
    gw_name: Optional[str],
    raw_r1: Optional[str],
    raw_r2: Optional[str],
    raw_subdir: str = "example_fastq",
) -> List[Pair]:
    """
    Return a single Pair: ("RAW", r1_path, r2_path)

    Prefer explicit raw_r1/raw_r2. Otherwise use gw_name + default layout:
        input_data/{exp}/{raw_subdir}/{gw_name}_R1_001.fastq
        input_data/{exp}/{raw_subdir}/{gw_name}_R2_001.fastq
    """
    if raw_r1 and raw_r2:
        r1 = Path(raw_r1).expanduser().resolve()
        r2 = Path(raw_r2).expanduser().resolve()
    else:
        if not gw_name:
            raise ValueError("RAW stage needs either (--raw_r1 and --raw_r2) OR --gw_name.")
        base = input_dir(experiment) / raw_subdir
        r1 = (base / f"{gw_name}_R1_001.fastq").resolve()
        r2 = (base / f"{gw_name}_R2_001.fastq").resolve()

    if not r1.exists():
        raise FileNotFoundError(f"RAW R1 not found: {r1}")
    if not r2.exists():
        raise FileNotFoundError(f"RAW R2 not found: {r2}")

    return [("RAW", r1, r2)]