# src/qc/stages/stage_after_demux.py
from __future__ import annotations
import re
from pathlib import Path
from typing import List, Tuple

from src.qc.io.paths import results_dir
from src.qc.io.manifest import Pair

def _pop_sort_key(label: str) -> Tuple[int, str]:
    m = re.match(r"^[Pp](\d+)$", label.strip())
    return (int(m.group(1)), label) if m else (10**9, label)

def build_demux_pairs(experiment: str) -> List[Pair]:
    demux_dir = results_dir(experiment) / "demultiplexing"
    if not demux_dir.exists():
        raise FileNotFoundError(f"Demultiplexing dir not found: {demux_dir} (run Script2 first)")

    pairs: List[Pair] = []
    for pop_dir in demux_dir.glob("P*"):
        if not pop_dir.is_dir():
            continue
        label = pop_dir.name  # P1, P2...
        r1 = (pop_dir / f"{label}_R1.fastq").resolve()
        r2 = (pop_dir / f"{label}_R2.fastq").resolve()
        if r1.exists() and r2.exists():
            pairs.append((label, r1, r2))

    pairs.sort(key=lambda x: _pop_sort_key(x[0]))

    if not pairs:
        raise FileNotFoundError(
            f"No population fastqs found under: {demux_dir} (expected P*/P*_R1.fastq & P*/P*_R2.fastq)"
        )

    return pairs