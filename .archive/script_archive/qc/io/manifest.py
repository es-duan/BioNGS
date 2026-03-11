# src/qc/io/manifest.py
from __future__ import annotations
import csv
from pathlib import Path
from typing import List, Tuple

Pair = Tuple[str, Path, Path]  # (label, r1, r2)

def write_manifest(pairs: List[Pair], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["label", "r1", "r2"])
        w.writeheader()
        for label, r1, r2 in pairs:
            w.writerow({"label": label, "r1": str(r1), "r2": str(r2)})

def read_manifest(out_csv: Path) -> List[Pair]:
    pairs: List[Pair] = []
    with out_csv.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for r in reader:
            label = (r.get("label") or "").strip()
            r1 = Path((r.get("r1") or "").strip()).expanduser().resolve()
            r2 = Path((r.get("r2") or "").strip()).expanduser().resolve()
            pairs.append((label, r1, r2))
    return pairs