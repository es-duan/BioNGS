from __future__ import annotations
from pathlib import Path

def repo_root() -> Path:
    # src/primerforge/utils/paths.py -> parents[3] == repo root
    return Path(__file__).resolve().parents[3]

def results_dir(experiment: str) -> Path:
    return repo_root() / "results" / experiment

def input_dir(experiment: str) -> Path:
    return repo_root() / "input_data" / experiment