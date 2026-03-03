# src/qc/io/paths.py
from __future__ import annotations
from pathlib import Path

# 兼容层：QC 子系统继续 import qc.io.paths，
# 但真正实现统一走 src/utils/paths.py
from src.utils.paths import repo_root, results_dir, input_dir  # noqa: F401