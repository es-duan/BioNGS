# src/pipeline/run_manifest.py
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional

def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()

def new_run_id() -> str:
    # 时间戳 + pid，足够区分同一天多次运行
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"{ts}_pid{os.getpid()}"

def write_json(path: Path, obj: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")
    tmp.replace(path)

def load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))

def init_run_manifest(
    manifest_path: Path,
    *,
    exp: str,
    run_tag: str,
    inputs: Dict[str, Any],
) -> Dict[str, Any]:
    obj: Dict[str, Any] = {
        "run_id": new_run_id(),
        "exp": exp,
        "run_tag": run_tag,
        "created_utc": utc_now_iso(),
        "inputs": inputs,
        "stages": {},
        "status": "running",
    }
    write_json(manifest_path, obj)
    return obj

def update_stage(
    manifest_path: Path,
    stage_name: str,
    payload: Dict[str, Any],
) -> None:
    obj = load_json(manifest_path)
    obj.setdefault("stages", {})
    obj["stages"][stage_name] = payload
    write_json(manifest_path, obj)

def finalize_manifest(manifest_path: Path, status: str) -> None:
    obj = load_json(manifest_path)
    obj["status"] = status
    obj["finished_utc"] = utc_now_iso()
    write_json(manifest_path, obj)