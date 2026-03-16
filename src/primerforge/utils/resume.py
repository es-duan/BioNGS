from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict, Optional


def file_fingerprint(path: str | Path) -> Dict[str, Any]:
    p = Path(path)
    if not p.exists():
        return {"path": str(p), "exists": False}
    st = p.stat()
    return {"path": str(p), "exists": True, "size": st.st_size, "mtime": int(st.st_mtime)}


def _stable(d: Dict[str, Any]) -> Dict[str, Any]:
    # stable dict for comparison (order-independent)
    return json.loads(json.dumps(d, sort_keys=True))


class ResumeManager:
    """
    Stores per-stage completion markers under:
      <exp_dir>/.resume/<run_tag>/<stage>.json
    """
    def __init__(self, exp_dir: Path, run_tag: str):
        self.exp_dir = Path(exp_dir)
        self.run_tag = run_tag

    @property
    def resume_dir(self) -> Path:
        d = self.exp_dir / ".resume" / self.run_tag
        d.mkdir(parents=True, exist_ok=True)
        return d

    def record_path(self, stage: str) -> Path:
        return self.resume_dir / f"{stage}.json"

    def is_done(self, stage: str, signature: Dict[str, Any]) -> bool:
        fp = self.record_path(stage)
        if not fp.exists():
            return False
        try:
            saved = json.loads(fp.read_text())
        except Exception:
            return False
        return saved.get("done") is True and saved.get("signature") == _stable(signature)

    def mark_done(self, stage: str, signature: Dict[str, Any], extra: Optional[Dict[str, Any]] = None) -> None:
        payload: Dict[str, Any] = {
            "stage": stage,
            "done": True,
            "time": int(time.time()),
            "signature": _stable(signature),
        }
        if extra:
            payload["extra"] = extra
        self.record_path(stage).write_text(json.dumps(payload, indent=2, sort_keys=True))

    def load(self, stage: str) -> dict:
        fp = self.record_path(stage)
        if not fp.exists():
            return {}
        try:
            return json.loads(fp.read_text())
        except Exception:
            return {}

    def clear_stage(self, stage: str) -> None:
        fp = self.record_path(stage)
        if fp.exists():
            fp.unlink()
