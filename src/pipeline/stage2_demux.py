# src/pipeline/stage2_demux.py
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict

from src.demux.demux_inline_index import demux_paired_inline_index, DemuxResult
from src.utils.fastq_io import FastqFormatError


def _write_json(path: Path, obj: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


def stage2_demux(
    *,
    exp: str,
    run_tag: str,
    r1: Path,
    r2: Path,
    multiplex_csv: Path,
    results_exp_dir: Path,
) -> Dict[str, Any]:
    """
    Stage 2:
    - inline index 分群 + strip index
    - unmatched 单独输出
    - 坏记录：立即报错停机（打印错误类型 + 首个出错 read id + 记录号）
    - 输出 metrics_demux.json
    - unmatched_rate 阈值提示/ABNORMAL 询问是否停机
    """
    outdir = results_exp_dir / "demultiplexing" / run_tag
    outdir.mkdir(parents=True, exist_ok=True)

    metrics_path = outdir / "metrics_demux.json"

    try:
        res: DemuxResult = demux_paired_inline_index(
            r1=r1,
            r2=r2,
            multiplex_csv=multiplex_csv,
            outdir=outdir,
            compress=True,
        )
    except FastqFormatError as e:
        # 你要求：错误类型 + 首个出错 read id + 行号/记录号，并写入 metrics 的 error 字段
        err = {
            "error_type": "FastqFormatError",
            "message": str(e),
            "first_bad_read_id": getattr(e, "read_id", None),
            "record_no": getattr(e, "record_no", None),
        }
        _write_json(metrics_path, {"stage": "stage2_demux", "status": "failed", "error": err})
        print("\nERROR (Stage 2 Demultiplex):", err["message"])
        print("First bad read id:", err["first_bad_read_id"])
        print("Record number:", err["record_no"])
        raise

    total = res.total_pairs
    unmatched = res.unmatched_pairs
    unmatched_rate = (unmatched / total) if total else 0.0

    metrics: Dict[str, Any] = {
        "stage": "stage2_demux",
        "exp": exp,
        "run_tag": run_tag,
        "status": "success",
        "inputs": {
            "r1": str(r1),
            "r2": str(r2),
            "multiplex_csv": str(multiplex_csv),
        },
        "outputs": {
            "outdir": str(outdir),
            "unmatched_r1": str(outdir / "unmatched_R1.fastq.gz"),
            "unmatched_r2": str(outdir / "unmatched_R2.fastq.gz"),
        },
        "unmatched_rate_definition": "unmatched pairs / total pairs",
        "total_pairs": total,
        "unmatched_pairs": unmatched,
        "unmatched_rate": unmatched_rate,
        "per_population_pairs": res.per_population_pairs,
        "index_lengths": res.index_lengths,
    }

    print("\n[Stage 2] Demultiplex metrics:")
    print(f"  total_pairs:     {total:,}")
    print(f"  unmatched_pairs: {unmatched:,}")
    print(f"  unmatched_rate:  {unmatched_rate:.3%}")

    # 阈值规则（你定义的）
    decision: Dict[str, Any] = {"abnormal": False, "stop_prompt": False}
    if unmatched_rate > 0.50:
        print("\nABNORMAL: unmatched_rate > 50%.")
        print("Recommended: STOP and check index/anchor/orientation/mismatch tolerance (likely CSV wrong or direction wrong).")
        ans = input("Stop now? (y/N) ").strip().lower()
        decision.update({"abnormal": True, "stop_prompt": True, "user_input": ans})
        if ans == "y":
            decision["stopped"] = True
            metrics["decision"] = decision
            _write_json(metrics_path, metrics)
            return {"status": "aborted", "metrics_path": str(metrics_path)}
        print("Proceeding despite ABNORMAL (user chose not to stop).")
    elif unmatched_rate > 0.35:
        print("\nSTRONG WARNING: unmatched_rate > 35%. Please inspect index/anchor/orientation.")
    elif unmatched_rate > 0.20:
        print("\nWARNING: unmatched_rate > 20%. Please inspect index/anchor/orientation.")

    metrics["decision"] = decision
    _write_json(metrics_path, metrics)

    return {"status": "success", "metrics_path": str(metrics_path)}