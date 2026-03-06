# src/pipeline/stage1_raw_qc.py
from __future__ import annotations

import json
import subprocess
import zipfile
from pathlib import Path
from typing import Any, Dict

from primerforge.qc.overview_raw import compute_raw_overview, write_overview_report

SHORT_LEN_RAW = 150  # Stage 1 固定 150bp


def _run_cmd(cmd: list[str]) -> None:
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(
            "Command failed:\n"
            f"  cmd: {' '.join(cmd)}\n"
            f"  stdout:\n{p.stdout}\n"
            f"  stderr:\n{p.stderr}\n"
        )


def _safe_write_json(path: Path, obj: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


def _unzip_fastqc(zip_path: Path, outdir: Path) -> Path:
    """
    解包 FastQC 的 zip 到 outdir。
    返回解包后的 *_fastqc/ 目录路径（其中包含 fastqc_data.txt, summary.txt, Images/ 等）。
    """
    outdir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "r") as z:
        z.extractall(outdir)

    # FastQC zip 内部顶层一般是 "<basename>_fastqc/"
    dirs = [p for p in outdir.iterdir() if p.is_dir() and p.name.endswith("_fastqc")]
    if len(dirs) == 1:
        return dirs[0]

    # 兜底：找 fastqc_data.txt
    candidates = list(outdir.rglob("fastqc_data.txt"))
    if not candidates:
        raise FileNotFoundError(f"fastqc_data.txt not found after unzip: {zip_path}")
    return candidates[0].parent


def stage1_raw_qc(
    *,
    exp: str,
    run_tag: str,
    r1: Path,
    r2: Path,
    results_exp_dir: Path,
) -> Dict[str, Any]:
    """
    Stage 1 Raw QC:
    - qc_details: 跑 FastQC 并解包到 results/<exp>/qc_details/00_raw/
    - qc_overview: 不使用 FastQC；直接读 FASTQ 计算 Q30(base-level)、length distribution、short_rate(<150)，输出 txt+图
    - 判定规则：
        short_rate>30% -> WARNING (print)
        Q30<70% -> WARNING (print)
        short_rate>30% AND Q30<70% -> ABNORMAL -> Stop now? (y/N)
        都没触发 -> print 满足要求
    """
    qc_details_dir = results_exp_dir / "qc_details" / "00_raw" / run_tag
    qc_overview_dir = results_exp_dir / "qc_overview" / "00_raw" / run_tag
    qc_details_dir.mkdir(parents=True, exist_ok=True)
    qc_overview_dir.mkdir(parents=True, exist_ok=True)

    # ---------- 1) FastQC (details layer) ----------
    _run_cmd(["fastqc", str(r1), str(r2), "-o", str(qc_details_dir)])

    # ---------- 2) unzip *_fastqc.zip into qc_details/00_raw ----------
    zips = list(qc_details_dir.glob("*_fastqc.zip"))
    if len(zips) < 2:
        raise FileNotFoundError(f"FastQC zip outputs not found in {qc_details_dir}")

    def pick_zip(tag: str) -> Path:
        hits = [p for p in zips if tag in p.name]
        if len(hits) == 1:
            return hits[0]
        # 兜底：包含原 fastq 名字
        for p in zips:
            if r1.name in p.name and tag == "R1":
                return p
            if r2.name in p.name and tag == "R2":
                return p
        # 再兜底：返回第一个
        return zips[0]

    z1 = pick_zip("R1")
    z2 = pick_zip("R2")
    d1 = _unzip_fastqc(z1, qc_details_dir)
    d2 = _unzip_fastqc(z2, qc_details_dir)

    # ---------- 3) qc_overview WITHOUT FastQC ----------
    # 直接读 FASTQ：base-level Q30、length histogram、short_rate(<150)
    metrics_r1 = compute_raw_overview(r1, short_len=SHORT_LEN_RAW)
    metrics_r2 = compute_raw_overview(r2, short_len=SHORT_LEN_RAW)

    write_overview_report(metrics_r1, qc_overview_dir / "R1")
    write_overview_report(metrics_r2, qc_overview_dir / "R2")

    # Stage1 合并口径：short_rate 取更坏端；Q30 取更坏端
    short_rate = max(metrics_r1["short_rate"], metrics_r2["short_rate"])
    q30_rate = min(metrics_r1["q30_rate"], metrics_r2["q30_rate"])

    # ---------- 4) metrics ----------
    metrics: Dict[str, Any] = {
        "stage": "stage1_raw_qc",
        "exp": exp,
        "run_tag": run_tag,
        "short_len_bp": SHORT_LEN_RAW,
        "definitions": {
            "short_rate": f"read-level; read is short if length < {SHORT_LEN_RAW}",
            "q30_rate": "base-level; fraction of bases with Phred >= 30",
        },
        "qc_details": {
            "dir": str(qc_details_dir),
            "r1_fastqc_dir": str(d1),
            "r2_fastqc_dir": str(d2),
        },
        "qc_overview": {
            "dir": str(qc_overview_dir),
            "r1_overview_dir": str(qc_overview_dir / "R1"),
            "r2_overview_dir": str(qc_overview_dir / "R2"),
        },
        "r1": {
            "total_reads": metrics_r1["total_reads"],
            "q30_rate": metrics_r1["q30_rate"],
            "short_rate_lt150": metrics_r1["short_rate"],
        },
        "r2": {
            "total_reads": metrics_r2["total_reads"],
            "q30_rate": metrics_r2["q30_rate"],
            "short_rate_lt150": metrics_r2["short_rate"],
        },
        "combined": {
            "short_rate_used": short_rate,
            "q30_rate_used": q30_rate,
        },
    }

    # ---------- 5) print + decision (WARNING/ABNORMAL/passed) ----------
    warn_short = short_rate > 0.30
    warn_q30 = q30_rate < 0.70

    print("\n[Stage 1] Raw QC metrics (computed from FASTQ):")

    print(f"  R1 short_rate(<150bp): {metrics_r1['short_rate']:.3%}")
    print(f"  R2 short_rate(<150bp): {metrics_r2['short_rate']:.3%}")

    print(f"  R1 Q30_rate(base-level): {metrics_r1['q30_rate']:.3%}")
    print(f"  R2 Q30_rate(base-level): {metrics_r2['q30_rate']:.3%}")

    print(f"\n  Final short_rate (max of R1/R2): {short_rate:.3%}")
    print(f"  Final Q30_rate (min of R1/R2):  {q30_rate:.3%}")

    if warn_short and warn_q30:
        print("\nABNORMAL: short_rate > 30% AND Q30 < 70% (Raw QC).")
        print("Recommended: stop and inspect raw data quality / sequencing issues.")
        ans = input("Stop now? (y/N) ").strip().lower()
        metrics["decision"] = {"abnormal": True, "stop_prompt": True, "user_input": ans}

        # 先写 metrics，再决定退出与否
        metrics_path = results_exp_dir / "metrics" / "metrics_raw_qc.json"
        _safe_write_json(metrics_path, metrics)

        if ans == "y":
            metrics["decision"]["stopped"] = True
            return {"status": "aborted", "metrics_path": str(metrics_path)}
        print("Proceeding despite ABNORMAL (user chose not to stop).")

    else:
        if warn_short:
            print("\nWARNING: short_rate > 30% (Raw QC).")
            print("Triggered by: max(R1_short_rate, R2_short_rate).")
            print("Please inspect qc_overview/00_raw and qc_details/00_raw. ")
        if warn_q30:
            print("\nWARNING: Q30 < 70% (Raw QC). Please inspect qc_overview/qc_details.")
        if (not warn_short) and (not warn_q30):
            print("\nRaw QC check passed.")
            print("Data meets quality requirement.")

        metrics["decision"] = {"abnormal": False, "stop_prompt": False}

    metrics_path = results_exp_dir / "metrics" / "metrics_raw_qc.json"
    _safe_write_json(metrics_path, metrics)

    return {"status": "success", "metrics_path": str(metrics_path)}