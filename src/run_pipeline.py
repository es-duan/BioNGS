# src/run_pipeline.py
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from src.utils.paths import results_dir
from src.pipeline.run_manifest import (
    init_run_manifest,
    update_stage,
    finalize_manifest,
    write_json,
)
from src.pipeline.stage0_validation import stage0_validate
from src.pipeline.stage1_raw_qc import stage1_raw_qc
from src.pipeline.stage2_demux import stage2_demux
from src.pipeline.stage2_5_demux_qc import stage2_5_demux_qc
from src.pipeline.stage3_umi_extract import stage3_umi_extract
from src.pipeline.stage3_5_umi_quality import stage3_5_umi_quality
from src.pipeline.stage4_clean_qc import stage4_clean_qc
from src.pipeline.stage5_preprocess import stage5_preprocess
from src.pipeline.stage6_post_trim_qc import stage6_post_trim_qc

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="BioNGS Pipeline")
    p.add_argument("--exp", required=True, help="Experiment name")
    p.add_argument("--run_tag", required=True, help="Run tag (e.g. run1)")
    p.add_argument("--r1", required=True, help="Path to R1 FASTQ")
    p.add_argument("--r2", required=True, help="Path to R2 FASTQ")
    p.add_argument("--multiplex_csv", required=True, help="Multiplex CSV file")
    p.add_argument("--umi_primers_csv", required=True, help="UMI primers CSV (columns: f,r)")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    exp = args.exp
    run_tag = args.run_tag
    r1 = Path(args.r1).resolve()
    r2 = Path(args.r2).resolve()
    multiplex_csv = Path(args.multiplex_csv).resolve()
    umi_primers_csv = Path(args.umi_primers_csv).resolve()

    exp_dir = results_dir(exp)
    metrics_dir = exp_dir / "metrics"
    metrics_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = exp_dir / "run_manifest.json"
    init_run_manifest(
        manifest_path,
        exp=exp,
        run_tag=run_tag,
        inputs={
            "r1": str(r1),
            "r2": str(r2),
            "multiplex_csv": str(multiplex_csv),
            "umi_primers_csv": str(umi_primers_csv),
        },
    )

    print("\n=== BioNGS Pipeline Started ===")
    print(f"Experiment: {exp}")
    print(f"Run tag:    {run_tag}\n")

    # ---------- Stage 0 ----------
    print("Running Stage 0: Input validation...")
    try:
        stage0_result = stage0_validate(r1=r1, r2=r2, multiplex_csv=multiplex_csv)
        write_json(metrics_dir / "metrics_stage0_validation.json", stage0_result)
        update_stage(manifest_path, "stage0_validation", {"status": "success", "result": stage0_result})
        print("Stage 0 completed successfully.\n")
    except Exception as e:
        update_stage(manifest_path, "stage0_validation", {"status": "failed", "error": str(e)})
        finalize_manifest(manifest_path, status="failed")
        print("\nStage 0 FAILED.")
        print("Error:", str(e))
        sys.exit(1)

    # ---------- Stage 1 ----------
    print("Running Stage 1: Raw QC (FastQC details + custom overview)...")
    try:
        stage1_result = stage1_raw_qc(exp=exp, run_tag=run_tag, r1=r1, r2=r2, results_exp_dir=exp_dir)
        update_stage(manifest_path, "stage1_raw_qc", stage1_result)

        if stage1_result.get("status") == "aborted":
            finalize_manifest(manifest_path, status="aborted_stage1")
            print("\nPipeline aborted by user at Stage 1.")
            sys.exit(2)

        print("\nStage 1 completed successfully.\n")
    except Exception as e:
        update_stage(manifest_path, "stage1_raw_qc", {"status": "failed", "error": str(e)})
        finalize_manifest(manifest_path, status="failed")
        print("\nStage 1 FAILED.")
        print("Error:", str(e))
        sys.exit(1)

    # ---------- Stage 2 ----------
    print("Running Stage 2: Demultiplex (inline index + strip index)...")
    try:
        stage2_result = stage2_demux(
            exp=exp,
            run_tag=run_tag,
            r1=r1,
            r2=r2,
            multiplex_csv=multiplex_csv,
            results_exp_dir=exp_dir,
        )
        update_stage(manifest_path, "stage2_demux", stage2_result)

        if stage2_result.get("status") == "aborted":
            finalize_manifest(manifest_path, status="aborted_stage2")
            print("\nPipeline aborted by user at Stage 2.")
            sys.exit(2)

        print("\nStage 2 completed successfully.\n")
    except Exception as e:
        update_stage(manifest_path, "stage2_demux", {"status": "failed", "error": str(e)})
        finalize_manifest(manifest_path, status="failed")
        print("\nStage 2 FAILED.")
        print("Error:", str(e))
        sys.exit(1)

    # ---------- Stage 2.5 ----------
    print("Running Stage 2.5: Demux QC (per-pop FastQC details + custom overview)...")
    try:
        stage2_5_result = stage2_5_demux_qc(
            exp=exp,
            run_tag=run_tag,
            results_exp_dir=exp_dir,
        )
        update_stage(manifest_path, "stage2_5_demux_qc", stage2_5_result)

        if stage2_5_result.get("status") == "aborted":
            finalize_manifest(manifest_path, status="aborted_stage2_5")
            print("\nPipeline aborted by user at Stage 2.5.")
            sys.exit(2)

        print("\nStage 2.5 completed successfully.\n")

    except Exception as e:
        update_stage(manifest_path, "stage2_5_demux_qc", {"status": "failed", "error": str(e)})
        finalize_manifest(manifest_path, status="failed")
        print("\nStage 2.5 FAILED.")
        print("Error:", str(e))
        sys.exit(1)

    # ---------- Stage 3 ----------
    print("Running Stage 3: UMI extraction (trim primers+UMI, tag header)...")
    try:
        stage3_result = stage3_umi_extract(
            exp=exp,
            run_tag=run_tag,
            results_exp_dir=exp_dir,
            umi_primers_csv=umi_primers_csv,
        )
        update_stage(manifest_path, "stage3_umi_extract", stage3_result)

        if stage3_result.get("status") == "aborted":
            finalize_manifest(manifest_path, status="aborted_stage3")
            print("\nPipeline aborted by user at Stage 3.")
            sys.exit(2)

        print("\nStage 3 completed successfully.\n")
    except Exception as e:
        update_stage(manifest_path, "stage3_umi_extract", {"status": "failed", "error": str(e)})
        finalize_manifest(manifest_path, status="failed")
        print("\nStage 3 FAILED.")
        print("Error:", str(e))
        sys.exit(1)

    # ---------- Stage 3.5 ----------
    print("Running Stage 3.5: UMI quality summary (txt + terminal print)...")
    try:
        stage3_5_result = stage3_5_umi_quality(
            exp=exp,
            run_tag=run_tag,
            results_exp_dir=exp_dir,
            stage3_result=stage3_result,
        )
        update_stage(manifest_path, "stage3_5_umi_quality", stage3_5_result)
        print("\nStage 3.5 completed successfully.\n")
    except Exception as e:
        update_stage(manifest_path, "stage3_5_umi_quality", {"status": "failed", "error": str(e)})
        finalize_manifest(manifest_path, status="failed")
        print("\nStage 3.5 FAILED.")
        print("Error:", str(e))
        sys.exit(1)

    # ---------- Stage 4 ----------
    print("Running Stage 4: Clean QC after UMI/primer stripping (per-pop details + overview)...")
    try:
        stage4_result = stage4_clean_qc(exp=exp, run_tag=run_tag, results_exp_dir=exp_dir)
        update_stage(manifest_path, "stage4_clean_qc", stage4_result)
        write_json(metrics_dir / "metrics_stage4_clean_qc.json", stage4_result)
        print("\nStage 4 completed successfully.\n")
    except Exception as e:
        update_stage(manifest_path, "stage4_clean_qc", {"status": "failed", "error": str(e)})
        finalize_manifest(manifest_path, status="failed")
        print("\nStage 4 FAILED.")
        print("Error:", str(e))
        sys.exit(1)

    # ---------- Stage 5 ----------
    print("Running Stage 5: Preprocessing (fastp tail-trim + interactive min_len preview + length filtering)...")
    try:
        stage5_result = stage5_preprocess(
            exp=exp,
            run_tag=run_tag,
            results_exp_dir=exp_dir,
            metrics_dir=metrics_dir,
        )
        update_stage(manifest_path, "stage5_preprocess", stage5_result)

        if stage5_result.get("status") == "aborted":
            finalize_manifest(manifest_path, status="aborted_stage5")
            print("\nPipeline aborted by user at Stage 5.")
            sys.exit(2)

        print("\nStage 5 completed successfully.\n")
    except Exception as e:
        update_stage(manifest_path, "stage5_preprocess", {"status": "failed", "error": str(e)})
        finalize_manifest(manifest_path, status="failed")
        print("\nStage 5 FAILED.")
        print("Error:", str(e))
        sys.exit(1)

    # ---------- Stage 6 ----------
    print("Running Stage 6: Post-preprocess QC (FastQC details + overview txt + freq plot)...")
    try:
        trimmed_dir = Path(stage5_result["trimmed_dir"])  # 你按自己stage5的返回字段改这个key
        stage6_result = stage6_post_trim_qc(
            exp=exp,
            run_tag=run_tag,
            trimmed_dir=trimmed_dir,
            results_exp_dir=exp_dir,
        )
        update_stage(manifest_path, "stage6_post_trim_qc", stage6_result)
        print("\nStage 6 completed successfully.\n")
    except Exception as e:
        update_stage(manifest_path, "stage6_post_trim_qc", {"status": "failed", "error": str(e)})
        finalize_manifest(manifest_path, status="failed")
        print("\nStage 6 FAILED.")
        print("Error:", str(e))
        sys.exit(1)


    finalize_manifest(manifest_path, status="stage_6_success")
    print("Pipeline finished (Stage 0–6).\n")
if __name__ == "__main__":
    main()