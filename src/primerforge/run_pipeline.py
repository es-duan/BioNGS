# src/run_pipeline.py
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from primerforge.utils.paths import results_dir
from primerforge.pipeline.run_manifest import (
    init_run_manifest,
    update_stage,
    finalize_manifest,
    write_json,
)
from primerforge.pipeline.stage0_validation import stage0_validate
from primerforge.pipeline.stage0_validation import print_stage0_summary
from primerforge.pipeline.stage1_raw_qc import stage1_raw_qc
from primerforge.pipeline.stage2_demux import stage2_demux
from primerforge.pipeline.stage2_5_demux_qc import stage2_5_demux_qc
from primerforge.pipeline.stage3_umi_extract import stage3_umi_extract
from primerforge.pipeline.stage3_5_umi_quality import stage3_5_umi_quality
from primerforge.pipeline.stage4_clean_qc import stage4_clean_qc
from primerforge.pipeline.stage5_preprocess import stage5_preprocess
from primerforge.pipeline.stage6_post_trim_qc import stage6_post_trim_qc
from primerforge.pipeline.stage7_family_strategy import stage7_family_strategy
from primerforge.pipeline.stage8_alignment import stage8_alignment
from primerforge.pipeline.stage8_5_alignment_qc import stage8_5_alignment_qc

from primerforge.utils.resume import ResumeManager, file_fingerprint
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="BioNGS Pipeline")
    p.add_argument("--exp", required=True, help="Experiment name")
    p.add_argument("--run_tag", required=True, help="Run tag (e.g. run1)")
    p.add_argument("--r1", required=True, help="Path to R1 FASTQ")
    p.add_argument("--r2", required=True, help="Path to R2 FASTQ")
    p.add_argument("--multiplex_csv", required=True, help="Multiplex CSV file")
    p.add_argument("--umi_primers_csv", required=True, help="UMI primers CSV (columns: f,r)")
    p.add_argument("--no_resume", action="store_true", help="Disable resume; always run all stages")
    p.add_argument("--rerun_from", default=None, help="Force rerun from this stage name (inclusive)")
    p.add_argument("--list_stages", action="store_true", help="Print stage names and exit")
    p.add_argument(
        "--reference_index",
        required=False,
        default="input_data/reference_docs/rpoB_bowtie_index/rpoBsequence_index",
        help="Bowtie2 index prefix"
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    exp = args.exp
    run_tag = args.run_tag
    r1 = Path(args.r1).resolve()
    r2 = Path(args.r2).resolve()
    multiplex_csv = Path(args.multiplex_csv).resolve()
    umi_primers_csv = Path(args.umi_primers_csv).resolve()
    reference_index = Path(args.reference_index).resolve()

    exp_dir = results_dir(exp)
    metrics_dir = exp_dir / "metrics"
    metrics_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = exp_dir / "run_manifest.json"
    resume_enabled = (not args.no_resume)
    resume = ResumeManager(exp_dir=exp_dir, run_tag=run_tag)
    base_sig = {
        "exp": exp,
        "run_tag": run_tag,
        "r1": file_fingerprint(r1),
        "r2": file_fingerprint(r2),
        "multiplex_csv": file_fingerprint(multiplex_csv),
        "umi_primers_csv": file_fingerprint(umi_primers_csv),
    }

    STAGE_ORDER = [
        "stage0_validation",
        "stage1_raw_qc",
        "stage2_demux",
        "stage2_5_demux_qc",
        "stage3_umi_extract",
        "stage3_5_umi_quality",
        "stage4_clean_qc",
        "stage5_preprocess",
        "stage6_post_trim_qc",
        "stage7_family_strategy",
        "stage8_alignment",
        "stage8_5_alignment_qc",
    ]

    if args.list_stages:
        print("\n".join(STAGE_ORDER))
        return

    # rerun-from: delete resume records for that stage and after
    if args.rerun_from:
        if args.rerun_from not in STAGE_ORDER:
            raise SystemExit(f"--rerun_from must be one of: {', '.join(STAGE_ORDER)}")
        start = STAGE_ORDER.index(args.rerun_from)
        for s in STAGE_ORDER[start:]:
            resume.clear_stage(s)
        print(f"[RESUME] cleared from {args.rerun_from} (inclusive)")
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

    STAGE = "stage0_validation"
    sig = {
        "stage": STAGE,
        "r1": file_fingerprint(r1),
        "r2": file_fingerprint(r2),
        "multiplex_csv": file_fingerprint(multiplex_csv),
    }

    try:

        if resume_enabled and resume.is_done(STAGE, sig):
            print("[RESUME] stage0_validation already done — skipping\n")
            stage0_result = {"status": "skipped_resume"}

        else:
            stage0_result = stage0_validate(r1=r1, r2=r2, multiplex_csv=multiplex_csv)

            print_stage0_summary(stage0_result, as_json=False)

            write_json(metrics_dir / "metrics_stage0_validation.json", stage0_result)

            if resume_enabled:
                resume.mark_done(STAGE, sig)

        update_stage(manifest_path, STAGE, {"status": "success", "result": stage0_result})

        print("Stage 0 completed successfully.\n")

    except Exception as e:
        update_stage(manifest_path, STAGE, {"status": "failed", "error": str(e)})
        finalize_manifest(manifest_path, status="failed")
        print("\nStage 0 FAILED.")
        print("Error:", str(e))
        sys.exit(1)

    # ---------- Stage 1 ----------
    print("Running Stage 1: Raw QC (FastQC details + custom overview)...")

    STAGE = "stage1_raw_qc"
    sig = {**base_sig, "stage": STAGE}

    try:
        if resume_enabled and resume.is_done(STAGE, sig):
            print("[RESUME] stage1_raw_qc already done — skipping\n")
            stage1_result = {"status": "skipped_resume"}
        else:
            stage1_result = stage1_raw_qc(
                exp=exp,
                run_tag=run_tag,
                r1=r1,
                r2=r2,
                results_exp_dir=exp_dir
            )

            if stage1_result.get("status") == "aborted":
                update_stage(manifest_path, STAGE, stage1_result)
                finalize_manifest(manifest_path, status="aborted_stage1")
                print("\nPipeline aborted by user at Stage 1.")
                sys.exit(2)

            # ✅ only mark done when success
            if resume_enabled:
                resume.mark_done(STAGE, sig)

        update_stage(manifest_path, STAGE, stage1_result)
        print("\nStage 1 completed successfully.\n")

    except Exception as e:
        update_stage(manifest_path, STAGE, {"status": "failed", "error": str(e)})
        finalize_manifest(manifest_path, status="failed")
        print("\nStage 1 FAILED.")
        print("Error:", str(e))
        sys.exit(1)
    # ---------- Stage 2 ----------
    print("Running Stage 2: Demultiplex (inline index + strip index)...")

    STAGE = "stage2_demux"
    sig = {**base_sig, "stage": STAGE}

    try:
        if resume_enabled and resume.is_done(STAGE, sig):
            print("[RESUME] stage2_demux already done — skipping\n")
            stage2_result = {"status": "skipped_resume"}
        else:
            stage2_result = stage2_demux(
                exp=exp,
                run_tag=run_tag,
                r1=r1,
                r2=r2,
                multiplex_csv=multiplex_csv,
                results_exp_dir=exp_dir,
            )

            if stage2_result.get("status") == "aborted":
                update_stage(manifest_path, STAGE, stage2_result)
                finalize_manifest(manifest_path, status="aborted_stage2")
                print("\nPipeline aborted by user at Stage 2.")
                sys.exit(2)

            #  only mark done when success
            if resume_enabled:
                resume.mark_done(STAGE, sig)

        update_stage(manifest_path, STAGE, stage2_result)
        print("\nStage 2 completed successfully.\n")

    except Exception as e:
        update_stage(manifest_path, STAGE, {"status": "failed", "error": str(e)})
        finalize_manifest(manifest_path, status="failed")
        print("\nStage 2 FAILED.")
        print("Error:", str(e))
        sys.exit(1)

    # ---------- Stage 2.5 ----------
    print("Running Stage 2.5: Demux QC (per-pop FastQC details + custom overview)...")

    STAGE = "stage2_5_demux_qc"
    sig = {**base_sig, "stage": STAGE}

    try:
        if resume_enabled and resume.is_done(STAGE, sig):
            print("[RESUME] stage2_5_demux_qc already done — skipping\n")
            stage2_5_result = {"status": "skipped_resume"}
        else:
            stage2_5_result = stage2_5_demux_qc(
                exp=exp,
                run_tag=run_tag,
                results_exp_dir=exp_dir,
            )

            if stage2_5_result.get("status") == "aborted":
                update_stage(manifest_path, STAGE, stage2_5_result)
                finalize_manifest(manifest_path, status="aborted_stage2_5")
                print("\nPipeline aborted by user at Stage 2.5.")
                sys.exit(2)

            #  only mark done when success
            if resume_enabled:
                resume.mark_done(STAGE, sig)

        update_stage(manifest_path, STAGE, stage2_5_result)
        print("\nStage 2.5 completed successfully.\n")

    except Exception as e:
        update_stage(manifest_path, STAGE, {"status": "failed", "error": str(e)})
        finalize_manifest(manifest_path, status="failed")
        print("\nStage 2.5 FAILED.")
        print("Error:", str(e))
        sys.exit(1)

    # ---------- Stage 3 ----------
    print("Running Stage 3: UMI extraction (trim primers + UMI, tag header)...")

    STAGE = "stage3_umi_extract"
    sig = {**base_sig, "stage": STAGE}

    try:
        if resume_enabled and resume.is_done(STAGE, sig):
            print("[RESUME] stage3_umi_extract already done — skipping\n")
            # 关键：从 resume 文件里恢复 stage3_result，供 stage3.5 使用
            rec = resume.load(STAGE)
            stage3_result = (rec.get("extra") or {}).get("stage3_result") or {"status": "skipped_resume"}
            stage3_result.setdefault("status", "skipped_resume")
        else:
            stage3_result = stage3_umi_extract(
                exp=exp,
                run_tag=run_tag,
                results_exp_dir=exp_dir,
                umi_primers_csv=umi_primers_csv,
            )

            update_stage(manifest_path, STAGE, stage3_result)

            if stage3_result.get("status") == "aborted":
                finalize_manifest(manifest_path, status="aborted_stage3")
                print("\nPipeline aborted by user at Stage 3.")
                sys.exit(2)

            # 写 resume 标记，并把 stage3_result 存进去，给 skip 时恢复用
            if resume_enabled:
                resume.mark_done(STAGE, sig, extra={"stage3_result": stage3_result})

        # 注意：如果是 skip，我们还没 update_stage，这里统一 update 一次
        update_stage(manifest_path, STAGE, stage3_result)
        print("\nStage 3 completed successfully.\n")

    except Exception as e:
        update_stage(manifest_path, STAGE, {"status": "failed", "error": str(e)})
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

    STAGE = "stage4_clean_qc"
    sig = {**base_sig, "stage": STAGE}

    try:
        if resume_enabled and resume.is_done(STAGE, sig):
            print("[RESUME] stage4_clean_qc already done — skipping\n")
            stage4_result = {"status": "skipped_resume"}
        else:
            stage4_result = stage4_clean_qc(
                exp=exp,
                run_tag=run_tag,
                results_exp_dir=exp_dir
            )

            #  写 metrics（只有真的跑时才写）
            write_json(metrics_dir / "metrics_stage4_clean_qc.json", stage4_result)

            #  only mark done when success
            if resume_enabled:
                resume.mark_done(STAGE, sig)

        update_stage(manifest_path, STAGE, stage4_result)
        print("\nStage 4 completed successfully.\n")

    except Exception as e:
        update_stage(manifest_path, STAGE, {"status": "failed", "error": str(e)})
        finalize_manifest(manifest_path, status="failed")
        print("\nStage 4 FAILED.")
        print("Error:", str(e))
        sys.exit(1)

    # ---------- Stage 5 ----------
    print("Running Stage 5: Preprocessing (fastp tail-trim + interactive min_len preview + length filtering)...")

    STAGE = "stage5_preprocess"
    sig = {**base_sig, "stage": STAGE}

    try:
        if resume_enabled and resume.is_done(STAGE, sig):
            print("[RESUME] stage5_preprocess already done — skipping\n")
            rec = resume.load(STAGE)
            stage5_result = (rec.get("extra") or {}).get("stage5_result") or {"status": "skipped_resume"}
            stage5_result.setdefault("status", "skipped_resume")
        else:
            stage5_result = stage5_preprocess(
                exp=exp,
                run_tag=run_tag,
                results_exp_dir=exp_dir,
                metrics_dir=metrics_dir,
            )

            update_stage(manifest_path, STAGE, stage5_result)

            if stage5_result.get("status") == "aborted":
                finalize_manifest(manifest_path, status="aborted_stage5")
                print("\nPipeline aborted by user at Stage 5.")
                sys.exit(2)

            if resume_enabled:
                resume.mark_done(STAGE, sig, extra={"stage5_result": stage5_result})

        update_stage(manifest_path, STAGE, stage5_result)
        print("\nStage 5 completed successfully.\n")

    except Exception as e:
        update_stage(manifest_path, STAGE, {"status": "failed", "error": str(e)})
        finalize_manifest(manifest_path, status="failed")
        print("\nStage 5 FAILED.")
        print("Error:", str(e))
        sys.exit(1)

    # ---------- Stage 6 ----------
    print("Running Stage 6: Post-preprocess QC (FastQC details + overview txt + freq plot)...")

    STAGE = "stage6_post_trim_qc"

    try:
        trimmed_dir = Path(stage5_result["trimmed_dir"])

        sig = {
            **base_sig,
            "stage": STAGE,
            "trimmed_dir": file_fingerprint(trimmed_dir),
        }

        if resume_enabled and resume.is_done(STAGE, sig):
            print("[RESUME] stage6_post_trim_qc already done — skipping\n")
            stage6_result = {"status": "skipped_resume"}
        else:
            stage6_result = stage6_post_trim_qc(
                exp=exp,
                run_tag=run_tag,
                trimmed_dir=trimmed_dir,
                results_exp_dir=exp_dir,
            )

            if resume_enabled:
                resume.mark_done(STAGE, sig)

        update_stage(manifest_path, STAGE, stage6_result)
        print("\nStage 6 completed successfully.\n")

    except Exception as e:
        update_stage(manifest_path, STAGE, {"status": "failed", "error": str(e)})
        finalize_manifest(manifest_path, status="failed")
        print("\nStage 6 FAILED.")
        print("Error:", str(e))
        sys.exit(1)

    # ---------- Stage 7 ----------
    print("Running Stage 7: Family strategy selection (reporting / collapse / consensus / family-resolved)...")

    STAGE = "stage7_family_strategy"
    sig = {**base_sig, "stage": STAGE}

    try:
        if resume_enabled and resume.is_done(STAGE, sig):
            print("[RESUME] stage7_family_strategy already done — skipping\n")
            stage7_result = {"status": "skipped_resume"}
        else:
            stage7_result = stage7_family_strategy(
                exp=exp,
                run_tag=run_tag,
                results_exp_dir=exp_dir,
                metrics_dir=metrics_dir,
                trimmed_version=None,
            )
            if resume_enabled:
                resume.mark_done(STAGE, sig)

        update_stage(manifest_path, STAGE, stage7_result)
        print("\nStage 7 completed successfully.\n")

    except Exception as e:
        update_stage(manifest_path, STAGE, {"status": "failed", "error": str(e)})
        finalize_manifest(manifest_path, status="failed")
        print("\nStage 7 FAILED.")
        print("Error:", str(e))
        sys.exit(1)

    # ---------- Stage 8 ----------
    print("Running Stage 8: Alignment (Bowtie2 + samtools)...")

    STAGE = "stage8_alignment"
    sig = {**base_sig, "stage": STAGE, "reference_index": str(reference_index)}

    try:
        if resume_enabled and resume.is_done(STAGE, sig):
            print("[RESUME] stage8_alignment already done — skipping\n")
            stage8_result = {"status": "skipped_resume"}
        else:
            stage8_result = stage8_alignment(
                exp=exp,
                run_tag=run_tag,
                results_exp_dir=exp_dir,
                reference_index=reference_index,
            )
            if resume_enabled:
                resume.mark_done(STAGE, sig)

        update_stage(manifest_path, STAGE, stage8_result)
        print("\nStage 8 completed successfully.\n")

    except Exception as e:
        update_stage(manifest_path, STAGE, {"status": "failed", "error": str(e)})
        finalize_manifest(manifest_path, status="failed")
        print("\nStage 8 FAILED.")
        print("Error:", str(e))
        sys.exit(1)

    # ---------- Stage 8.5 ----------
    print("Running Stage 8.5: Alignment QC (flagstat / coverage / depth)...")

    STAGE = "stage8_5_alignment_qc"
    sig = {**base_sig, "stage": STAGE}

    try:
        if resume_enabled and resume.is_done(STAGE, sig):
            print("[RESUME] stage8_5_alignment_qc already done — skipping\n")
            stage8_5_result = {"status": "skipped_resume"}
        else:
            stage8_5_result = stage8_5_alignment_qc(
                exp=exp,
                run_tag=run_tag,
                results_exp_dir=exp_dir,
            )
            if resume_enabled:
                resume.mark_done(STAGE, sig)

        update_stage(manifest_path, STAGE, stage8_5_result)
        print("\nStage 8.5 completed successfully.\n")

    except Exception as e:
        update_stage(manifest_path, STAGE, {"status": "failed", "error": str(e)})
        finalize_manifest(manifest_path, status="failed")
        print("\nStage 8.5 FAILED.")
        print("Error:", str(e))
        sys.exit(1)

    finalize_manifest(manifest_path, status="stage_8_5_success")
    print("Pipeline finished (Stage 0–8.5).\n")

if __name__ == "__main__":
    main()