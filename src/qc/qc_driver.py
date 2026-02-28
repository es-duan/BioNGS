# src/qc/qc_driver.py
import argparse
import json
import subprocess
import sys
from pathlib import Path
from datetime import datetime


def repo_root() -> Path:
    # this file: <root>/src/qc/qc_driver.py
    return Path(__file__).resolve().parents[2]


def run(cmd: list[str]) -> None:
    print(f"[qc_driver] RUN: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def write_manifest_csv(out_csv: Path, rows: list[dict]) -> None:
    # rows: [{"label": "...", "r1": "...", "r2": "..."}]
    import csv

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["label", "r1", "r2"])
        w.writeheader()
        for r in rows:
            w.writerow(r)


def build_raw_manifest(exp: str, gw_name: str) -> list[dict]:
    root = repo_root()
    input_dir = root / "input_data" / exp

    # recursive find GW_name_R1/R2
    r1 = sorted(input_dir.rglob(f"{gw_name}_R1*.fastq")) + sorted(input_dir.rglob(f"{gw_name}_R1*.fq"))
    r2 = sorted(input_dir.rglob(f"{gw_name}_R2*.fastq")) + sorted(input_dir.rglob(f"{gw_name}_R2*.fq"))
    if not r1 or not r2:
        raise FileNotFoundError(f"Cannot find raw fastq for GW_name={gw_name} under {input_dir}")

    return [{"label": "RAW", "r1": str(r1[0]), "r2": str(r2[0])}]


def build_demux_manifest(demux_dir: Path) -> list[dict]:
    # demux_dir contains P1/P1_R1.fastq, P1_R2.fastq, ...
    rows = []
    for pop_dir in sorted([p for p in demux_dir.iterdir() if p.is_dir() and p.name.startswith("P")]):
        r1 = pop_dir / f"{pop_dir.name}_R1.fastq"
        r2 = pop_dir / f"{pop_dir.name}_R2.fastq"
        if r1.exists() and r2.exists():
            rows.append({"label": pop_dir.name, "r1": str(r1), "r2": str(r2)})
    if not rows:
        raise FileNotFoundError(f"No demux population fastq found under: {demux_dir}")
    return rows


def run_fancy_qc(stage_name: str, manifest_csv: Path, out_dir: Path) -> None:
    root = repo_root()
    ENTRY_QC = root / "src" / "qc" / "entry_qc.py"
    if not ENTRY_QC.exists():
        raise FileNotFoundError(f"entry_qc.py not found: {ENTRY_QC}")

    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable,
        str(ENTRY_QC),
        "--manifest",
        str(manifest_csv),
        "--outdir",
        str(out_dir),
        "--title",
        stage_name,
    ]
    run(cmd)


def run_basic_overview_after_demux(exp: str, gw_name: str, demux_dir: Path, overview_run_dir: Path, min_len: int) -> None:
    """
    Run the 'idiot-proof' overview plots (read distribution + read length histogram)
    and save into results/<exp>/qc_overview/run_<min_len>/ (not in demultiplexing).
    """
    root = repo_root()
    CHECK = root / "src" / "check_index_quality.py"
    if not CHECK.exists():
        raise FileNotFoundError(f"check_index_quality.py not found: {CHECK}")

    overview_run_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable,
        str(CHECK),
        exp,
        "--demux_dir",
        str(demux_dir),
        "--outdir",
        str(overview_run_dir),
        "--min_len",
        str(min_len),
        "--gw_name",
        gw_name,
    ]
    run(cmd)


def warn_and_maybe_rerun(exp: str, gw_name: str, base_run_min_len: int, threshold: float, interactive: bool) -> int | None:
    """Return new min_len if user wants rerun, else None."""
    root = repo_root()
    metrics_path = root / "results" / exp / "qc_overview" / f"run_{base_run_min_len}" / "metrics.json"
    if not metrics_path.exists():
        print(f"[qc_driver] NOTE: metrics.json not found at {metrics_path}; skip warning check.")
        return None

    metrics = json.loads(metrics_path.read_text())
    discard_rate = metrics.get("overall", {}).get("discard_rate", 0.0)

    if discard_rate <= threshold:
        print(f"[qc_driver] discard_rate={discard_rate:.4f} <= {threshold:.2f}; no warning.")
        return None

    print("\n" + "=" * 70)
    print("[qc_driver] ⚠️  HIGH DISCARD WARNING")
    print(f"  Experiment: {exp} | GW_name: {gw_name} | min_len: {base_run_min_len}")
    print(f"  discard_rate(short/total) = {discard_rate:.2%}  (threshold = {threshold:.0%})")
    print("  This may indicate min_len is overly restrictive or data characteristics are unexpected.")
    print("\n  Suggested re-run command example:")
    print(f"    python -m src.qc.qc_driver {exp} --gw_name {gw_name} --min_len 120 --interactive 0")
    print("=" * 70 + "\n")

    if not interactive:
        return None

    ans = input("Discard rate is high. Re-run with a different min_len? [y/N]: ").strip().lower()
    if ans not in ("y", "yes"):
        return None

    while True:
        s = input("Enter new min_len (integer, e.g., 120): ").strip()
        try:
            v = int(s)
            if v <= 0:
                raise ValueError
            return v
        except Exception:
            print("Invalid min_len. Please input a positive integer.")

def run_umi_steps(exp: str, demux_dir: Path, run_tag: str) -> None:
    root = repo_root()

    umi_script = root / "src" / "demultiplex_UMI.py"
    umi_qc_script = root / "src" / "check_UMI_quality.py"

    if not umi_script.exists():
        raise FileNotFoundError(f"demultiplex_UMI.py not found: {umi_script}")
    if not umi_qc_script.exists():
        raise FileNotFoundError(f"check_UMI_quality.py not found: {umi_qc_script}")

    # 统一 run 层级
    umi_outdir = root / "results" / exp / "umi" / run_tag
    umi_qc_outdir = root / "results" / exp / "UMI_quality" / run_tag

    umi_outdir.mkdir(parents=True, exist_ok=True)
    umi_qc_outdir.mkdir(parents=True, exist_ok=True)

    print(f"[qc_driver] Running UMI dict build for {run_tag}")

    # 必须传 umi_outdir
    run([
        sys.executable,
        str(umi_script),
        exp,
        "--demux_dir", str(demux_dir),
        "--umi_outdir", str(umi_outdir),
    ])

    print(f"[qc_driver] Running UMI QC for {run_tag}")

    run([
        sys.executable,
        str(umi_qc_script),
        exp,
        "--demux_dir", str(umi_outdir),  # 扫描刚才生成的 dict
        "--outdir", str(umi_qc_outdir),
    ])

    print(f"[qc_driver] UMI dicts: {umi_outdir}")
    print(f"[qc_driver] UMI quality outputs: {umi_qc_outdir}")

def run_one(exp: str, gw_name: str, min_len: int) -> None:
    root = repo_root()
    run_tag = f"run_{min_len}"

    # -------------------------
    # (1) demultiplexing: ONLY sequencing outputs
    # -------------------------
    demux_dir = root / "results" / exp / "demultiplexing" / run_tag
    demux_dir.mkdir(parents=True, exist_ok=True)

    # -------------------------
    # (2) qc_overview: idiot-proof summary + metrics (per run)
    # -------------------------
    overview_run_dir = root / "results" / exp / "qc_overview" / run_tag
    overview_run_dir.mkdir(parents=True, exist_ok=True)
    metrics_json = overview_run_dir / "metrics.json"

    # -------------------------
    # (3) qc_details: fancy reports (raw shared + after-demux per run)
    # -------------------------
    qc_details_root = root / "results" / exp / "qc_details"
    qc_raw_dir = qc_details_root / "00_raw"
    qc_after_demux_dir = qc_details_root / "02_after_demux" / run_tag

    # -------------------------
    # (4) manifests (keep separate from demultiplexing)
    # -------------------------
    manifest_dir = root / "results" / exp / "manifests" / run_tag
    manifest_dir.mkdir(parents=True, exist_ok=True)
    raw_manifest_csv = manifest_dir / "raw_manifest.csv"
    demux_manifest_csv = manifest_dir / "demux_manifest.csv"

    # ---- (QC RAW) fancy report (does NOT depend on min_len, but fine to regenerate)
    raw_rows = build_raw_manifest(exp, gw_name)
    write_manifest_csv(raw_manifest_csv, raw_rows)
    run_fancy_qc(stage_name=f"Fancy QC (RAW) – {exp}", manifest_csv=raw_manifest_csv, out_dir=qc_raw_dir)

    # ---- Script1: create demux folder structure under demultiplexing/run_<min_len>
    run([sys.executable, str(root / "src" / "demultiplex_folders.py"), exp, "--output_dir", str(demux_dir)])

    # ---- Script2: demux with min_len; write metrics.json into qc_overview (NOT demultiplexing)
    run(
        [
            sys.executable,
            str(root / "src" / "demultiplex_index.py"),
            exp,
            "--min_len",
            str(min_len),
            "--output_dir",
            str(demux_dir),
            "--metrics_json",
            str(metrics_json),
        ]
    )

    # ---- (QC OVERVIEW AFTER DEMUX) idiot-proof plots into qc_overview/run_<min_len>
    run_basic_overview_after_demux(exp=exp, gw_name=gw_name, demux_dir=demux_dir, overview_run_dir=overview_run_dir, min_len=min_len)

    # ---- (QC DETAILS AFTER DEMUX) fancy report for each P*
    demux_rows = build_demux_manifest(demux_dir)
    write_manifest_csv(demux_manifest_csv, demux_rows)
    run_fancy_qc(
        stage_name=f"Fancy QC (After Demux, min_len={min_len}) – {exp}",
        manifest_csv=demux_manifest_csv,
        out_dir=qc_after_demux_dir,
    )

    print(f"\n[qc_driver] DONE run={run_tag}.")
    print(f"  Demux outputs (data only): {demux_dir}")
    print(f"  QC overview: {overview_run_dir}")
    print(f"  QC details:  {qc_after_demux_dir}\n")


def main():
    ap = argparse.ArgumentParser(description="QC driver with run_150 baseline + warning + interactive rerun")
    ap.add_argument("experiment", help="Experiment name, e.g. example")
    ap.add_argument("--gw_name", required=True, help="GW_name, e.g. P22R1")
    ap.add_argument("--min_len", type=int, default=150, help="Baseline min_len (default 150)")
    ap.add_argument("--warn_discard_rate", type=float, default=0.30, help="Warn if discard_rate > this (default 0.30)")
    ap.add_argument("--interactive", type=int, default=1, help="1=ask rerun prompt, 0=non-interactive")

    args = ap.parse_args()

    # 1) baseline run
    run_one(args.experiment, args.gw_name, args.min_len)

    # 2) warning + optional rerun
    new_min_len = warn_and_maybe_rerun(
        exp=args.experiment,
        gw_name=args.gw_name,
        base_run_min_len=args.min_len,
        threshold=args.warn_discard_rate,
        interactive=bool(args.interactive),
    )

    # 3) decide final run tag
    final_min_len = new_min_len if new_min_len is not None else args.min_len

    # if user chose rerun, generate it first
    if new_min_len is not None:
        run_one(args.experiment, args.gw_name, new_min_len)

    # 4) run UMI steps ONCE on the final run
    root = repo_root()
    final_run_tag = f"run_{final_min_len}"
    final_demux_dir = root / "results" / args.experiment / "demultiplexing" / final_run_tag
    run_umi_steps(args.experiment, final_demux_dir, final_run_tag)


if __name__ == "__main__":
    main()