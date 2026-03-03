from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Dict, Any, List, Tuple
import gzip

from src.utils.decision import print_level, ask_stop_if_abnormal
from src.utils.fastq_tools import (
    count_fastq_reads,
    write_fastq_sample_pairs,
    preview_pair_short_rate_from_trimmed,
    filter_pairs_by_min_len,
)
from src.qc.overview_raw import compute_raw_overview, write_overview_report


def _run_cmd(cmd: List[str]) -> None:
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed ({p.returncode}): {' '.join(cmd)}\n{p.stdout}")


def _write_json(path: Path, obj: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)


def _find_pop_pairs_from_umi_extracted(results_exp_dir: Path, run_tag: str) -> List[Tuple[str, Path, Path]]:
    """
    输入：results/<exp>/umi_extracted/<run_tag>/<pop>/*_R1.fastq(.gz), *_R2.fastq(.gz)
    """
    root = results_exp_dir / "umi_extracted" / run_tag
    if not root.exists():
        raise FileNotFoundError(f"umi_extracted dir not found: {root}")

    pairs: List[Tuple[str, Path, Path]] = []
    for pop_dir in sorted([p for p in root.iterdir() if p.is_dir()]):
        pop = pop_dir.name
        r1s = sorted(list(pop_dir.glob("*_R1.fastq")) + list(pop_dir.glob("*_R1.fastq.gz")))
        r2s = sorted(list(pop_dir.glob("*_R2.fastq")) + list(pop_dir.glob("*_R2.fastq.gz")))
        if not r1s or not r2s:
            continue
        pairs.append((pop, r1s[0], r2s[0]))

    if not pairs:
        raise FileNotFoundError(f"No pop FASTQ pairs found under: {root}")
    return pairs


def _parse_fastp_json(fastp_json: Path) -> Dict[str, Any]:
    d = json.loads(fastp_json.read_text(encoding="utf-8"))
    # 兼容 fastp 常见结构
    summ = d.get("summary", {})
    bf = summ.get("before_filtering", {})
    af = summ.get("after_filtering", {})
    # fastp 提供 q30_bases/total_bases
    def q30_rate(x: Dict[str, Any]) -> float:
        tb = float(x.get("total_bases", 0) or 0)
        q30 = float(x.get("q30_bases", 0) or 0)
        return (q30 / tb) if tb else 0.0

    total_reads_after = float(af.get("total_reads", 0) or 0)
    total_bases_after = float(af.get("total_bases", 0) or 0)
    mean_len = (total_bases_after / total_reads_after) if total_reads_after else 0.0

    # trimmed bases：不同版本字段可能在 filtering_result 或者 read1/read2
    fr = d.get("filtering_result", {})
    trimmed_bases = fr.get("trimmed_bases", None)

    return {
        "q30_before": q30_rate(bf),
        "q30_after": q30_rate(af),
        "total_reads_before": bf.get("total_reads", None),
        "total_reads_after": af.get("total_reads", None),
        "total_bases_before": bf.get("total_bases", None),
        "total_bases_after": af.get("total_bases", None),
        "mean_len_after": mean_len,
        "trimmed_bases": trimmed_bases,
    }


def _fastp_trim_only(
    r1_in: Path,
    r2_in: Path,
    r1_out: Path,
    r2_out: Path,
    html_out: Path,
    json_out: Path,
    q: int,
) -> None:
    """
    只做 3' trimming（trim-only），不做长度过滤。
    fastp 用 cut_tail + window=1 近似 “从3'往前剪直到 >=Q”
    """
    r1_out.parent.mkdir(parents=True, exist_ok=True)
    r2_out.parent.mkdir(parents=True, exist_ok=True)
    html_out.parent.mkdir(parents=True, exist_ok=True)
    json_out.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "fastp",
        "-i", str(r1_in),
        "-I", str(r2_in),
        "-o", str(r1_out),
        "-O", str(r2_out),
        "--cut_tail",
        "--cut_tail_window_size", "1",
        "--cut_tail_mean_quality", str(q),
        "--disable_length_filtering",
        "--json", str(json_out),
        "--html", str(html_out),
    ]
    _run_cmd(cmd)


def stage5_preprocess(
    exp: str,
    run_tag: str,
    results_exp_dir: Path,
    metrics_dir: Path,
) -> Dict[str, Any]:
    """
    Stage 5:
      5.1 candidate Q evaluation (sampled) -> user choose q_threshold -> generate trim-only FASTQ (versioned v1/v2...)
      5.2 qc_overview on trim-only + user choose min_len -> preview short_rate/kept_rate -> warning logic + (ABNORMAL -> ask stop)
          if ABNORMAL and user chooses NOT stop -> go back to 5.1 (new version)
      5.3 execute length filtering -> final trimmed FASTQ.gz
    """
    pop_pairs = _find_pop_pairs_from_umi_extracted(results_exp_dir, run_tag)

    # 输出根目录（版本化）
    trim_only_root = results_exp_dir / "trim_only" / run_tag
    trimmed_root = results_exp_dir / "trimmed" / run_tag
    trim_reports_root = results_exp_dir / "trim_reports" / run_tag
    qc_overview_root = results_exp_dir / "qc_overview" / "03_trim_only" / run_tag

    # 版本号：v1, v2...（不覆盖）
    version = 1
    while (trim_only_root / f"v{version}").exists():
        version += 1

    candidates = [30, 29, 28, 25, 20]

    while True:
        vtag = f"v{version}"
        v_trim_only = trim_only_root / vtag
        v_reports = trim_reports_root / vtag
        v_qc_over = qc_overview_root / vtag
        v_eval_root = results_exp_dir / "trim_eval" / run_tag / vtag

        print(f"\nRunning Stage 5.1: fastp tail-trim ONLY (candidate Q eval) ... [{vtag}]")
        print("Candidate Q evaluation (sampled 30% reads, capped at 200,000 pairs)")

        eval_rows: List[Dict[str, Any]] = []
        per_pop_eval: Dict[str, Any] = {}

        for pop, r1, r2 in pop_pairs:
            # 计算抽样 pairs 数：min(30%, 200000)
            total_pairs = min(count_fastq_reads(r1), count_fastq_reads(r2))
            sample_pairs = min(int(total_pairs * 0.30), 200_000)
            if sample_pairs <= 0:
                sample_pairs = min(total_pairs, 10_000)  # 兜底

            sample_dir = v_eval_root / pop / "sample"
            s_r1 = sample_dir / f"{pop}.sample.R1.fastq"
            s_r2 = sample_dir / f"{pop}.sample.R2.fastq"
            sample_info = write_fastq_sample_pairs(r1, r2, s_r1, s_r2, sample_pairs)

            per_pop_eval[pop] = {
                "total_pairs_est": total_pairs,
                "sample": sample_info,
                "candidates": {},
            }

            # 对每个Q跑 fastp（sample）
            for q in candidates:
                out_dir = v_eval_root / pop / f"Q{q}"
                out_r1 = out_dir / f"{pop}.Q{q}.trim.R1.fastq"
                out_r2 = out_dir / f"{pop}.Q{q}.trim.R2.fastq"
                out_html = out_dir / "fastp.html"
                out_json = out_dir / "fastp.json"

                _fastp_trim_only(s_r1, s_r2, out_r1, out_r2, out_html, out_json, q=q)
                m = _parse_fastp_json(out_json)

                # trimmed_bases% 需要 before/after bases
                tb = float(m.get("total_bases_before") or 0)
                trimmed_bases = float(m.get("trimmed_bases") or 0)
                trimmed_pct = (trimmed_bases / tb) if tb else 0.0

                row = {
                    "pop": pop,
                    "Q": q,
                    "Q30_after": m["q30_after"],
                    "trimmed_bases_pct": trimmed_pct,
                    "mean_len": m["mean_len_after"],
                }
                eval_rows.append(row)

                per_pop_eval[pop]["candidates"][f"Q{q}"] = {
                    "fastp_json": str(out_json),
                    "fastp_html": str(out_html),
                    "metrics": m,
                    "trimmed_bases_pct": trimmed_pct,
                }

        # 终端汇总打印（按你的格式，先合并显示：每个pop单独块，避免混）
        for pop in sorted(per_pop_eval.keys()):
            print(f"\n[Pop {pop}] Candidate Q evaluation (sampled)")
            print("Q   Q30_after   trimmed_bases%   mean_len")
            rows = [r for r in eval_rows if r["pop"] == pop]
            rows.sort(key=lambda x: -x["Q"])
            for r in rows:
                print(f"{r['Q']:<3d} {r['Q30_after']:<10.3f} {r['trimmed_bases_pct']:<14.3%} {r['mean_len']:<7.1f}")

        # 记录 eval json
        eval_json = {
            "stage": "stage5_1_candidate_eval",
            "exp": exp,
            "run_tag": run_tag,
            "version": vtag,
            "candidates": candidates,
            "per_population": per_pop_eval,
        }
        _write_json(metrics_dir / f"metrics_stage5_1_candidate_eval_{vtag}.json", eval_json)

        # 让用户选择 Q（默认 30）
        raw = input("\nEnter q_threshold for trimming (recommended 30): ").strip()
        q_sel = 30 if raw == "" else int(raw)
        print(f"Selected q_threshold = {q_sel} [{vtag}]")

        # 跑全量 trim-only
        stage5_1_outputs: Dict[str, Any] = {"version": vtag, "q_threshold": q_sel, "per_population": {}}
        for pop, r1, r2 in pop_pairs:
            out_dir = v_trim_only / pop
            out_r1 = out_dir / f"{pop}.trim.R1.fastq.gz"
            out_r2 = out_dir / f"{pop}.trim.R2.fastq.gz"
            rep_dir = v_reports / pop
            html_out = rep_dir / "fastp.html"
            json_out = rep_dir / "fastp.json"

            _fastp_trim_only(r1, r2, out_r1, out_r2, html_out, json_out, q=q_sel)
            m = _parse_fastp_json(json_out)

            stage5_1_outputs["per_population"][pop] = {
                "in_r1": str(r1),
                "in_r2": str(r2),
                "trim_only_r1": str(out_r1),
                "trim_only_r2": str(out_r2),
                "fastp_html": str(html_out),
                "fastp_json": str(json_out),
                "q30_before": m["q30_before"],
                "q30_after": m["q30_after"],
            }

        _write_json(metrics_dir / f"metrics_stage5_1_trim_only_{vtag}.json", stage5_1_outputs)

        # 5.2：qc_overview（trim-only 的长度分布+summary），然后让用户输入 min_len 做 preview
        print(f"\nRunning Stage 5.2: qc_overview on trim-only + length filtering preview ... [{vtag}]")

        for pop, _r1, _r2 in pop_pairs:
            t1 = Path(stage5_1_outputs["per_population"][pop]["trim_only_r1"])
            t2 = Path(stage5_1_outputs["per_population"][pop]["trim_only_r2"])

            out_pop = v_qc_over / pop
            # 只做 overview：txt + length histogram
            m1 = compute_raw_overview(t1, short_len=150)
            m2 = compute_raw_overview(t2, short_len=150)
            write_overview_report(m1, out_pop / "R1", title=f"{pop} - Trim-only R1 length distribution ({vtag})")
            write_overview_report(m2, out_pop / "R2", title=f"{pop} - Trim-only R2 length distribution ({vtag})")

        print(f"[Stage 5.2] Trim-only qc_overview written to: {v_qc_over}")

        def _open_fastq(p: Path):
            return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p, "r")

        def compute_pair_short_rates_table(
                r1_fastq: Path,
                r2_fastq: Path,
                cutoffs: List[int] = [180, 170, 160, 150, 140, 130, 120, 110],
        ) -> Dict[int, Dict[str, float]]:
            """
            Pair-level 口径：
            - 对每个 read pair，取 min(len(R1), len(R2)) = Lmin
            - 对于 cutoff=min_len，short if Lmin < min_len
            返回：每个 cutoff 的 short_rate / kept_rate / est_dropped_pairs
            """
            total_pairs = 0
            short_counts = {c: 0 for c in cutoffs}

            with _open_fastq(r1_fastq) as f1, _open_fastq(r2_fastq) as f2:
                while True:
                    h1 = f1.readline()
                    h2 = f2.readline()
                    if not h1 and not h2:
                        break
                    if (not h1) or (not h2):
                        raise RuntimeError("R1/R2 FASTQ records mismatch (EOF not aligned).")

                    s1 = f1.readline().strip()
                    p1 = f1.readline()
                    q1 = f1.readline().strip()

                    s2 = f2.readline().strip()
                    p2 = f2.readline()
                    q2 = f2.readline().strip()

                    # 最小必要 sanity（坏记录直接停机）
                    if not s1 or not q1 or not s2 or not q2 or not p1 or not p2:
                        raise RuntimeError(f"Bad FASTQ record encountered near pair #{total_pairs + 1}.")

                    Lmin = min(len(s1), len(s2))
                    total_pairs += 1

                    for c in cutoffs:
                        if Lmin < c:
                            short_counts[c] += 1

            out: Dict[int, Dict[str, float]] = {}
            for c in cutoffs:
                short = short_counts[c]
                short_rate = short / total_pairs if total_pairs else 0.0
                kept_rate = 1.0 - short_rate
                out[c] = {
                    "short_rate": short_rate,
                    "kept_rate": kept_rate,
                    "est_dropped_pairs": float(short),
                    "total_pairs": float(total_pairs),
                }
            return out

        def print_pair_short_rates_table(
                table: Dict[int, Dict[str, float]],
                title: str = "Length filtering preview table (pair-level, from trim-only FASTQ)",
        ) -> None:
            print("\n" + title)
            print("(Short pair: either mate length < min_len. See qc_overview report for details.)\n")
            print(f"{'min_len':<8} {'short_rate(<min_len)':<22} {'kept_rate':<10} {'est_dropped_pairs':<16}")
            for c in sorted(table.keys(), reverse=True):
                short_rate = table[c]["short_rate"]
                kept_rate = table[c]["kept_rate"]
                dropped = int(table[c]["est_dropped_pairs"])
                print(f"{c:<8} {short_rate:<22.1%} {kept_rate:<10.1%} {dropped:<16,}")
            print("")
        # ✅ 先给用户一个 min_len 候选表格（在输入 min_len 之前）
        cutoffs = [180, 170, 160, 150, 140, 130, 120, 110]

        for pop, _r1, _r2 in pop_pairs:
            t1 = Path(stage5_1_outputs["per_population"][pop]["trim_only_r1"])
            t2 = Path(stage5_1_outputs["per_population"][pop]["trim_only_r2"])

            table = compute_pair_short_rates_table(t1, t2, cutoffs=cutoffs)
            print_pair_short_rates_table(
                table,
                title=f"[Pop {pop}] Length filtering preview (trim-only, pair-level)"
            )

        print("Tip: detailed length distributions are in qc_overview (trim-only). Please inspect reports before choosing min_len.")

        # 用户输入 min_len（做 preview）
        raw2 = input("\nEnter min_len for length filtering (e.g. 150): ").strip()
        min_len = 150 if raw2 == "" else int(raw2)

        preview_all: Dict[str, Any] = {"version": vtag, "min_len": min_len, "per_population": {}}
        # 总体口径：把所有pop加总（pair-level）
        total_pairs = 0
        total_short = 0
        for pop, _r1, _r2 in pop_pairs:
            t1 = Path(stage5_1_outputs["per_population"][pop]["trim_only_r1"])
            t2 = Path(stage5_1_outputs["per_population"][pop]["trim_only_r2"])
            pv = preview_pair_short_rate_from_trimmed(t1, t2, min_len=min_len)
            preview_all["per_population"][pop] = pv
            total_pairs += pv["total_pairs"]
            total_short += pv["short_pairs"]

        short_rate_global = (total_short / total_pairs) if total_pairs else 0.0
        kept_rate_global = 1.0 - short_rate_global

        print("\nLength filtering preview (vX)")
        print("----------------------------------")
        print(f"min_len = {min_len}")
        print(f"short_rate (<{min_len}) = {short_rate_global:.1%}")
        print(f"kept_rate = {kept_rate_global:.1%}")
        print(f"estimated_dropped_reads = {total_short:,}")

        # 规则：>30 warning, >40 strong, >50 abnormal + ask stop
        decision = {"abnormal": False, "stop_prompt": False}
        if short_rate_global > 0.50:
            decision["abnormal"] = True
            decision["stop_prompt"] = True
            msg = (
                f"short_rate {short_rate_global:.1%} (>{50:.0f}%). "
                f"Recommended: re-select q_threshold first because trimming changes length distribution."
            )
            stop = ask_stop_if_abnormal(msg)
            decision["user_input"] = "y" if stop else "n"
            if stop:
                preview_all["decision"] = decision
                _write_json(metrics_dir / f"metrics_stage5_2_len_preview_{vtag}.json", preview_all)
                return {"status": "aborted", "stage": "stage5_2_len_preview", "version": vtag, "min_len": min_len}
            # 用户选 N：回到 5.1，版本号 +1，不覆盖旧版本
            print("Proceeding despite ABNORMAL (user chose not to stop). Returning to Stage 5.1 for new Q selection...")
            preview_all["decision"] = decision
            _write_json(metrics_dir / f"metrics_stage5_2_len_preview_{vtag}.json", preview_all)
            version += 1
            continue
        elif short_rate_global > 0.40:
            print_level(f"short_rate {short_rate_global:.1%} (>40%). Please consider lowering q_threshold or revisiting trimming.", "STRONG WARNING")
        elif short_rate_global > 0.30:
            print_level(f"short_rate {short_rate_global:.1%} (>30%). Please inspect length distribution.", "WARNING")
        else:
            print("Parameter check completed.")
            print("This does not trigger a stop condition (50%).")
            print("Proceeding to length filtering...")

        preview_all["global"] = {
            "total_pairs": total_pairs,
            "short_pairs": total_short,
            "short_rate": short_rate_global,
            "kept_rate": kept_rate_global,
        }
        preview_all["decision"] = decision
        _write_json(metrics_dir / f"metrics_stage5_2_len_preview_{vtag}.json", preview_all)

        # 5.3：真正 length filtering
        print(f"\nRunning Stage 5.3: execute length filtering ... [{vtag}]")
        stage5_3: Dict[str, Any] = {
            "version": vtag,
            "q_threshold": q_sel,
            "min_len": min_len,
            "per_population": {},
            "global": {},
        }

        g_total = 0
        g_kept = 0
        g_drop = 0

        for pop, _r1, _r2 in pop_pairs:
            t1 = Path(stage5_1_outputs["per_population"][pop]["trim_only_r1"])
            t2 = Path(stage5_1_outputs["per_population"][pop]["trim_only_r2"])
            out_dir = trimmed_root / vtag / pop
            out_r1 = out_dir / f"{pop}.trimQ{q_sel}.minlen{min_len}.R1.fastq.gz"
            out_r2 = out_dir / f"{pop}.trimQ{q_sel}.minlen{min_len}.R2.fastq.gz"

            filt = filter_pairs_by_min_len(t1, t2, out_r1, out_r2, min_len=min_len)
            stage5_3["per_population"][pop] = {
                "trim_only_r1": str(t1),
                "trim_only_r2": str(t2),
                "out_r1": str(out_r1),
                "out_r2": str(out_r2),
                "filter_metrics": filt,
            }

            g_total += filt["total_pairs"]
            g_kept += filt["kept_pairs"]
            g_drop += filt["dropped_pairs"]

        stage5_3["global"] = {
            "total_pairs": g_total,
            "kept_pairs": g_kept,
            "dropped_pairs": g_drop,
            "kept_rate": (g_kept / g_total) if g_total else 0.0,
            "short_rate": (g_drop / g_total) if g_total else 0.0,
        }

        _write_json(metrics_dir / f"metrics_stage5_3_len_filter_{vtag}.json", stage5_3)

        print("\nStage 5 completed successfully.")
        print(f"Outputs:")
        print(f"  trim-only: {v_trim_only}")
        print(f"  trimmed:   {trimmed_root / vtag}")
        print(f"  reports:   {v_reports}")
        print(f"  overview:  {v_qc_over}")

        return {
            "status": "success",
            "version": vtag,
            "q_threshold": q_sel,
            "min_len": min_len,
            "trim_only_dir": str(v_trim_only),
            "trimmed_dir": str(trimmed_root / vtag),
            "trim_reports_dir": str(v_reports),
            "qc_overview_dir": str(v_qc_over),
            "global": stage5_3["global"],
        }