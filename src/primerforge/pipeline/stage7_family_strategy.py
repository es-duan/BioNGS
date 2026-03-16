from __future__ import annotations

import csv
import gzip
import json
import random
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

from primerforge.utils.fastq_io import FastqRecord, read_fastq, write_fastq_record, normalize_read_id


UMI_RE = re.compile(r"(?:^|\s|[|;])UMI=([A-Za-z\-]+)")
DEFAULT_RANDOM_SEED = 42
CONSENSUS_MIN_FAMILY_SIZE = 3
CONSENSUS_MIN_FRACTION = 0.70
CONSENSUS_MAX_N_FRACTION = 0.10


@dataclass
class PairedRecord:
    r1: FastqRecord
    r2: FastqRecord
    read_id: str
    pair_length: int


@dataclass
class FamilySummary:
    total_families: int
    size_1: int
    size_2: int
    size_ge_3: int
    mean_family_size: float
    median_family_size: float
    total_pairs: int


class FamilyStrategyError(RuntimeError):
    pass


def _write_json(path: Path, obj: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")
    tmp.replace(path)


def _open_text_maybe_gz(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8")
    return path.open(mode, encoding="utf-8")


def _extract_umi(header: str) -> str:
    m = UMI_RE.search(header)
    if not m:
        raise FamilyStrategyError(f"UMI tag not found in header: {header}")
    return m.group(1)


def _find_trimmed_pairs(results_exp_dir: Path, run_tag: str, version: Optional[str] = None) -> Tuple[str, List[Tuple[str, Path, Path]]]:
    trimmed_root = results_exp_dir / "trimmed" / run_tag
    if not trimmed_root.exists():
        raise FileNotFoundError(f"trimmed dir not found: {trimmed_root}")

    if version is None:
        versions = sorted([p.name for p in trimmed_root.iterdir() if p.is_dir() and p.name.startswith("v")])
        if not versions:
            raise FileNotFoundError(f"No trimmed versions found under: {trimmed_root}")
        version = versions[-1]

    version_dir = trimmed_root / version
    if not version_dir.exists():
        raise FileNotFoundError(f"trimmed version dir not found: {version_dir}")

    pairs: List[Tuple[str, Path, Path]] = []
    for pop_dir in sorted([p for p in version_dir.iterdir() if p.is_dir()]):
        pop = pop_dir.name
        r1s = sorted(pop_dir.glob("*.R1.fastq.gz")) + sorted(pop_dir.glob("*.R1.fastq"))
        r2s = sorted(pop_dir.glob("*.R2.fastq.gz")) + sorted(pop_dir.glob("*.R2.fastq"))
        if not r1s or not r2s:
            continue
        pairs.append((pop, r1s[0], r2s[0]))

    if not pairs:
        raise FileNotFoundError(f"No trimmed FASTQ pairs found under: {version_dir}")
    return version, pairs


def _rebuild_families(r1_path: Path, r2_path: Path) -> Dict[str, List[PairedRecord]]:
    families: Dict[str, List[PairedRecord]] = {}
    it1 = read_fastq(r1_path)
    it2 = read_fastq(r2_path)

    pair_no = 0
    while True:
        try:
            rec1 = next(it1)
            got1 = True
        except StopIteration:
            got1 = False
            rec1 = None
        try:
            rec2 = next(it2)
            got2 = True
        except StopIteration:
            got2 = False
            rec2 = None

        if not got1 and not got2:
            break
        if got1 != got2:
            raise FamilyStrategyError(f"R1/R2 record count mismatch between {r1_path} and {r2_path}")
        assert rec1 is not None and rec2 is not None
        pair_no += 1

        rid1 = normalize_read_id(rec1.header)
        rid2 = normalize_read_id(rec2.header)
        if rid1 != rid2:
            raise FamilyStrategyError(
                f"Paired read id mismatch at pair #{pair_no}: R1={rid1}, R2={rid2}"
            )

        umi1 = _extract_umi(rec1.header)
        umi2 = _extract_umi(rec2.header)
        if umi1 != umi2:
            raise FamilyStrategyError(
                f"UMI mismatch within pair #{pair_no}: R1={umi1}, R2={umi2}"
            )

        families.setdefault(umi1, []).append(
            PairedRecord(
                r1=rec1,
                r2=rec2,
                read_id=rid1,
                pair_length=len(rec1.seq) + len(rec2.seq),
            )
        )

    return families


def _compute_family_summary(families: Dict[str, List[PairedRecord]]) -> FamilySummary:
    sizes = sorted(len(v) for v in families.values())
    total_families = len(sizes)
    total_pairs = sum(sizes)
    if total_families == 0:
        return FamilySummary(0, 0, 0, 0, 0.0, 0.0, 0)

    size_1 = sum(1 for s in sizes if s == 1)
    size_2 = sum(1 for s in sizes if s == 2)
    size_ge_3 = sum(1 for s in sizes if s >= 3)
    mean_size = total_pairs / total_families
    if total_families % 2 == 1:
        median_size = float(sizes[total_families // 2])
    else:
        median_size = (sizes[total_families // 2 - 1] + sizes[total_families // 2]) / 2.0
    return FamilySummary(total_families, size_1, size_2, size_ge_3, mean_size, median_size, total_pairs)


def _print_family_report(pop: str, summary: FamilySummary) -> None:
    print("\n" + "=" * 48)
    print(f"UMI Family Report - Pop {pop}")
    print("=" * 48)
    print(f"Total families:                {summary.total_families}")
    print(f"Total read pairs:              {summary.total_pairs}")
    print(f"Singleton families (size=1):   {summary.size_1}")
    print(f"Doubleton families (size=2):   {summary.size_2}")
    print(f"Families with size>=3:         {summary.size_ge_3}")
    print(f"Mean family size:              {summary.mean_family_size:.3f}")
    print(f"Median family size:            {summary.median_family_size:.3f}")
    print("=" * 48)
    print("Choose one downstream strategy:")
    print("1 - UMI-aware reporting and direct downstream processing")
    print("2 - Representative-family collapse")
    print("3 - UMI family consensus (family size >= 3 only)")
    print("4 - UMI-family-resolved alignment and variant calling")


def _prompt_strategy_choice() -> str:
    valid = {"1": "reporting", "2": "collapse", "3": "consensus", "4": "family_resolved"}
    while True:
        ans = input("Enter selection [1/2/3/4]: ").strip()
        if ans in valid:
            return valid[ans]
        print("Invalid selection. Please enter 1, 2, 3, or 4.")


def _write_family_table(path: Path, families: Dict[str, List[PairedRecord]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["umi", "family_size"])
        for umi in sorted(families.keys()):
            w.writerow([umi, len(families[umi])])


def _make_header(base_header: str, umi: str, family_size: int, mode: str) -> str:
    core = base_header if base_header.startswith("@") else f"@{base_header}"
    return f"{core} | UMI={umi} | family_size={family_size} | mode={mode}"


def _run_reporting_strategy(
    pop: str,
    in_r1: Path,
    in_r2: Path,
    strategy_dir: Path,
    families: Dict[str, List[PairedRecord]],
) -> Dict[str, Any]:
    info = {
        "strategy": "reporting",
        "pop": pop,
        "out_r1": str(in_r1),
        "out_r2": str(in_r2),
        "note": "Direct downstream processing from Stage 5 trimmed FASTQ.",
    }
    _write_json(strategy_dir / "strategy_output.json", info)
    return info


def _select_representative_pair(records: List[PairedRecord], rng: random.Random) -> PairedRecord:
    eligible = [r for r in records if r.r1 and r.r2]
    if not eligible:
        raise FamilyStrategyError("No paired-end reads available for collapse.")
    max_len = max(r.pair_length for r in eligible)
    longest = [r for r in eligible if r.pair_length == max_len]
    return rng.choice(longest)


def _run_collapse_strategy(
    pop: str,
    strategy_dir: Path,
    families: Dict[str, List[PairedRecord]],
    seed: int = DEFAULT_RANDOM_SEED,
) -> Dict[str, Any]:
    rng = random.Random(seed)
    out_r1 = strategy_dir / f"{pop}.collapsed.R1.fastq.gz"
    out_r2 = strategy_dir / f"{pop}.collapsed.R2.fastq.gz"
    decisions_csv = strategy_dir / "collapse_decisions.csv"

    kept = 0
    with _open_text_maybe_gz(out_r1, "wt") as o1, _open_text_maybe_gz(out_r2, "wt") as o2, decisions_csv.open(
        "w", newline="", encoding="utf-8"
    ) as f:
        w = csv.writer(f)
        w.writerow(["umi", "family_size", "selected_read_id", "pair_length", "random_seed"])
        for umi in sorted(families.keys()):
            chosen = _select_representative_pair(families[umi], rng)
            fam_size = len(families[umi])
            r1 = FastqRecord(
                header=_make_header(chosen.r1.header, umi, fam_size, "collapse"),
                seq=chosen.r1.seq,
                plus=chosen.r1.plus,
                qual=chosen.r1.qual,
            )
            r2 = FastqRecord(
                header=_make_header(chosen.r2.header, umi, fam_size, "collapse"),
                seq=chosen.r2.seq,
                plus=chosen.r2.plus,
                qual=chosen.r2.qual,
            )
            write_fastq_record(o1, r1)
            write_fastq_record(o2, r2)
            w.writerow([umi, fam_size, chosen.read_id, chosen.pair_length, seed])
            kept += 1

    info = {
        "strategy": "collapse",
        "pop": pop,
        "out_r1": str(out_r1),
        "out_r2": str(out_r2),
        "collapse_decisions_csv": str(decisions_csv),
        "families_emitted": kept,
        "random_seed": seed,
        "selection_rule": "paired-end only; longest combined read length; tie broken by fixed-seed random choice",
    }
    _write_json(strategy_dir / "strategy_output.json", info)
    return info


def _consensus_base(chars: List[str]) -> Tuple[str, float]:
    if not chars:
        return "N", 0.0
    counts = Counter(chars)
    base, count = counts.most_common(1)[0]
    frac = count / len(chars)
    if frac >= CONSENSUS_MIN_FRACTION:
        return base, frac
    return "N", frac


def _consensus_qual(frac: float, base: str) -> str:
    if base == "N":
        q = 10
    elif frac >= 0.90:
        q = 40
    elif frac >= 0.80:
        q = 30
    else:
        q = 20
    return chr(q + 33)


def _build_consensus_record(records: List[FastqRecord], umi: str, family_size: int, mate_label: str) -> FastqRecord:
    max_len = max(len(r.seq) for r in records)
    seq_out: List[str] = []
    qual_out: List[str] = []

    for i in range(max_len):
        chars = [r.seq[i] for r in records if i < len(r.seq)]
        base, frac = _consensus_base(chars)
        seq_out.append(base)
        qual_out.append(_consensus_qual(frac, base))

    header = f"@{umi}_{mate_label} | UMI={umi} | family_size={family_size} | mode=consensus"
    return FastqRecord(header=header, seq="".join(seq_out), plus="+", qual="".join(qual_out))


def _run_consensus_strategy(
    pop: str,
    strategy_dir: Path,
    families: Dict[str, List[PairedRecord]],
) -> Dict[str, Any]:
    out_r1 = strategy_dir / f"{pop}.consensus.R1.fastq.gz"
    out_r2 = strategy_dir / f"{pop}.consensus.R2.fastq.gz"
    decisions_csv = strategy_dir / "consensus_decisions.csv"
    excluded_csv = strategy_dir / "excluded_families.csv"

    emitted = 0
    excluded = 0
    with _open_text_maybe_gz(out_r1, "wt") as o1, _open_text_maybe_gz(out_r2, "wt") as o2, decisions_csv.open(
        "w", newline="", encoding="utf-8"
    ) as d_f, excluded_csv.open("w", newline="", encoding="utf-8") as e_f:
        d_w = csv.writer(d_f)
        e_w = csv.writer(e_f)
        d_w.writerow(["umi", "family_size", "consensus_r1_len", "consensus_r2_len", "n_fraction_r1", "n_fraction_r2"])
        e_w.writerow(["umi", "family_size", "reason"])

        for umi in sorted(families.keys()):
            fam = families[umi]
            fam_size = len(fam)
            if fam_size < CONSENSUS_MIN_FAMILY_SIZE:
                e_w.writerow([umi, fam_size, "family_size_below_threshold"])
                excluded += 1
                continue

            r1_records = [p.r1 for p in fam]
            r2_records = [p.r2 for p in fam]
            c1 = _build_consensus_record(r1_records, umi, fam_size, "R1")
            c2 = _build_consensus_record(r2_records, umi, fam_size, "R2")

            n_frac_1 = c1.seq.count("N") / len(c1.seq) if c1.seq else 1.0
            n_frac_2 = c2.seq.count("N") / len(c2.seq) if c2.seq else 1.0
            if n_frac_1 > CONSENSUS_MAX_N_FRACTION or n_frac_2 > CONSENSUS_MAX_N_FRACTION:
                e_w.writerow([umi, fam_size, "consensus_too_ambiguous"])
                excluded += 1
                continue

            write_fastq_record(o1, c1)
            write_fastq_record(o2, c2)
            d_w.writerow([umi, fam_size, len(c1.seq), len(c2.seq), f"{n_frac_1:.6f}", f"{n_frac_2:.6f}"])
            emitted += 1

    info = {
        "strategy": "consensus",
        "pop": pop,
        "out_r1": str(out_r1),
        "out_r2": str(out_r2),
        "consensus_decisions_csv": str(decisions_csv),
        "excluded_families_csv": str(excluded_csv),
        "families_emitted": emitted,
        "families_excluded": excluded,
        "min_family_size": CONSENSUS_MIN_FAMILY_SIZE,
        "consensus_min_fraction": CONSENSUS_MIN_FRACTION,
        "consensus_max_n_fraction": CONSENSUS_MAX_N_FRACTION,
    }
    _write_json(strategy_dir / "strategy_output.json", info)
    return info


def _sanitize_umi_for_path(umi: str) -> str:
    return umi.replace("/", "_")


def _run_family_resolved_strategy(
    pop: str,
    strategy_dir: Path,
    families: Dict[str, List[PairedRecord]],
) -> Dict[str, Any]:
    family_fastq_root = strategy_dir / "family_fastq"
    summary_csv = strategy_dir / "family_resolved_summary.csv"

    emitted = 0
    with summary_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["umi", "family_size", "r1_fastq", "r2_fastq"])
        for umi in sorted(families.keys()):
            fam_dir = family_fastq_root / _sanitize_umi_for_path(umi)
            fam_dir.mkdir(parents=True, exist_ok=True)
            r1_out = fam_dir / f"{_sanitize_umi_for_path(umi)}_R1.fastq.gz"
            r2_out = fam_dir / f"{_sanitize_umi_for_path(umi)}_R2.fastq.gz"
            with _open_text_maybe_gz(r1_out, "wt") as o1, _open_text_maybe_gz(r2_out, "wt") as o2:
                for pair in families[umi]:
                    write_fastq_record(o1, pair.r1)
                    write_fastq_record(o2, pair.r2)
            w.writerow([umi, len(families[umi]), str(r1_out), str(r2_out)])
            emitted += 1

    info = {
        "strategy": "family_resolved",
        "pop": pop,
        "family_fastq_root": str(family_fastq_root),
        "family_resolved_summary_csv": str(summary_csv),
        "families_emitted": emitted,
        "note": "Prepared per-family FASTQ files for family-resolved alignment and variant calling.",
    }
    _write_json(strategy_dir / "strategy_output.json", info)
    return info


def stage7_family_strategy(
    exp: str,
    run_tag: str,
    results_exp_dir: Path,
    metrics_dir: Path,
    trimmed_version: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Stage 7:
      - Rebuild UMI families from Stage 5 trimmed FASTQ.
      - Print family-size report in terminal.
      - Ask user to choose exactly one strategy:
          1) reporting
          2) collapse
          3) consensus
          4) family_resolved
      - Emit only the selected branch.
    """
    version, pop_pairs = _find_trimmed_pairs(results_exp_dir, run_tag, trimmed_version)
    stage_root = results_exp_dir / "family_strategy" / run_tag / version
    stage_root.mkdir(parents=True, exist_ok=True)

    payload: Dict[str, Any] = {
        "stage": "stage7_family_strategy",
        "exp": exp,
        "run_tag": run_tag,
        "trimmed_version": version,
        "per_population": {},
    }

    for pop, r1, r2 in pop_pairs:
        pop_root = stage_root / pop
        pop_root.mkdir(parents=True, exist_ok=True)

        families = _rebuild_families(r1, r2)
        summary = _compute_family_summary(families)
        _write_family_table(pop_root / "family_size_table.csv", families)

        family_report = {
            "pop": pop,
            "trimmed_r1": str(r1),
            "trimmed_r2": str(r2),
            "total_families": summary.total_families,
            "total_pairs": summary.total_pairs,
            "size_1": summary.size_1,
            "size_2": summary.size_2,
            "size_ge_3": summary.size_ge_3,
            "mean_family_size": summary.mean_family_size,
            "median_family_size": summary.median_family_size,
        }
        _write_json(pop_root / "family_report.json", family_report)
        _print_family_report(pop, summary)
        strategy = _prompt_strategy_choice()

        strategy_dir = pop_root / strategy

        print(f"\nSelected strategy for pop {pop}: {strategy}\n")
        print(f"Output directory:")
        print(f"  {strategy_dir.resolve()}\n")

        selection = {"selected_strategy": strategy}
        _write_json(pop_root / "family_strategy_selection.json", selection)

        strategy_dir = pop_root / strategy
        strategy_dir.mkdir(parents=True, exist_ok=True)

        if strategy == "reporting":
            info = _run_reporting_strategy(pop, r1, r2, strategy_dir, families)
        elif strategy == "collapse":
            info = _run_collapse_strategy(pop, strategy_dir, families)
        elif strategy == "consensus":
            info = _run_consensus_strategy(pop, strategy_dir, families)
        elif strategy == "family_resolved":
            info = _run_family_resolved_strategy(pop, strategy_dir, families)

        payload["per_population"][pop] = {
            "family_report": family_report,
            "selection": selection,
            "strategy_output": info,
            "family_size_table_csv": str(pop_root / "family_size_table.csv"),
            "family_report_json": str(pop_root / "family_report.json"),
        }

    _write_json(metrics_dir / f"metrics_stage7_family_strategy_{version}.json", payload)
    print(f"Stage7 metrics written to:")
    print(metrics_dir / f"metrics_stage7_family_strategy_{version}.json")
    return payload
