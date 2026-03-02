# src/demultiplex_index.py
"""
Script 2: Sort NGS reads by DNA index

Enhancements:
- --min_len (default 150)
- --output_dir (default results/<exp>/demultiplexing)
- --metrics_json (machine-readable stats; optional)
"""

import os
import csv
import argparse
import glob
import json
from datetime import datetime
from Bio import SeqIO


def find_multiplexing_csv(experiment_name):
    input_dir = os.path.join("input_data", experiment_name)
    csv_files = glob.glob(os.path.join(input_dir, "*multiplexing_info*.csv"))
    if csv_files:
        return csv_files[0]
    raise FileNotFoundError(f"No multiplexing info CSV file found in {input_dir}")


def find_fastq_files(experiment_name, gw_name):
    input_dir = os.path.join("input_data", experiment_name)

    r1_patterns = [
        os.path.join(input_dir, "**", f"{gw_name}_R1*.fastq"),
        os.path.join(input_dir, "**", f"{gw_name}_R1*.fq"),
    ]
    r2_patterns = [
        os.path.join(input_dir, "**", f"{gw_name}_R2*.fastq"),
        os.path.join(input_dir, "**", f"{gw_name}_R2*.fq"),
    ]

    r1_files, r2_files = [], []
    for pattern in r1_patterns:
        r1_files.extend(glob.glob(pattern, recursive=True))
    for pattern in r2_patterns:
        r2_files.extend(glob.glob(pattern, recursive=True))

    if not r1_files or not r2_files:
        return None
    return r1_files[0], r2_files[0]


def load_populations(csv_path):
    populations = []
    with open(csv_path, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            populations.append(
                {
                    "GW_name": row["GW_name"],
                    "Population": row["Population"],
                    "Time": row["Time"],
                    "R1_index": row["R1_index"],
                    "R2_index": row["R2_index"],
                }
            )
    return populations



def process_population_reads(population_info, r1_fastq, r2_fastq, output_dir, min_len):

    gw_name = population_info["GW_name"]
    population = f"P{population_info['Population']}"
    r1_index = population_info["R1_index"]
    r2_index = population_info["R2_index"]

    counts = {"total": 0, "short": 0, "unmatched": 0, "matched": 0}

    pop_folder = os.path.join(output_dir, population)
    os.makedirs(pop_folder, exist_ok=True)

    def fastq_count(path: str) -> int:
        with open(path, "r") as fh:
            return sum(1 for _ in fh) // 4

    n1 = fastq_count(r1_fastq)
    n2 = fastq_count(r2_fastq)

    if n1 != n2:
        raise RuntimeError(
            f"[PAIRING ERROR] R1/R2 read counts differ: R1={n1}, R2={n2}"
        )

    pop_r1_file = os.path.join(pop_folder, f"{population}_R1.fastq")
    pop_r2_file = os.path.join(pop_folder, f"{population}_R2.fastq")

    short_r1_path = os.path.join(output_dir, f"{gw_name}_short_reads_R1.fastq")
    short_r2_path = os.path.join(output_dir, f"{gw_name}_short_reads_R2.fastq")
    unmatched_r1_path = os.path.join(output_dir, f"{gw_name}_unmatched_reads_R1.fastq")
    unmatched_r2_path = os.path.join(output_dir, f"{gw_name}_unmatched_reads_R2.fastq")

    pop_r1_handle = open(pop_r1_file, "w")
    pop_r2_handle = open(pop_r2_file, "w")
    short_r1_handle = open(short_r1_path, "w")
    short_r2_handle = open(short_r2_path, "w")
    unmatched_r1_handle = open(unmatched_r1_path, "w")
    unmatched_r2_handle = open(unmatched_r2_path, "w")

    with open(r1_fastq, "r") as r1_handle, open(r2_fastq, "r") as r2_handle:
        r1_records = SeqIO.parse(r1_handle, "fastq")
        r2_records = SeqIO.parse(r2_handle, "fastq")

        for r1_record, r2_record in zip(r1_records, r2_records):
            counts["total"] += 1

            r1_id = r1_record.id.split()[0]
            r2_id = r2_record.id.split()[0]
            if r1_id != r2_id:
                raise RuntimeError(
                    f"[PAIRING ERROR] R1/R2 read ID mismatch detected for {gw_name} ({population}).\n"
                    f"R1: {r1_id}\n"
                    f"R2: {r2_id}\n"
                    "Paired FASTQ files are not synchronized (possible single-end filtering, file mix-up, or truncation)."
                )

            # Length filter (user-configurable)
            if len(r1_record.seq) < min_len or len(r2_record.seq) < min_len:
                SeqIO.write(r1_record, short_r1_handle, "fastq")
                SeqIO.write(r2_record, short_r2_handle, "fastq")
                counts["short"] += 1
                continue

            read_r1_index = str(r1_record.seq[:8])
            read_r2_index = str(r2_record.seq[:8])

            if read_r1_index == r1_index and read_r2_index == r2_index:
                SeqIO.write(r1_record, pop_r1_handle, "fastq")
                SeqIO.write(r2_record, pop_r2_handle, "fastq")
                counts["matched"] += 1
            else:
                SeqIO.write(r1_record, unmatched_r1_handle, "fastq")
                SeqIO.write(r2_record, unmatched_r2_handle, "fastq")
                counts["unmatched"] += 1

            if counts["total"] % 10000 == 0:
                print(f"  Processed {counts['total']} reads...")

    pop_r1_handle.close()
    pop_r2_handle.close()
    short_r1_handle.close()
    short_r2_handle.close()
    unmatched_r1_handle.close()
    unmatched_r2_handle.close()

    return counts


def demultiplex_all_populations(experiment_name, csv_path, output_dir, min_len):
    populations = load_populations(csv_path)

    print(f"Found {len(populations)} population(s) in CSV\n")

    overall_stats = []
    overall_totals = {"total": 0, "short": 0, "unmatched": 0, "matched": 0}

    for pop_info in populations:
        gw_name = pop_info["GW_name"]
        print(f"Processing population: {gw_name}")
        print(f"  Expected indexes - R1: {pop_info['R1_index']}, R2: {pop_info['R2_index']}")

        fastq_files = find_fastq_files(experiment_name, gw_name)
        if fastq_files is None:
            print(f"  Warning: No fastq files found for {gw_name}, skipping...\n")
            continue

        r1_fastq, r2_fastq = fastq_files
        print(f"  Found R1: {r1_fastq}")
        print(f"  Found R2: {r2_fastq}")

        counts = process_population_reads(pop_info, r1_fastq, r2_fastq, output_dir, min_len)

        population = f"P{pop_info['Population']}"
        print(f"  Summary for {population} (GW_name: {gw_name}):")
        print(f"    Total reads: {counts['total']}")
        print(f"    Matched reads: {counts['matched']} ({counts['matched']/counts['total']*100:.2f}%)")
        print(f"    Short reads: {counts['short']} ({counts['short']/counts['total']*100:.2f}%)")
        print(f"    Unmatched reads: {counts['unmatched']} ({counts['unmatched']/counts['total']*100:.2f}%)\n")

        overall_stats.append({"population": population, "gw_name": gw_name, "counts": counts})

        for k in overall_totals:
            overall_totals[k] += counts[k]

    print("=" * 60)
    print("OVERALL DEMULTIPLEXING SUMMARY")
    print("=" * 60)
    for stat in overall_stats:
        population = stat["population"]
        gw_name = stat["gw_name"]
        counts = stat["counts"]
        print(
            f"{population} (GW_name: {gw_name}):\n"
            f"  Total: {counts['total']}, Matched: {counts['matched']}, "
            f"Short: {counts['short']}, Unmatched: {counts['unmatched']}"
        )
    print("=" * 60)

    # overall discard rate = short / total (your current definition)
    discard_rate = (overall_totals["short"] / overall_totals["total"]) if overall_totals["total"] else 0.0

    return {
        "experiment": experiment_name,
        "min_len": min_len,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "overall": {
            **overall_totals,
            "discard_rate": discard_rate,
            "short_rate": discard_rate,
            "unmatched_rate": (overall_totals["unmatched"] / overall_totals["total"]) if overall_totals["total"] else 0.0,
            "matched_rate": (overall_totals["matched"] / overall_totals["total"]) if overall_totals["total"] else 0.0,
        },
        "per_population": overall_stats,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Demultiplex paired-end NGS reads based on DNA indexes for all populations"
    )
    parser.add_argument("experiment_name", help='Name of the experiment (e.g., "example")')
    parser.add_argument("--min_len", type=int, default=150, help="Minimum read length. Default: 150")
    parser.add_argument(
        "--output_dir",
        default=None,
        help="Output directory. Default: results/<experiment_name>/demultiplexing",
    )
    parser.add_argument(
        "--metrics_json",
        default=None,
        help="Optional path to write machine-readable metrics JSON.",
    )

    args = parser.parse_args()

    try:
        csv_path = find_multiplexing_csv(args.experiment_name)
        print(f"Found multiplexing CSV: {csv_path}")
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return

    output_dir = args.output_dir or os.path.join("results", args.experiment_name, "demultiplexing")
    if not os.path.exists(output_dir):
        print(f"Error: Output directory not found: {output_dir}")
        print("Please run demultiplex_folders.py first to create the directory structure.")
        return

    print(f"Output directory: {output_dir}")
    print(f"Using min_len = {args.min_len}\n")

    metrics = demultiplex_all_populations(args.experiment_name, csv_path, output_dir, args.min_len)

    if args.metrics_json:
        os.makedirs(os.path.dirname(args.metrics_json), exist_ok=True)
        with open(args.metrics_json, "w") as f:
            json.dump(metrics, f, indent=2)
        print(f"\nWrote metrics JSON: {args.metrics_json}")


if __name__ == "__main__":
    main()