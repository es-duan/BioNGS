# src/demultiplex_folders.py
"""
Script 1: Generate folders to store population demultiplexed files

Input: Experiment name (looks for CSV in input_data/{experiment_name}/)
Output: Create a folder and empty R1 and R2 fastq files for each population
"""

import os
import csv
import argparse
import glob


def find_multiplexing_csv(experiment_name):
    input_dir = os.path.join("input_data", experiment_name)

    if not os.path.exists(input_dir):
        raise FileNotFoundError(f"Experiment directory not found: {input_dir}")

    csv_files = glob.glob(os.path.join(input_dir, "*multiplexing_info*.csv"))
    if csv_files:
        return csv_files[0]

    raise FileNotFoundError(f"No CSV file found in {input_dir}")


def create_demultiplex_folders(experiment_name, output_dir=None):
    csv_path = find_multiplexing_csv(experiment_name)
    print(f"Found multiplexing CSV: {csv_path}")

    if output_dir is None:
        output_dir = os.path.join("results", experiment_name, "demultiplexing")

    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}\n")

    with open(csv_path, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            gw_name = row["GW_name"]
            population = f"P{row['Population']}"

            pop_folder = os.path.join(output_dir, population)
            os.makedirs(pop_folder, exist_ok=True)

            r1_file = os.path.join(pop_folder, f"{population}_R1.fastq")
            r2_file = os.path.join(pop_folder, f"{population}_R2.fastq")

            open(r1_file, "w").close()
            open(r2_file, "w").close()

            print(f"Created folder and files for population: {population} (GW_name: {gw_name})")
            print(f"  R1 index: {row['R1_index']}, R2 index: {row['R2_index']}")

    print(f"\nAll population folders created in: {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate folders and empty fastq files for each population in an experiment"
    )
    parser.add_argument("experiment_name", help='Name of the experiment (e.g., "example")')
    parser.add_argument(
        "--output_dir",
        default=None,
        help="Optional output directory. Default: results/<experiment_name>/demultiplexing",
    )
    args = parser.parse_args()

    create_demultiplex_folders(args.experiment_name, output_dir=args.output_dir)


if __name__ == "__main__":
    main()