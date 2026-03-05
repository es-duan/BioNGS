"""
Script 1: Generate folders to store population demultiplexed files

Input: Experiment name (looks for CSV in input_data/{experiment_name}/)
Output: Create a folder and empty R1 and R2 fastq files for each population
Description: This is the first script of the pipeline. The user provides an experiment name,
             and the script finds the multiplexing CSV file in the corresponding input_data
             folder and creates the output directory structure in results/{experiment_name}/.
"""

import os
import csv
import argparse
import glob
import sys
import atexit
from datetime import datetime


class TeeStream:
    """Write stream output to both terminal and a log file."""

    def __init__(self, terminal_stream, log_file):
        self.terminal_stream = terminal_stream
        self.log_file = log_file

    def write(self, message):
        self.terminal_stream.write(message)
        self.log_file.write(message)

    def flush(self):
        self.terminal_stream.flush()
        self.log_file.flush()


def setup_terminal_logging(experiment_name, script_name):
    """Save all terminal output for this script run to a text file."""
    log_dir = os.path.join("results", experiment_name, "logs")
    os.makedirs(log_dir, exist_ok=True)

    log_path = os.path.join(log_dir, f"{script_name}_terminal_output.txt")
    log_file = open(log_path, 'w', encoding='utf-8')

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_file.write("\n" + "=" * 80 + "\n")
    log_file.write(f"Run started: {timestamp}\n")
    log_file.write(f"Script: {script_name}\n")
    log_file.write("=" * 80 + "\n")

    original_stdout = sys.stdout
    original_stderr = sys.stderr
    sys.stdout = TeeStream(original_stdout, log_file)
    sys.stderr = TeeStream(original_stderr, log_file)

    def cleanup_logging():
        end_timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log_file.write("\n" + "=" * 80 + "\n")
        log_file.write(f"Run ended: {end_timestamp}\n")
        log_file.write("=" * 80 + "\n")
        log_file.flush()
        sys.stdout = original_stdout
        sys.stderr = original_stderr
        log_file.close()

    atexit.register(cleanup_logging)
    print(f"Logging terminal output to: {log_path}")


def find_multiplexing_csv(experiment_name):
    """
    Find the multiplexing CSV file for a given experiment.
    
    Parameters:
    -----------
    experiment_name : str
        Name of the experiment
        
    Returns:
    --------
    str
        Path to the multiplexing CSV file
        
    Raises:
    -------
    FileNotFoundError
        If no suitable CSV file is found
    """
    input_dir = os.path.join("input_data", experiment_name)
    
    if not os.path.exists(input_dir):
        raise FileNotFoundError(f"Experiment directory not found: {input_dir}")
    
    # Look for CSV files with 'multiplexing' in the name
    csv_files = glob.glob(os.path.join(input_dir, "*multiplexing_info*.csv"))
    
    if csv_files:
        return csv_files[0]
    
    raise FileNotFoundError(f"No CSV file found in {input_dir}")


def create_demultiplex_folders(experiment_name):
    """
    Create folders and empty fastq files for each population in the experiment.
    
    Parameters:
    -----------
    experiment_name : str
        Name of the experiment (used to locate input CSV and create output directory)
    """
    # Find the CSV file
    csv_path = find_multiplexing_csv(experiment_name)
    print(f"Found multiplexing CSV: {csv_path}")
    
    # Create output directory
    output_dir = os.path.join("results", experiment_name, "demultiplexing")
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}\n")
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Read the CSV file
    with open(csv_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        
        for row in reader:
            gw_name = row['GW_name']
            population = f"P{row['Population']}"
            
            # Create folder for this population
            pop_folder = os.path.join(output_dir, population)
            os.makedirs(pop_folder, exist_ok=True)
            
            # Create empty R1 and R2 fastq files
            r1_file = os.path.join(pop_folder, f"{population}_R1.fastq")
            r2_file = os.path.join(pop_folder, f"{population}_R2.fastq")
            
            # Create empty files
            open(r1_file, 'w').close()
            open(r2_file, 'w').close()
            
            print(f"Created folder and files for population: {population} (GW_name: {gw_name})")
            print(f"  R1 index: {row['R1_index']}, R2 index: {row['R2_index']}")
    
    print(f"\nAll population folders created in: {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate folders and empty fastq files for each population in an experiment'
    )
    parser.add_argument(
        'experiment_name',
        help='Name of the experiment (e.g., "example"). Script will look for CSV in input_data/{experiment_name}/ and create outputs in results/{experiment_name}/'
    )
    
    args = parser.parse_args()

    setup_terminal_logging(args.experiment_name, "demultiplex_folders")
    
    create_demultiplex_folders(args.experiment_name)


if __name__ == '__main__':
    main()
