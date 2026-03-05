#!/usr/bin/env python3
"""
Script 6: Prepare UMI sequences for alignment

This script takes the UMI dictionaries created by Script 4 (demultiplex_UMI.py),
trims the sequences using stored trim positions, and writes them to FASTQ files
organized by UMI pair for downstream alignment.

Input: UMI dictionary pickle files from results/{experiment}/demultiplexing/P{Population}/
Output: Trimmed sequence files organized by UMI pair in results/{experiment}/alignment/fastq

Dependencies: Script 4 (demultiplex_UMI.py) must be run first
Usage: python alignment_prep.py <experiment_name>
"""

import os
import pickle
import glob
import argparse
import sys
import atexit
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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


def find_umi_libraries(experiment_name, output_base='results'):
    """
    Find all UMI library pickle files for an experiment.
    
    Parameters:
    -----------
    experiment_name : str
        Name of the experiment
    output_base : str
        Base output directory
        
    Returns:
    --------
    dict
        Dictionary mapping population names to UMI library paths
    """
    demultiplex_dir = os.path.join(output_base, experiment_name, 'demultiplexing')
    
    if not os.path.exists(demultiplex_dir):
        raise FileNotFoundError(f"Demultiplexing directory not found: {demultiplex_dir}")
    
    libraries = {}
    pattern = os.path.join(demultiplex_dir, 'P*', '*_UMI_dict.pkl')
    
    for lib_path in sorted(glob.glob(pattern)):
        # Extract population name from path (e.g., P1)
        population = os.path.basename(os.path.dirname(lib_path)).replace('P', '')
        libraries[f"P{population}"] = lib_path
    
    return libraries


def load_umi_dict(library_path):
    """
    Load a UMI dictionary from a pickle file.
    
    Parameters:
    -----------
    library_path : str
        Path to the pickle file
        
    Returns:
    --------
    dict
        UMI dictionary: {(forward_UMI, reverse_UMI): {'R1': [SeqRecord], 'R2': [SeqRecord], 'R1_trim_pos': pos, ...}}
    """
    with open(library_path, 'rb') as f:
        umi_dict = pickle.load(f)
    return umi_dict


def trim_and_prepare_sequences(umi_dict, population, output_dir, min_reads_per_umi=1):
    """
    Trim sequences using stored trim positions and write to FASTQ files by UMI pair.
    
    Parameters:
    -----------
    umi_dict : dict
        UMI dictionary with full untrimmed sequences and trim positions
    population : str
        Population name (e.g., "P1")
    output_dir : str
        Output directory for alignment prep files
    min_reads_per_umi : int
        Minimum number of reads required per UMI pair to be included (default: 2)
        
    Returns:
    --------
    dict
        Statistics about trimming and file generation
    """
    print(f"  Preparing sequences for {population}...")
    
    stats = {
        'total_umis': len(umi_dict),
        'total_reads_trimmed': 0,
        'umi_files_created': 0,
        'umis_skipped_below_threshold': 0
    }
    
    # Create population-specific output directory
    pop_output_dir = os.path.join(output_dir, population)
    os.makedirs(pop_output_dir, exist_ok=True)
    
    for umi_pair, data in umi_dict.items():
        forward_umi, reverse_umi = umi_pair
        
        # Only process UMI pairs that meet the minimum read threshold
        num_sequences = len(data['R1'])
        if num_sequences < min_reads_per_umi:
            stats['umis_skipped_below_threshold'] += 1
            continue
        
        # Get trim positions (should be consistent across all sequences in this UMI pair)
        r1_trim_pos = data['R1_trim_pos']
        r2_trim_pos = data['R2_trim_pos']
        
        # Create FASTQ records with trimmed sequences
        r1_records = []
        r2_records = []

        # Backward-compatible handling for older dictionaries that stored string sequences + IDs
        has_legacy_ids = ('R1_ids' in data and 'R2_ids' in data)

        if has_legacy_ids:
            for r1_seq, r2_seq, r1_id, r2_id in zip(data['R1'], data['R2'], data['R1_ids'], data['R2_ids']):
                trimmed_r1_seq = r1_seq[r1_trim_pos:]
                trimmed_r2_seq = r2_seq[r2_trim_pos:]

                r1_record = SeqRecord(
                    Seq(trimmed_r1_seq),
                    id=r1_id,
                    description=f"UMI:{forward_umi}_{reverse_umi}",
                    letter_annotations={"phred_quality": [40] * len(trimmed_r1_seq)}
                )
                r2_record = SeqRecord(
                    Seq(trimmed_r2_seq),
                    id=r2_id,
                    description=f"UMI:{forward_umi}_{reverse_umi}",
                    letter_annotations={"phred_quality": [40] * len(trimmed_r2_seq)}
                )

                r1_records.append(r1_record)
                r2_records.append(r2_record)
                stats['total_reads_trimmed'] += 1
        else:
            for r1_record, r2_record in zip(data['R1'], data['R2']):
                # Records are SeqRecord objects; slicing preserves per-base quality annotations
                trimmed_r1_record = r1_record[r1_trim_pos:]
                trimmed_r2_record = r2_record[r2_trim_pos:]

                # Keep IDs/quality from source and annotate description with UMI pair
                trimmed_r1_record.description = f"UMI:{forward_umi}_{reverse_umi}"
                trimmed_r2_record.description = f"UMI:{forward_umi}_{reverse_umi}"

                r1_records.append(trimmed_r1_record)
                r2_records.append(trimmed_r2_record)
                stats['total_reads_trimmed'] += 1
        
        # Write trimmed sequences to FASTQ files
        umi_safe_name = f"{forward_umi}_{reverse_umi}"
        r1_output_path = os.path.join(pop_output_dir, f"{umi_safe_name}_R1.fastq")
        r2_output_path = os.path.join(pop_output_dir, f"{umi_safe_name}_R2.fastq")
        
        SeqIO.write(r1_records, r1_output_path, "fastq")
        SeqIO.write(r2_records, r2_output_path, "fastq")
        
        stats['umi_files_created'] += 1
    
    return stats


def process_all_populations_for_alignment(experiment_name, min_reads_per_umi=1):
    """
    Process all populations to prepare sequences for alignment.
    
    Parameters:
    -----------
    experiment_name : str
        Name of the experiment (e.g., "example")
    min_reads_per_umi : int
        Minimum number of reads required per UMI pair to be included (default: 2)
        
    Returns:
    --------
    bool
        True if processing was successful
    """
    # Find UMI libraries
    print("Searching for UMI libraries...")
    try:
        libraries = find_umi_libraries(experiment_name)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return False
    
    if not libraries:
        print("Error: No UMI libraries found")
        return False
    
    print(f"Found {len(libraries)} population(s)\n")
    
    # Create output directory
    output_dir = os.path.join('results', experiment_name, 'alignment', 'fastq')
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}\n")
    
    overall_stats = {}
    
    # Process each population
    for population, lib_path in libraries.items():
        print(f"Processing {population}...")
        try:
            umi_dict = load_umi_dict(lib_path)
            stats = trim_and_prepare_sequences(umi_dict, population, output_dir, min_reads_per_umi)
            overall_stats[population] = stats
            
            print(f"  Summary for {population}:")
            print(f"    Total UMI pairs: {stats['total_umis']}")
            print(f"    UMI pairs with >= {min_reads_per_umi} reads: {stats['umi_files_created']}")
            print(f"    UMI pairs with < {min_reads_per_umi} reads (skipped): {stats['umis_skipped_below_threshold']}")
            print(f"    Reads trimmed and written: {stats['total_reads_trimmed']}\n")
        except Exception as e:
            print(f"  Error processing {population}: {e}\n")
            continue
    
    # Print overall summary
    print("=" * 60)
    print("ALIGNMENT PREPARATION SUMMARY")
    print("=" * 60)
    for population, stats in sorted(overall_stats.items()):
        print(f"{population}:")
        print(f"  Total UMI pairs: {stats['total_umis']}")
        print(f"  UMI pairs with >= {min_reads_per_umi} reads: {stats['umi_files_created']}")
        print(f"  UMI pairs with < {min_reads_per_umi} reads (skipped): {stats['umis_skipped_below_threshold']}")
        print(f"  Reads trimmed: {stats['total_reads_trimmed']}")
    print("=" * 60)
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description='Prepare UMI sequences for alignment by trimming and organizing into FASTQ files'
    )
    parser.add_argument(
        'experiment_name',
        help='Name of the experiment (e.g., "example"). Script will process UMI dictionaries from results/{experiment_name}/demultiplexing/'
    )
    parser.add_argument(
        '--min-reads-per-umi',
        type=int,
        default=2,
        help='Minimum number of reads required per UMI pair before it is included (default: 2)'
    )
    
    args = parser.parse_args()

    if args.min_reads_per_umi < 1:
        parser.error('--min-reads-per-umi must be an integer >= 1')
    
    setup_terminal_logging(args.experiment_name, "alignment_prep")
    
    print(f"\n{'='*80}")
    print(f"Alignment Preparation: {args.experiment_name}")
    print(f"Minimum reads per UMI: {args.min_reads_per_umi}")
    print(f"{'='*80}\n")
    
    success = process_all_populations_for_alignment(args.experiment_name, args.min_reads_per_umi)
    
    if success:
        print(f"\n{'='*80}")
        print("Alignment preparation complete!")
        print(f"{'='*80}\n")
    else:
        print("\nAlignment preparation failed")


if __name__ == '__main__':
    main()
