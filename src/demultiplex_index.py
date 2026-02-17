"""
Script 2: Sort NGS reads by DNA index

Input: Experiment name and fastq file basename
Output: Populate the empty fastq files with reads that belong to each population based on
        matched index sequences; create additional fastq files for short reads (<150 bp)
        and unmatched reads
Dependencies: Script 1 (demultiplex_folders.py) must be run first
Description: Using biopython, loop through each pair of input fastq files, detect the
             forward (first 8 bp of R1) and reverse (first 8 bp of R2) indexes, and
             match those with populations listed in the multiplexing_info.csv
"""

import os
import csv
import argparse
import glob
from Bio import SeqIO


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
    """
    input_dir = os.path.join("input_data", experiment_name)
    
    # Look for CSV files with 'multiplexing' in the name
    csv_files = glob.glob(os.path.join(input_dir, "*multiplexing*.csv"))
    
    if csv_files:
        return csv_files[0]
    
    # If no multiplexing CSV found, look for any CSV file
    csv_files = glob.glob(os.path.join(input_dir, "*.csv"))
    
    if csv_files:
        return csv_files[0]
    
    raise FileNotFoundError(f"No CSV file found in {input_dir}")


def find_fastq_files(experiment_name, fastq_basename):
    """
    Find R1 and R2 fastq files for a given experiment and basename.
    
    Parameters:
    -----------
    experiment_name : str
        Name of the experiment
    fastq_basename : str
        Base name of the fastq files (e.g., "P22R1")
        
    Returns:
    --------
    tuple
        (r1_path, r2_path) paths to the R1 and R2 fastq files
    """
    input_dir = os.path.join("input_data", experiment_name)
    
    # Search for R1 and R2 files recursively
    r1_patterns = [
        os.path.join(input_dir, "**", f"{fastq_basename}_R1*.fastq"),
        os.path.join(input_dir, "**", f"{fastq_basename}_R1*.fq"),
    ]
    r2_patterns = [
        os.path.join(input_dir, "**", f"{fastq_basename}_R2*.fastq"),
        os.path.join(input_dir, "**", f"{fastq_basename}_R2*.fq"),
    ]
    
    r1_files = []
    r2_files = []
    
    for pattern in r1_patterns:
        r1_files.extend(glob.glob(pattern, recursive=True))
    for pattern in r2_patterns:
        r2_files.extend(glob.glob(pattern, recursive=True))
    
    if not r1_files:
        raise FileNotFoundError(f"No R1 fastq file found for {fastq_basename} in {input_dir}")
    if not r2_files:
        raise FileNotFoundError(f"No R2 fastq file found for {fastq_basename} in {input_dir}")
    
    return r1_files[0], r2_files[0]


def load_population_indexes(csv_path):
    """
    Load population names and their corresponding indexes from CSV.
    
    Parameters:
    -----------
    csv_path : str
        Path to CSV file with population and index information
        
    Returns:
    --------
    dict
        Dictionary with (R1_index, R2_index) tuples as keys and population info as values
    """
    index_map = {}
    
    with open(csv_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            key = (row['R1_index'], row['R2_index'])
            index_map[key] = {
                'GW_name': row['GW_name'],
                'Population': row['Population'],
                'Time': row['Time']
            }
    
    return index_map


def demultiplex_reads(r1_fastq, r2_fastq, csv_path, output_dir):
    """
    Demultiplex paired-end reads based on their index sequences.
    
    Parameters:
    -----------
    r1_fastq : str
        Path to R1 input fastq file
    r2_fastq : str
        Path to R2 input fastq file
    csv_path : str
        Path to CSV file with population and index information
    output_dir : str
        Directory containing population folders with empty fastq files
    """
    # Load the index mapping
    index_map = load_population_indexes(csv_path)
    
    # Initialize counters
    counts = {
        'total': 0,
        'short': 0,
        'unmatched': 0,
        'matched': {}
    }
    
    # Extract base name for short/unmatched files
    base_name = os.path.splitext(os.path.basename(r1_fastq))[0]
    base_name = base_name.replace('_R1_001', '').replace('_R1', '')
    
    # Create short reads and unmatched reads files
    short_r1_path = os.path.join(output_dir, f"{base_name}_short_reads_R1.fastq")
    short_r2_path = os.path.join(output_dir, f"{base_name}_short_reads_R2.fastq")
    unmatched_r1_path = os.path.join(output_dir, f"{base_name}_unmatched_reads_R1.fastq")
    unmatched_r2_path = os.path.join(output_dir, f"{base_name}_unmatched_reads_R2.fastq")
    
    # Open file handles for population files
    pop_file_handles = {}
    for (r1_idx, r2_idx), pop_info in index_map.items():
        gw_name = pop_info['GW_name']
        pop_folder = os.path.join(output_dir, gw_name)
        r1_file = os.path.join(pop_folder, f"{gw_name}_R1.fastq")
        r2_file = os.path.join(pop_folder, f"{gw_name}_R2.fastq")
        
        pop_file_handles[(r1_idx, r2_idx)] = {
            'R1': open(r1_file, 'a'),
            'R2': open(r2_file, 'a')
        }
        counts['matched'][gw_name] = 0
    
    # Open short and unmatched files
    short_r1_handle = open(short_r1_path, 'w')
    short_r2_handle = open(short_r2_path, 'w')
    unmatched_r1_handle = open(unmatched_r1_path, 'w')
    unmatched_r2_handle = open(unmatched_r2_path, 'w')
    
    # Parse both fastq files simultaneously
    with open(r1_fastq, 'r') as r1_handle, open(r2_fastq, 'r') as r2_handle:
        r1_records = SeqIO.parse(r1_handle, 'fastq')
        r2_records = SeqIO.parse(r2_handle, 'fastq')
        
        for r1_record, r2_record in zip(r1_records, r2_records):
            counts['total'] += 1
            
            # Verify that R1 and R2 match by checking headers
            r1_id = r1_record.id.split()[0]
            r2_id = r2_record.id.split()[0]
            
            if r1_id != r2_id:
                print(f"Warning: R1 and R2 IDs do not match!")
                print(f"  R1: {r1_id}")
                print(f"  R2: {r2_id}")
                continue
            
            # Check if read is too short (< 150 bp)
            if len(r1_record.seq) < 150 or len(r2_record.seq) < 150:
                SeqIO.write(r1_record, short_r1_handle, 'fastq')
                SeqIO.write(r2_record, short_r2_handle, 'fastq')
                counts['short'] += 1
                continue
            
            # Extract the first 8 base pairs as indexes
            r1_index = str(r1_record.seq[:8])
            r2_index = str(r2_record.seq[:8])
            
            # Check if indexes match a population
            index_key = (r1_index, r2_index)
            
            if index_key in pop_file_handles:
                # Write to the appropriate population files
                SeqIO.write(r1_record, pop_file_handles[index_key]['R1'], 'fastq')
                SeqIO.write(r2_record, pop_file_handles[index_key]['R2'], 'fastq')
                
                gw_name = index_map[index_key]['GW_name']
                counts['matched'][gw_name] += 1
            else:
                # Write to unmatched files
                SeqIO.write(r1_record, unmatched_r1_handle, 'fastq')
                SeqIO.write(r2_record, unmatched_r2_handle, 'fastq')
                counts['unmatched'] += 1
            
            # Print progress every 10000 reads
            if counts['total'] % 10000 == 0:
                print(f"Processed {counts['total']} reads...")
    
    # Close all file handles
    for handles in pop_file_handles.values():
        handles['R1'].close()
        handles['R2'].close()
    
    short_r1_handle.close()
    short_r2_handle.close()
    unmatched_r1_handle.close()
    unmatched_r2_handle.close()
    
    # Print summary
    print("\n" + "="*60)
    print("DEMULTIPLEXING SUMMARY")
    print("="*60)
    print(f"Total reads processed: {counts['total']}")
    print(f"Short reads (< 150 bp): {counts['short']} ({counts['short']/counts['total']*100:.2f}%)")
    print(f"Unmatched reads: {counts['unmatched']} ({counts['unmatched']/counts['total']*100:.2f}%)")
    print(f"\nMatched reads by population:")
    for gw_name, count in counts['matched'].items():
        print(f"  {gw_name}: {count} ({count/counts['total']*100:.2f}%)")
    print("="*60)


def main():
    parser = argparse.ArgumentParser(
        description='Demultiplex paired-end NGS reads based on DNA indexes'
    )
    parser.add_argument(
        'experiment_name',
        help='Name of the experiment (e.g., "example")'
    )
    parser.add_argument(
        'fastq_basename',
        help='Base name of the fastq files (e.g., "P22R1" for P22R1_R1_001.fastq and P22R1_R2_001.fastq)'
    )
    
    args = parser.parse_args()
    
    # Find the CSV file
    try:
        csv_path = find_multiplexing_csv(args.experiment_name)
        print(f"Found multiplexing CSV: {csv_path}")
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return
    
    # Find the fastq files
    try:
        r1_fastq, r2_fastq = find_fastq_files(args.experiment_name, args.fastq_basename)
        print(f"Found R1 fastq: {r1_fastq}")
        print(f"Found R2 fastq: {r2_fastq}")
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return
    
    # Set output directory
    output_dir = os.path.join("Outputs", args.experiment_name, "demultiplexing")
    
    if not os.path.exists(output_dir):
        print(f"Error: Output directory not found: {output_dir}")
        print("Please run demultiplex_folders.py first to create the directory structure.")
        return
    
    print(f"Output directory: {output_dir}\n")
    
    demultiplex_reads(r1_fastq, r2_fastq, csv_path, output_dir)


if __name__ == '__main__':
    main()
