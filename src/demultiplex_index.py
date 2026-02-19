"""
Script 2: Sort NGS reads by DNA index

Input: Experiment name
Output: For each population in the multiplexing CSV, process its corresponding fastq files
        (matched by GW_name) and extract reads with matching indexes; create additional 
        files for short reads (<150 bp) and unmatched reads
Dependencies: Script 1 (demultiplex_folders.py) must be run first
Description: Using biopython, for each population in the CSV, find fastq files matching
             its GW_name, detect the forward (first 8 bp of R1) and reverse (first 8 bp 
             of R2) indexes, and extract reads that match that population's indexes.
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
    
    # Look for CSV files with 'multiplexing_info' in the name
    csv_files = glob.glob(os.path.join(input_dir, "*multiplexing_info*.csv"))
    
    if csv_files:
        return csv_files[0]
    
    raise FileNotFoundError(f"No multiplexing info CSV file found in {input_dir}")


def find_fastq_files(experiment_name, gw_name):
    """
    Find R1 and R2 fastq files for a given experiment and GW_name.
    
    Parameters:
    -----------
    experiment_name : str
        Name of the experiment
    gw_name : str
        GW_name of the population (e.g., "P22R1")
        
    Returns:
    --------
    tuple or None
        (r1_path, r2_path) paths to the R1 and R2 fastq files, or None if not found
    """
    input_dir = os.path.join("input_data", experiment_name)
    
    # Search for R1 and R2 files recursively
    r1_patterns = [
        os.path.join(input_dir, "**", f"{gw_name}_R1*.fastq"),
        os.path.join(input_dir, "**", f"{gw_name}_R1*.fq"),
    ]
    r2_patterns = [
        os.path.join(input_dir, "**", f"{gw_name}_R2*.fastq"),
        os.path.join(input_dir, "**", f"{gw_name}_R2*.fq"),
    ]
    
    r1_files = []
    r2_files = []
    
    for pattern in r1_patterns:
        r1_files.extend(glob.glob(pattern, recursive=True))
    for pattern in r2_patterns:
        r2_files.extend(glob.glob(pattern, recursive=True))
    
    if not r1_files or not r2_files:
        return None
    
    return r1_files[0], r2_files[0]


def load_populations(csv_path):
    """
    Load all populations and their information from CSV.
    
    Parameters:
    -----------
    csv_path : str
        Path to CSV file with population and index information
        
    Returns:
    --------
    list
        List of dictionaries containing population information
    """
    populations = []
    
    with open(csv_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            populations.append({
                'GW_name': row['GW_name'],
                'Population': row['Population'],
                'Time': row['Time'],
                'R1_index': row['R1_index'],
                'R2_index': row['R2_index']
            })
    
    return populations


def process_population_reads(population_info, r1_fastq, r2_fastq, output_dir):
    """
    Process paired-end reads for a single population, extracting reads that match
    the population's indexes.
    
    Parameters:
    -----------
    population_info : dict
        Dictionary with GW_name, Population, Time, R1_index, R2_index
    r1_fastq : str
        Path to R1 input fastq file
    r2_fastq : str
        Path to R2 input fastq file
    output_dir : str
        Directory containing population folders
        
    Returns:
    --------
    dict
        Dictionary with counts of total, short, matched, and unmatched reads
    """
    gw_name = population_info['GW_name']
    population = f"P{population_info['Population']}"
    r1_index = population_info['R1_index']
    r2_index = population_info['R2_index']
    
    # Initialize counters
    counts = {
        'total': 0,
        'short': 0,
        'unmatched': 0,
        'matched': 0
    }
    
    # Set up output files
    pop_folder = os.path.join(output_dir, population)
    pop_r1_file = os.path.join(pop_folder, f"{population}_R1.fastq")
    pop_r2_file = os.path.join(pop_folder, f"{population}_R2.fastq")
    
    short_r1_path = os.path.join(output_dir, f"{gw_name}_short_reads_R1.fastq")
    short_r2_path = os.path.join(output_dir, f"{gw_name}_short_reads_R2.fastq")
    unmatched_r1_path = os.path.join(output_dir, f"{gw_name}_unmatched_reads_R1.fastq")
    unmatched_r2_path = os.path.join(output_dir, f"{gw_name}_unmatched_reads_R2.fastq")
    
    # Open file handles
    pop_r1_handle = open(pop_r1_file, 'w')
    pop_r2_handle = open(pop_r2_file, 'w')
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
            read_r1_index = str(r1_record.seq[:8])
            read_r2_index = str(r2_record.seq[:8])
            
            # Check if indexes match this population
            if read_r1_index == r1_index and read_r2_index == r2_index:
                # Write to the population files
                SeqIO.write(r1_record, pop_r1_handle, 'fastq')
                SeqIO.write(r2_record, pop_r2_handle, 'fastq')
                counts['matched'] += 1
            else:
                # Write to unmatched files
                SeqIO.write(r1_record, unmatched_r1_handle, 'fastq')
                SeqIO.write(r2_record, unmatched_r2_handle, 'fastq')
                counts['unmatched'] += 1
            
            # Print progress every 10000 reads
            if counts['total'] % 10000 == 0:
                print(f"  Processed {counts['total']} reads...")
    
    # Close all file handles
    pop_r1_handle.close()
    pop_r2_handle.close()
    short_r1_handle.close()
    short_r2_handle.close()
    unmatched_r1_handle.close()
    unmatched_r2_handle.close()
    
    return counts


def demultiplex_all_populations(experiment_name, csv_path, output_dir):
    """
    Process all populations in the CSV, finding and demultiplexing their fastq files.
    
    Parameters:
    -----------
    experiment_name : str
        Name of the experiment
    csv_path : str
        Path to the multiplexing CSV file
    output_dir : str
        Output directory for demultiplexed files
    """
    # Load all populations
    populations = load_populations(csv_path)
    
    print(f"Found {len(populations)} population(s) in CSV\n")
    
    overall_stats = []
    
    for pop_info in populations:
        gw_name = pop_info['GW_name']
        print(f"Processing population: {gw_name}")
        print(f"  Expected indexes - R1: {pop_info['R1_index']}, R2: {pop_info['R2_index']}")
        
        # Find fastq files for this population
        fastq_files = find_fastq_files(experiment_name, gw_name)
        
        if fastq_files is None:
            print(f"  Warning: No fastq files found for {gw_name}, skipping...")
            print()
            continue
        
        r1_fastq, r2_fastq = fastq_files
        print(f"  Found R1: {r1_fastq}")
        print(f"  Found R2: {r2_fastq}")
        
        # Process the reads
        counts = process_population_reads(pop_info, r1_fastq, r2_fastq, output_dir)
        
        # Print summary for this population
        population = f"P{pop_info['Population']}"
        print(f"  Summary for {population} (GW_name: {gw_name}):")
        print(f"    Total reads: {counts['total']}")
        print(f"    Matched reads: {counts['matched']} ({counts['matched']/counts['total']*100:.2f}%)")
        print(f"    Short reads: {counts['short']} ({counts['short']/counts['total']*100:.2f}%)")
        print(f"    Unmatched reads: {counts['unmatched']} ({counts['unmatched']/counts['total']*100:.2f}%)")
        print()
        
        overall_stats.append({
            'population': population,
            'gw_name': gw_name,
            'counts': counts
        })
    
    # Print overall summary
    print("="*60)
    print("OVERALL DEMULTIPLEXING SUMMARY")
    print("="*60)
    for stat in overall_stats:
        population = stat['population']
        gw_name = stat['gw_name']
        counts = stat['counts']
        print(f"{population} (GW_name: {gw_name}):")
        print(f"  Total: {counts['total']}, Matched: {counts['matched']}, "
              f"Short: {counts['short']}, Unmatched: {counts['unmatched']}")
    print("="*60)


def main():
    parser = argparse.ArgumentParser(
        description='Demultiplex paired-end NGS reads based on DNA indexes for all populations'
    )
    parser.add_argument(
        'experiment_name',
        help='Name of the experiment (e.g., "example")'
    )
    
    args = parser.parse_args()
    
    # Find the CSV file
    try:
        csv_path = find_multiplexing_csv(args.experiment_name)
        print(f"Found multiplexing CSV: {csv_path}")
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
    
    # Process all populations
    demultiplex_all_populations(args.experiment_name, csv_path, output_dir)


if __name__ == '__main__':
    main()
