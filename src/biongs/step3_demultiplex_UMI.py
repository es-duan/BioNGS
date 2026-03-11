"""
Step 3: Create a dictionary of UMI reads for each population

Input: Demultiplexed population fastq files, UMI primer sequences from CSV
Output: Dictionary where keys are UMI pairs and values are lists of R1 and R2 SeqRecord
        objects that match that UMI. Dictionaries are saved in respective demultiplexing/population
        folders as pickle files named P{population_number}_UMI_dict
Dependencies: Step 2 (step2_demultiplex_index.py) must be run first
Description: For each population, detect UMI sequences by aligning to primer sequences
             (the UMI is the 10 N bps on the primer). Forward UMIs are detected from R1,
             reverse UMIs from R2. For each unique UMI pair, store matched full SeqRecord
             objects plus trim positions; trimming is done in step4_alignment_prep.py.
"""

import os
import csv
import pickle
import argparse
import glob
import sys
import atexit
from datetime import datetime
from Bio import SeqIO


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


def parse_primer_for_umi(primer_sequence):
    """
    Parse a primer sequence to find UMI position and flanking sequences.
    
    Parameters:
    -----------
    primer_sequence : str
        Primer sequence with 10 N's indicating the UMI position
        
    Returns:
    --------
    dict
        Dictionary with keys:
        - 'before': sequence before the N's
        - 'after': sequence after the N's
        - 'umi_start': position where the 10 N's start
        - 'umi_end': position where the 10 N's end
    """
    primer_upper = primer_sequence.upper()
    
    # Find the consecutive N's
    n_count = 0
    n_start = -1
    
    for i, base in enumerate(primer_upper):
        if base == 'N':
            if n_count == 0:
                n_start = i
            n_count += 1
        else:
            if n_count > 0:
                # Found end of N stretch
                break
    
    if n_count != 10:
        raise ValueError(f"Primer must contain exactly 10 consecutive N's for UMI, found {n_count}")
    
    n_end = n_start + 10
    
    return {
        'before': primer_upper[:n_start],
        'after': primer_upper[n_end:],
        'umi_start': n_start,
        'umi_end': n_end,
        'full_primer': primer_upper
    }


def extract_umi_from_sequence(sequence, primer_info):
    """
    Extract the UMI from a sequence by aligning primer regions.
    
    Parameters:
    -----------
    sequence : str or Bio.Seq.Seq
        The DNA sequence to search
    primer_info : dict
        Dictionary with primer structure info from parse_primer_for_umi
        
    Returns:
    --------
    str or None
        The extracted UMI sequence (10 bp), or None if primer not found
    """
    sequence_str = str(sequence).upper()
    primer_before = primer_info['before']
    primer_after = primer_info['after']
    umi_start = primer_info['umi_start']
    umi_end = primer_info['umi_end']
    
    # Try to find the primer parts in the sequence
    # Look for the "before" part first
    if primer_before:
        before_pos = sequence_str.find(primer_before)
        if before_pos == -1:
            return None
        
        # The UMI should start right after the "before" part
        umi_pos = before_pos + len(primer_before)
    else:
        # If no "before" part, assume UMI starts at position 0
        umi_pos = 0
    
    # Check if the "after" part matches where expected
    if primer_after:
        expected_after_pos = umi_pos + 10
        if expected_after_pos + len(primer_after) > len(sequence_str):
            return None
        
        actual_after = sequence_str[expected_after_pos:expected_after_pos + len(primer_after)]
        if actual_after != primer_after:
            return None
    else:
        # If no "after" part, just check that we have enough sequence
        if umi_pos + 10 > len(sequence_str):
            return None
    
    # Extract the UMI
    umi = sequence_str[umi_pos:umi_pos + 10]
    
    # Calculate the trim position (position after entire primer: "before" + UMI + "after")
    trim_pos = umi_pos + 10 + len(primer_after)
    
    return umi, trim_pos


def load_primers_from_csv(primer_csv_path):
    """
    Load UMI primer sequences from CSV file.
    
    Parameters:
    -----------
    primer_csv_path : str
        Path to the UMI primers CSV file
        
    Returns:
    --------
    dict
        Dictionary with keys 'forward' and 'reverse' containing primer info
    """
    if not os.path.exists(primer_csv_path):
        raise FileNotFoundError(f"Primer CSV not found: {primer_csv_path}")
    
    primers = {}
    
    with open(primer_csv_path, 'r', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            # Strip whitespace from values to handle any encoding issues
            row = {k.strip(): v.strip() if isinstance(v, str) else v for k, v in row.items()}
            # Assuming the CSV has 'f' for forward and 'r' for reverse
            if 'f' in row:
                primers['forward'] = parse_primer_for_umi(row['f'])
            if 'r' in row:
                primers['reverse'] = parse_primer_for_umi(row['r'])
    
    if 'forward' not in primers or 'reverse' not in primers:
        raise ValueError("Primer CSV must contain 'f' (forward) and 'r' (reverse) columns")
    
    return primers


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
    csv_files = glob.glob(os.path.join(input_dir, "*multiplexing_info*.csv"))
    
    if csv_files:
        return csv_files[0]
    
    raise FileNotFoundError(f"No multiplexing CSV file found in {input_dir}")


def find_umi_primers_csv(experiment_name):
    """
    Find the UMI primers CSV file for a given experiment.
    
    Parameters:
    -----------
    experiment_name : str
        Name of the experiment
        
    Returns:
    --------
    str
        Path to the UMI primers CSV file
    """
    input_dir = os.path.join("input_data", experiment_name)
    
    # Look for CSV files with 'UMI' in the name
    csv_files = glob.glob(os.path.join(input_dir, "*UMI_primers*.csv"))
    
    if csv_files:
        return csv_files[0]
    
    raise FileNotFoundError(f"No UMI primers CSV file found in {input_dir}. "
                           f"Expected file with 'UMI' or 'primer' in the name.")


def create_umi_dict(population_folder, population_num, gw_name, forward_primer_info, reverse_primer_info):
    """
    Create a UMI dictionary for a single population by processing its fastq files.
    
    Parameters:
    -----------
    population_folder : str
        Path to the population folder (e.g., results/example/demultiplexing/P1/)
    population_num : str
        Population number (e.g., "1")
    gw_name : str
        GW_name of the population (e.g., "P22R1")
    forward_primer_info : dict
        Parsed forward primer info from parse_primer_for_umi
    reverse_primer_info : dict
        Parsed reverse primer info from parse_primer_for_umi
        
    Returns:
    --------
    dict
        Library with (forward_UMI, reverse_UMI) tuples as keys and
        {'R1': [SeqRecord], 'R2': [SeqRecord], 'R1_trim_pos': pos, 'R2_trim_pos': pos} as values
    """
    population = f"P{population_num}"
    
    # Define fastq file paths
    r1_fastq = os.path.join(population_folder, f"{population}_R1.fastq")
    r2_fastq = os.path.join(population_folder, f"{population}_R2.fastq")
    
    # Define unmatched UMI file paths
    unmatched_r1_fastq = os.path.join(population_folder, f"{population}_unmatched_UMI_R1.fastq")
    unmatched_r2_fastq = os.path.join(population_folder, f"{population}_unmatched_UMI_R2.fastq")
    
    if not os.path.exists(r1_fastq) or not os.path.exists(r2_fastq):
        print(f"  Warning: Fastq files not found for {population}")
        return {}
    
    # Initialize UMI library
    umi_library = {}
    stats = {
        'total_reads': 0,
        'reads_with_umi': 0,
        'unmatched_reads': 0,
        'unique_umis': 0
    }
    
    # Parse both fastq files simultaneously
    print(f"  Processing {population}...")
    
    with open(r1_fastq, 'r') as r1_handle, open(r2_fastq, 'r') as r2_handle, \
         open(unmatched_r1_fastq, 'w') as unmatched_r1_handle, \
         open(unmatched_r2_fastq, 'w') as unmatched_r2_handle:
        r1_records = SeqIO.parse(r1_handle, 'fastq')
        r2_records = SeqIO.parse(r2_handle, 'fastq')
        
        for r1_record, r2_record in zip(r1_records, r2_records):
            stats['total_reads'] += 1
            
            # Verify that R1 and R2 match by checking headers
            r1_id = r1_record.id.split()[0]
            r2_id = r2_record.id.split()[0]
            
            if r1_id != r2_id:
                print(f"    Warning: R1 and R2 IDs do not match for {r1_id}")
                continue
            
            # Extract UMIs from both reads
            r1_result = extract_umi_from_sequence(r1_record.seq, forward_primer_info)
            r2_result = extract_umi_from_sequence(r2_record.seq, reverse_primer_info)
            
            if r1_result is None or r2_result is None:
                # Save unmatched sequences to separate fastq files
                SeqIO.write(r1_record, unmatched_r1_handle, 'fastq')
                SeqIO.write(r2_record, unmatched_r2_handle, 'fastq')
                stats['unmatched_reads'] += 1
                continue
            
            forward_umi, r1_trim_pos = r1_result
            reverse_umi, r2_trim_pos = r2_result
            
            stats['reads_with_umi'] += 1
            
            # Create UMI pair key
            umi_pair = (forward_umi, reverse_umi)
            
            # Initialize UMI entry if it doesn't exist
            if umi_pair not in umi_library:
                umi_library[umi_pair] = {
                    'R1': [],
                    'R2': [],
                    'R1_trim_pos': r1_trim_pos,
                    'R2_trim_pos': r2_trim_pos
                }
                stats['unique_umis'] += 1
            
            # Store full SeqRecord objects (untrimmed) with original ID and quality
            # Trimming will be done in alignment_prep.py
            umi_library[umi_pair]['R1'].append(r1_record)
            umi_library[umi_pair]['R2'].append(r2_record)
            
            # Print progress every 10000 reads
            if stats['total_reads'] % 10000 == 0:
                print(f"    Processed {stats['total_reads']} reads...")
    
    # Validate that R1 and R2 lists have matching IDs for each UMI
    for umi_pair, data in umi_library.items():
        r1_ids = [record.id.split()[0] for record in data['R1']]
        r2_ids = [record.id.split()[0] for record in data['R2']]
        if r1_ids != r2_ids:
            print(f"    Warning: R1 and R2 IDs don't match for UMI {umi_pair}")
    
    # Print summary
    print(f"    Summary for {population}:")
    print(f"      Total reads: {stats['total_reads']}")
    print(f"      Reads with UMI: {stats['reads_with_umi']}")
    print(f"      Unmatched reads: {stats['unmatched_reads']}")
    print(f"      Unique UMI pairs: {stats['unique_umis']}")
    if stats['total_reads'] > 0 and stats['unique_umis'] > 0:
        print(f"      Mean reads per UMI: {stats['reads_with_umi'] / stats['unique_umis']:.1f}")
    
    return umi_library


def process_all_populations(experiment_name):
    """
    Process all populations in an experiment to create UMI libraries.
    
    Parameters:
    -----------
    experiment_name : str
        Name of the experiment (e.g., "example")
        
    Returns:
    --------
    bool
        True if at least one UMI library was successfully saved, False otherwise
    """
    # Base output directory
    output_base = os.path.join("results", experiment_name, "demultiplexing")
    
    if not os.path.exists(output_base):
        print(f"Error: Output directory not found: {output_base}")
        print("Please run demultiplex_folders.py and demultiplex_index.py first.")
        return
    
    # Find and load multiplexing CSV
    try:
        multiplexing_csv_path = find_multiplexing_csv(experiment_name)
        print(f"Found multiplexing CSV: {multiplexing_csv_path}")
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return
    
    # Find and load UMI primers CSV
    try:
        umi_primers_csv_path = find_umi_primers_csv(experiment_name)
        print(f"Found UMI primers CSV: {umi_primers_csv_path}")
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return
    
    # Load primers from CSV
    try:
        primers = load_primers_from_csv(umi_primers_csv_path)
        print(f"  Forward primer: {primers['forward']['full_primer']}")
        print(f"  Reverse primer: {primers['reverse']['full_primer']}\n")
    except Exception as e:
        print(f"Error loading primers: {e}")
        return
    
    # Read populations from CSV
    populations = []
    with open(multiplexing_csv_path, 'r', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            populations.append({
                'GW_name': row['GW_name'],
                'Population': row['Population'],
            })
    
    print(f"Found {len(populations)} population(s) in CSV")
    print(f"Output directory: {output_base}\n")
    
    # Track if any libraries were successfully saved
    saved_count = 0
    
    # Process each population
    for pop_info in populations:
        gw_name = pop_info['GW_name']
        population_num = pop_info['Population']
        population = f"P{population_num}"
        
        population_folder = os.path.join(output_base, population)
        
        if not os.path.exists(population_folder):
            print(f"Warning: Population folder not found: {population_folder}")
            continue
        
        # Create UMI dictionary
        umi_dict = create_umi_dict(
            population_folder,
            population_num,
            gw_name,
            primers['forward'],
            primers['reverse']
        )
        
        # Save dictionary as pickle file
        lib_filename = f"{population}_UMI_dict.pkl"
        lib_path = os.path.join(population_folder, lib_filename)
        
        with open(lib_path, 'wb') as f:
            pickle.dump(umi_dict, f)
        
        print(f"  Saved UMI dictionary: {lib_filename}\n")
        saved_count += 1
    
    return saved_count > 0


def main():
    parser = argparse.ArgumentParser(
        description='Create UMI libraries for each population from demultiplexed fastq files'
    )
    parser.add_argument(
        'experiment_name',
        help='Name of the experiment (e.g., "example"). Script will auto-detect multiplexing and UMI primers CSVs from input_data/{experiment_name}/'
    )
    
    args = parser.parse_args()

    setup_terminal_logging(args.experiment_name, "demultiplex_UMI")
    
    success = process_all_populations(args.experiment_name)
    
    if success:
        print("="*60)
        print("UMI dictionary creation complete!")
        print("="*60)


if __name__ == '__main__':
    main()
