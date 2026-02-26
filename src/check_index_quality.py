#!/usr/bin/env python3
"""
Script 3: Check quality of indexing/NGS data

This script analyzes the quality of demultiplexing results from Script 2 (demultiplex_index.py).
It generates visualizations of read length distributions and read counts across populations.

Input: Input fastq files, populated demultiplexed population fastq files, 
       short_read fastq files, and unmatched_read fastq files
Output: Histogram of read lengths from input fastq files; bar plot with the number 
        of reads in each population fastq, short_read, and unmatched_read files
Dependencies: Script 2 (demultiplex_index.py) must be run first
Usage: python check_index_quality.py <experiment_name>
"""

import os
import csv
import glob
import argparse
from collections import defaultdict
from Bio import SeqIO
import altair as alt
import pandas as pd


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


def find_input_fastq_files(experiment_name):
    """
    Find all input fastq files for an experiment.
    
    Parameters:
    -----------
    experiment_name : str
        Name of the experiment
        
    Returns:
    --------
    list
        List of tuples containing (gw_name, r1_path, r2_path)
    """
    input_dir = os.path.join("input_data", experiment_name)
    
    # Find all R1 files
    r1_files = []
    r1_files.extend(glob.glob(os.path.join(input_dir, "**", "*_R1*.fastq"), recursive=True))
    r1_files.extend(glob.glob(os.path.join(input_dir, "**", "*_R1*.fq"), recursive=True))
    
    fastq_pairs = []
    
    for r1_file in r1_files:
        # Extract base name (remove _R1_001.fastq or similar)
        basename = os.path.basename(r1_file)
        gw_name = basename.split('_R1')[0]
        
        # Find corresponding R2 file
        r2_file = r1_file.replace('_R1', '_R2')
        
        if os.path.exists(r2_file):
            fastq_pairs.append((gw_name, r1_file, r2_file))
    
    return fastq_pairs


def count_reads_in_fastq(fastq_path):
    """
    Count the number of reads in a fastq file.
    
    Parameters:
    -----------
    fastq_path : str
        Path to the fastq file
        
    Returns:
    --------
    int
        Number of reads in the file
    """
    if not os.path.exists(fastq_path):
        return 0
    
    count = 0
    with open(fastq_path, 'r') as f:
        for record in SeqIO.parse(f, 'fastq'):
            count += 1
    return count


def get_read_lengths(fastq_path, max_reads=10000):
    """
    Get read lengths from a fastq file.
    
    Parameters:
    -----------
    fastq_path : str
        Path to the fastq file
    max_reads : int
        Maximum number of reads to sample (for large files)
        
    Returns:
    --------
    list
        List of read lengths
    """
    if not os.path.exists(fastq_path):
        return []
    
    lengths = []
    with open(fastq_path, 'r') as f:
        for i, record in enumerate(SeqIO.parse(f, 'fastq')):
            if i >= max_reads:
                break
            lengths.append(len(record.seq))
    return lengths


def load_populations_from_csv(csv_path):
    """
    Load population information from multiplexing CSV.
    
    Parameters:
    -----------
    csv_path : str
        Path to multiplexing CSV file
        
    Returns:
    --------
    list
        List of dictionaries with population info
    """
    populations = []
    
    with open(csv_path, 'r', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            # Strip whitespace from all values
            row = {k.strip(): v.strip() if isinstance(v, str) else v for k, v in row.items()}
            populations.append(row)
    
    return populations


def analyze_demultiplexing_results(experiment_name):
    """
    Analyze demultiplexing results for an experiment.
    
    Parameters:
    -----------
    experiment_name : str
        Name of the experiment
        
    Returns:
    --------
    dict
        Dictionary with analysis results
    """
    output_base = os.path.join("Outputs", experiment_name, "demultiplexing")
    
    if not os.path.exists(output_base):
        raise FileNotFoundError(f"Demultiplexing output not found: {output_base}")
    
    # Load populations from CSV
    csv_path = find_multiplexing_csv(experiment_name)
    populations = load_populations_from_csv(csv_path)
    
    # Find input fastq files
    input_files = find_input_fastq_files(experiment_name)
    
    results = {}
    
    # For each input file (GW_name), analyze its demultiplexing results
    for gw_name, r1_input, r2_input in input_files:
        print(f"\nAnalyzing results for {gw_name}...")
        
        file_results = {
            'populations': {},
            'short_reads': 0,
            'unmatched_reads': 0,
            'input_total': 0
        }
        
        # Count reads in input files
        print(f"  Counting input reads...")
        input_r1_count = count_reads_in_fastq(r1_input)
        file_results['input_total'] = input_r1_count
        
        # Count reads in population files
        for pop in populations:
            if pop['GW_name'] == gw_name:
                pop_num = pop['Population']
                pop_folder = os.path.join(output_base, f"P{pop_num}")
                pop_r1 = os.path.join(pop_folder, f"P{pop_num}_R1.fastq")
                
                if os.path.exists(pop_r1):
                    count = count_reads_in_fastq(pop_r1)
                    file_results['populations'][f"P{pop_num}"] = count
                    print(f"    P{pop_num}: {count:,} reads")
        
        # Count short reads
        short_r1 = os.path.join(output_base, f"{gw_name}_short_reads_R1.fastq")
        if os.path.exists(short_r1):
            count = count_reads_in_fastq(short_r1)
            file_results['short_reads'] = count
            print(f"    Short reads: {count:,} reads")
        
        # Count unmatched reads
        unmatched_r1 = os.path.join(output_base, f"{gw_name}_unmatched_reads_R1.fastq")
        if os.path.exists(unmatched_r1):
            count = count_reads_in_fastq(unmatched_r1)
            file_results['unmatched_reads'] = count
            print(f"    Unmatched reads: {count:,} reads")
        
        # Calculate total distributed reads
        total_distributed = (sum(file_results['populations'].values()) + 
                            file_results['short_reads'] + 
                            file_results['unmatched_reads'])
        print(f"  Total input: {file_results['input_total']:,}")
        print(f"  Total distributed: {total_distributed:,}")
        
        if file_results['input_total'] > 0:
            recovery_rate = (total_distributed / file_results['input_total']) * 100
            print(f"  Recovery rate: {recovery_rate:.1f}%")
        
        results[gw_name] = file_results
    
    return results, input_files


def create_read_length_histogram(experiment_name, input_files, output_dir):
    """
    Create histogram of read lengths from input fastq files.
    
    Parameters:
    -----------
    experiment_name : str
        Name of the experiment
    input_files : list
        List of tuples (gw_name, r1_path, r2_path)
    output_dir : str
        Directory to save plots
    """
    print("\nGenerating read length histogram...")
    
    for gw_name, r1_path, r2_path in input_files:
        print(f"  Sampling read lengths from {gw_name}...")
        
        # Sample reads from both R1 and R2
        r1_lengths = get_read_lengths(r1_path, max_reads=10000)
        r2_lengths = get_read_lengths(r2_path, max_reads=10000)
        
        # Create DataFrames
        r1_df = pd.DataFrame({
            'Read Length (bp)': r1_lengths,
            'Read Type': 'R1'
        })
        r2_df = pd.DataFrame({
            'Read Length (bp)': r2_lengths,
            'Read Type': 'R2'
        })
        
        # Calculate statistics
        r1_mean = r1_df['Read Length (bp)'].mean() if len(r1_df) > 0 else 0
        r1_median = r1_df['Read Length (bp)'].median() if len(r1_df) > 0 else 0
        r2_mean = r2_df['Read Length (bp)'].mean() if len(r2_df) > 0 else 0
        r2_median = r2_df['Read Length (bp)'].median() if len(r2_df) > 0 else 0
        
        # R1 histogram
        if len(r1_df) > 0:
            hist_r1 = alt.Chart(r1_df).mark_bar(
                opacity=0.7,
                color='steelblue'
            ).encode(
                x=alt.X('Read Length (bp):Q', bin=alt.Bin(maxbins=50), title='Read Length (bp)'),
                y=alt.Y('count()', title='Frequency')
            ).properties(
                width=400,
                height=300,
                title=f'{gw_name} - R1 Read Lengths'
            )
            
            # Add mean and median lines
            mean_line_r1 = alt.Chart(pd.DataFrame({'x': [r1_mean]})).mark_rule(
                color='red',
                strokeDash=[5, 5],
                size=2
            ).encode(
                x='x:Q'
            )
            
            median_line_r1 = alt.Chart(pd.DataFrame({'x': [r1_median]})).mark_rule(
                color='orange',
                strokeDash=[5, 5],
                size=2
            ).encode(
                x='x:Q'
            )
            
            # Add text annotations
            text_r1 = alt.Chart(pd.DataFrame({
                'x': [r1_mean, r1_median],
                'y': [0, 0],
                'label': [f'Mean: {r1_mean:.1f} bp', f'Median: {r1_median:.1f} bp'],
                'color': ['red', 'orange']
            })).mark_text(
                align='left',
                dx=5,
                dy=-5,
                fontSize=11
            ).encode(
                x='x:Q',
                text='label:N',
                color=alt.Color('color:N', scale=None)
            )
            
            chart_r1 = hist_r1 + mean_line_r1 + median_line_r1
        else:
            chart_r1 = alt.Chart(pd.DataFrame({'x': [0]})).mark_text(
                text='No data available'
            ).properties(width=400, height=300)
        
        # R2 histogram
        if len(r2_df) > 0:
            hist_r2 = alt.Chart(r2_df).mark_bar(
                opacity=0.7,
                color='seagreen'
            ).encode(
                x=alt.X('Read Length (bp):Q', bin=alt.Bin(maxbins=50), title='Read Length (bp)'),
                y=alt.Y('count()', title='Frequency')
            ).properties(
                width=400,
                height=300,
                title=f'{gw_name} - R2 Read Lengths'
            )
            
            # Add mean and median lines
            mean_line_r2 = alt.Chart(pd.DataFrame({'x': [r2_mean]})).mark_rule(
                color='red',
                strokeDash=[5, 5],
                size=2
            ).encode(
                x='x:Q'
            )
            
            median_line_r2 = alt.Chart(pd.DataFrame({'x': [r2_median]})).mark_rule(
                color='orange',
                strokeDash=[5, 5],
                size=2
            ).encode(
                x='x:Q'
            )
            
            chart_r2 = hist_r2 + mean_line_r2 + median_line_r2
        else:
            chart_r2 = alt.Chart(pd.DataFrame({'x': [0]})).mark_text(
                text='No data available'
            ).properties(width=400, height=300)
        
        # Combine charts horizontally
        combined_chart = alt.hconcat(chart_r1, chart_r2).configure_axis(
            labelFontSize=11,
            titleFontSize=12
        ).configure_title(
            fontSize=14
        )
        
        # Save as PNG
        output_path_png = os.path.join(output_dir, f'{gw_name}_read_length_histogram.png')
        try:
            combined_chart.save(output_path_png, scale_factor=2.0)
            print(f"  Saved: {output_path_png}")
        except Exception as e:
            print(f"  Warning: Could not save PNG. Trying HTML fallback...")
            print(f"         Error: {e}")
            # Fallback to HTML if PNG export fails
            output_path_html = os.path.join(output_dir, f'{gw_name}_read_length_histogram.html')
            combined_chart.save(output_path_html)
            print(f"  Saved: {output_path_html}")


def create_read_distribution_barplot(experiment_name, results, output_dir):
    """
    Create bar plot showing distribution of reads across populations and categories.
    
    Parameters:
    -----------
    experiment_name : str
        Name of the experiment
    results : dict
        Analysis results from analyze_demultiplexing_results
    output_dir : str
        Directory to save plots
    """
    print("\nGenerating read distribution bar plots...")
    
    for gw_name, file_results in results.items():
        print(f"  Creating plot for {gw_name}...")
        
        # Prepare data
        categories = []
        counts = []
        colors = []
        category_types = []
        
        # Add populations
        for pop_name in sorted(file_results['populations'].keys()):
            categories.append(pop_name)
            counts.append(file_results['populations'][pop_name])
            colors.append('steelblue')
            category_types.append('Population')
        
        # Add short reads
        categories.append('Short Reads')
        counts.append(file_results['short_reads'])
        colors.append('orange')
        category_types.append('QC')
        
        # Add unmatched reads
        categories.append('Unmatched Reads')
        counts.append(file_results['unmatched_reads'])
        colors.append('crimson')
        category_types.append('QC')
        
        # Calculate percentages
        total = file_results['input_total']
        percentages = [(count / total * 100) if total > 0 else 0 for count in counts]
        
        # Create DataFrame
        df = pd.DataFrame({
            'Category': categories,
            'Count': counts,
            'Color': colors,
            'Type': category_types,
            'Percentage': percentages,
            'Label': [f"{count:,}" for count in counts],
            'PctLabel': [f"{pct:.1f}%" for pct in percentages]
        })
        
        # Create bar chart
        bars = alt.Chart(df).mark_bar(
            opacity=0.7,
            stroke='black',
            strokeWidth=1
        ).encode(
            x=alt.X('Category:N', title='Category', sort=categories),
            y=alt.Y('Count:Q', title='Number of Reads'),
            color=alt.Color('Color:N', scale=None),
            tooltip=[
                alt.Tooltip('Category:N', title='Category'),
                alt.Tooltip('Count:Q', title='Count', format=','),
                alt.Tooltip('Percentage:Q', title='Percentage', format='.1f')
            ]
        )
        
        # Add count labels on top of bars
        text_top = alt.Chart(df).mark_text(
            align='center',
            baseline='bottom',
            dy=-5,
            fontSize=11
        ).encode(
            x=alt.X('Category:N', sort=categories),
            y=alt.Y('Count:Q'),
            text='Label:N'
        )
        
        # Add percentage labels in the middle of bars
        text_middle = alt.Chart(df).mark_text(
            align='center',
            baseline='middle',
            color='white',
            fontSize=10,
            fontWeight='bold'
        ).encode(
            x=alt.X('Category:N', sort=categories),
            y=alt.Y('Count:Q', stack='zero'),
            text='PctLabel:N'
        ).transform_filter(
            alt.datum.Count > 0
        )
        
        # Combine layers
        chart = (bars + text_top + text_middle).properties(
            width=600,
            height=400,
            title={
                'text': f'Read Distribution for {gw_name}',
                'subtitle': f'Total Input: {file_results["input_total"]:,} reads'
            }
        ).configure_axis(
            labelFontSize=11,
            titleFontSize=12,
            labelAngle=-45
        ).configure_title(
            fontSize=14
        )
        
        # Save as PNG
        output_path_png = os.path.join(output_dir, f'{gw_name}_read_distribution.png')
        try:
            chart.save(output_path_png, scale_factor=2.0)
            print(f"  Saved: {output_path_png}")
        except Exception as e:
            print(f"  Warning: Could not save PNG. Trying HTML fallback...")
            print(f"         Error: {e}")
            # Fallback to HTML if PNG export fails
            output_path_html = os.path.join(output_dir, f'{gw_name}_read_distribution.html')
            chart.save(output_path_html)
            print(f"  Saved: {output_path_html}")


def generate_summary_report(experiment_name, results, output_dir):
    """
    Generate a text summary report of the analysis.
    
    Parameters:
    -----------
    experiment_name : str
        Name of the experiment
    results : dict
        Analysis results from analyze_demultiplexing_results
    output_dir : str
        Directory to save the report
    """
    output_path = os.path.join(output_dir, 'index_quality_summary.txt')
    
    with open(output_path, 'w') as f:
        f.write(f"Index Quality Analysis Report\n")
        f.write(f"Experiment: {experiment_name}\n")
        f.write(f"=" * 80 + "\n\n")
        
        for gw_name, file_results in results.items():
            f.write(f"\n{gw_name}\n")
            f.write("-" * 80 + "\n")
            f.write(f"Total input reads: {file_results['input_total']:,}\n\n")
            
            f.write("Population distribution:\n")
            for pop_name in sorted(file_results['populations'].keys()):
                count = file_results['populations'][pop_name]
                pct = (count / file_results['input_total'] * 100) if file_results['input_total'] > 0 else 0
                f.write(f"  {pop_name}: {count:,} reads ({pct:.1f}%)\n")
            
            f.write("\nQuality control reads:\n")
            short_pct = (file_results['short_reads'] / file_results['input_total'] * 100) if file_results['input_total'] > 0 else 0
            f.write(f"  Short reads (<150 bp): {file_results['short_reads']:,} reads ({short_pct:.1f}%)\n")
            
            unmatched_pct = (file_results['unmatched_reads'] / file_results['input_total'] * 100) if file_results['input_total'] > 0 else 0
            f.write(f"  Unmatched reads: {file_results['unmatched_reads']:,} reads ({unmatched_pct:.1f}%)\n")
            
            total_distributed = (sum(file_results['populations'].values()) + 
                                file_results['short_reads'] + 
                                file_results['unmatched_reads'])
            recovery_rate = (total_distributed / file_results['input_total'] * 100) if file_results['input_total'] > 0 else 0
            f.write(f"\nTotal distributed: {total_distributed:,} reads\n")
            f.write(f"Recovery rate: {recovery_rate:.1f}%\n")
            f.write("\n")
    
    print(f"  Saved summary report: {output_path}")


def main():
    """Main function to run the analysis."""
    parser = argparse.ArgumentParser(
        description='Check quality of indexing/NGS data after demultiplexing'
    )
    parser.add_argument(
        'experiment_name',
        help='Name of the experiment (e.g., "example")'
    )
    
    args = parser.parse_args()
    experiment_name = args.experiment_name
    
    print(f"\n{'='*80}")
    print(f"Index Quality Analysis: {experiment_name}")
    print(f"{'='*80}")
    
    # Create output directory
    output_dir = os.path.join("Outputs", experiment_name, "index_quality")
    os.makedirs(output_dir, exist_ok=True)
    print(f"\nOutput directory: {output_dir}")
    
    # Analyze demultiplexing results
    results, input_files = analyze_demultiplexing_results(experiment_name)
    
    # Create visualizations
    create_read_length_histogram(experiment_name, input_files, output_dir)
    create_read_distribution_barplot(experiment_name, results, output_dir)
    
    # Generate summary report
    generate_summary_report(experiment_name, results, output_dir)
    
    print(f"\n{'='*80}")
    print("Analysis complete!")
    print(f"{'='*80}\n")


if __name__ == "__main__":
    main()
