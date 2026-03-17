#!/usr/bin/env python3
"""
Step 3a: Check quality of UMI data

This step analyzes the quality of UMI dictionaries created by Step 3 (step3_demultiplex_UMI.py).
It generates visualizations comparing the number of UMIs per population and the distribution 
of reads per UMI.

Input: UMI dictionary pickle files from results/{experiment}/demultiplexing/P{Population}/
Output: Bar plots saved to results/{experiment}/UMI_quality/

Dependencies: Step 3 (step3_demultiplex_UMI.py) must be run first
Usage: python step3a_check_UMI_quality.py <experiment_name>
"""

import os
import pickle
import glob
import argparse
from collections import Counter
import altair as alt
import pandas as pd
from Bio import SeqIO


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
        UMI dictionary: {(forward_UMI, reverse_UMI): {'R1': [...], 'R2': [...]}}
    """
    with open(library_path, 'rb') as f:
        umi_dict = pickle.load(f)
    return umi_dict


def analyze_umi_dict(umi_dict):
    """
    Analyze a UMI dictionary and extract quality metrics.
    
    Parameters:
    -----------
    umi_dict : dict
        UMI dictionary structure
        
    Returns:
    --------
    dict
        Dictionary containing:
        - num_umis: Total number of unique UMI pairs
        - reads_per_umi: List of read counts for each UMI
        - total_reads: Total number of reads in dictionary
        - min_reads_per_umi: Minimum reads among UMI pairs
        - max_reads_per_umi: Maximum reads among UMI pairs
    """
    reads_per_umi = []
    
    for umi_pair, sequences in umi_dict.items():
        num_reads = len(sequences['R1'])
        reads_per_umi.append(num_reads)
    
    return {
        'num_umis': len(umi_dict),
        'reads_per_umi': reads_per_umi,
        'total_reads': sum(reads_per_umi) if reads_per_umi else 0,
        'min_reads_per_umi': min(reads_per_umi) if reads_per_umi else 0,
        'max_reads_per_umi': max(reads_per_umi) if reads_per_umi else 0,
    }


def count_unmatched_reads(population_folder, population):
    """
    Count the number of unmatched UMI reads for a population.
    
    Parameters:
    -----------
    population_folder : str
        Path to the population folder (e.g., demultiplexing/P1/)
    population : str
        Population name (e.g., "P1")
        
    Returns:
    --------
    int
        Number of unmatched reads (reads without UMI match in either R1 or R2)
    """
    unmatched_r1_fastq = os.path.join(population_folder, f"{population}_unmatched_UMI_R1.fastq")
    
    if not os.path.exists(unmatched_r1_fastq):
        return 0
    
    count = 0
    try:
        with open(unmatched_r1_fastq, 'r') as f:
            for record in SeqIO.parse(f, 'fastq'):
                count += 1
    except Exception as e:
        print(f"    Warning: Could not count unmatched reads for {population}: {e}")
        return 0
    
    return count


def create_umi_count_plot(libraries_data, output_dir):
    """
    Create a bar plot comparing the number of UMIs per population.
    
    Parameters:
    -----------
    libraries_data : dict
        Dictionary mapping population names to analysis results
    output_dir : str
        Directory to save the plot
    """
    populations = sorted(libraries_data.keys())
    umi_counts = [libraries_data[pop]['num_umis'] for pop in populations]
    
    # Create DataFrame
    df = pd.DataFrame({
        'Population': populations,
        'UMI_Count': umi_counts,
        'Label': [f'{count}' for count in umi_counts]
    })
    
    # Create bar chart
    bars = alt.Chart(df).mark_bar(
        opacity=0.7,
        color='steelblue',
        stroke='black',
        strokeWidth=1
    ).encode(
        x=alt.X('Population:N', title='Population', sort=populations),
        y=alt.Y('UMI_Count:Q', title='Number of Unique UMI Pairs'),
        tooltip=[
            alt.Tooltip('Population:N', title='Population'),
            alt.Tooltip('UMI_Count:Q', title='UMI Count', format=',')
        ]
    )
    
    # Add value labels on top of bars
    text = alt.Chart(df).mark_text(
        align='center',
        baseline='bottom',
        dy=-5,
        fontSize=11
    ).encode(
        x=alt.X('Population:N', sort=populations),
        y=alt.Y('UMI_Count:Q'),
        text='Label:N'
    )
    
    # Combine layers
    chart = (bars + text).properties(
        width=500,
        height=400,
        title='Number of Unique UMI Pairs per Population'
    ).configure_axis(
        labelFontSize=11,
        titleFontSize=12
    ).configure_title(
        fontSize=14
    )
    
    # Save as PNG
    output_path = os.path.join(output_dir, 'UMI_count_per_population.png')
    try:
        chart.save(output_path, scale_factor=2.0)
        print(f"  Saved: {output_path}")
    except Exception as e:
        print(f"  Warning: Could not save PNG. Error: {e}")


def create_reads_per_umi_plot(libraries_data, output_dir):
    """
    Create violin plots comparing the distribution of reads per UMI across populations.
    
    Parameters:
    -----------
    libraries_data : dict
        Dictionary mapping population names to analysis results
    output_dir : str
        Directory to save the plot
    """
    populations = sorted(libraries_data.keys())
    
    # Create DataFrame with all data points
    data_rows = []
    for pop in populations:
        reads_per_umi = libraries_data[pop]['reads_per_umi']
        for count in reads_per_umi:
            data_rows.append({'Population': pop, 'Reads_per_UMI': count})
    
    df = pd.DataFrame(data_rows)
    
    # Create violin plots by faceting one panel per population
    violin = alt.Chart(df).transform_density(
        'Reads_per_UMI',
        as_=['Reads_per_UMI', 'Density'],
        groupby=['Population']
    ).mark_area(
        orient='horizontal',
        opacity=0.7,
        color='steelblue'
    ).encode(
        y=alt.Y('Reads_per_UMI:Q', title='Number of Reads per UMI'),
        x=alt.X('Density:Q', stack='center', title=None, axis=None),
        tooltip=[
            alt.Tooltip('Population:N', title='Population'),
            alt.Tooltip('Reads_per_UMI:Q', title='Reads per UMI', format='.0f'),
            alt.Tooltip('Density:Q', title='Density', format='.4f')
        ]
    ).properties(
        width=120,
        height=400
    )

    chart = violin.facet(
        column=alt.Column('Population:N', title='Population', sort=populations)
    ).properties(
        title='Distribution of Reads per UMI Pair by Population (Violin Plot)'
    ).configure_axis(
        labelFontSize=11,
        titleFontSize=12
    ).configure_header(
        labelFontSize=11,
        titleFontSize=12
    ).configure_title(
        fontSize=14
    )
    
    # Save as PNG
    output_path = os.path.join(output_dir, 'reads_per_UMI_distribution.png')
    try:
        chart.save(output_path, scale_factor=2.0)
        print(f"  Saved: {output_path}")
    except Exception as e:
        print(f"  Warning: Could not save PNG. Error: {e}")


def create_umi_match_status_plot(libraries_data, output_dir, experiment_name):
    """
    Create a faceted bar plot comparing matched vs unmatched UMI reads per population.
    
    Parameters:
    -----------
    libraries_data : dict
        Dictionary mapping population names to analysis results
    output_dir : str
        Directory to save the plot
    experiment_name : str
        Name of the experiment
    """
    populations = sorted(libraries_data.keys())
    demultiplex_base = os.path.join('results', experiment_name, 'demultiplexing')
    
    # Collect data for all populations
    data_rows = []
    for pop in populations:
        # Get matched UMI count
        matched_count = libraries_data[pop]['num_umis']
        
        # Get unmatched read count
        population_folder = os.path.join(demultiplex_base, pop)
        unmatched_count = count_unmatched_reads(population_folder, pop)
        
        # Add rows for the bar chart
        data_rows.append({'Population': pop, 'Status': 'Matched UMIs', 'Count': matched_count})
        data_rows.append({'Population': pop, 'Status': 'Unmatched Reads', 'Count': unmatched_count})
    
    df = pd.DataFrame(data_rows)
    
    # Create bar chart with faceting by population
    bars = alt.Chart(df).mark_bar(
        opacity=0.7,
        stroke='black',
        strokeWidth=1
    ).encode(
        x=alt.X('Status:N', title='UMI Status', sort=['Matched UMIs', 'Unmatched Reads']),
        y=alt.Y('Count:Q', title='Count'),
        color=alt.Color('Status:N', scale=alt.Scale(scheme='set2')),
        tooltip=[
            alt.Tooltip('Population:N', title='Population'),
            alt.Tooltip('Status:N', title='Status'),
            alt.Tooltip('Count:Q', title='Count', format=',')
        ]
    ).properties(
        width=200,
        height=300
    )
    
    # Add value labels on top of bars
    text = alt.Chart(df).mark_text(
        align='center',
        baseline='bottom',
        dy=-5,
        fontSize=10
    ).encode(
        x=alt.X('Status:N', sort=['Matched UMIs', 'Unmatched Reads']),
        y=alt.Y('Count:Q'),
        text=alt.Text('Count:Q', format=',')
    )
    
    # Facet by population
    chart = (bars + text).facet(
        column=alt.Column('Population:N', title='Population', sort=populations)
    ).properties(
        title='UMI Match Status per Population'
    ).configure_axis(
        labelFontSize=11,
        titleFontSize=12
    ).configure_header(
        labelFontSize=11,
        titleFontSize=12
    ).configure_title(
        fontSize=14
    )
    
    # Save as PNG
    output_path = os.path.join(output_dir, 'UMI_match_status.png')
    try:
        chart.save(output_path, scale_factor=2.0)
        print(f"  Saved: {output_path}")
    except Exception as e:
        print(f"  Warning: Could not save PNG. Error: {e}")


def create_summary_statistics_plot(libraries_data, output_dir):
    """
    Create a summary statistics table as a text file.
    
    Parameters:
    -----------
    libraries_data : dict
        Dictionary mapping population names to analysis results
    output_dir : str
        Directory to save the statistics
    """
    populations = sorted(libraries_data.keys())
    
    output_path = os.path.join(output_dir, 'UMI_quality_summary.txt')
    with open(output_path, 'w') as f:
        f.write("UMI Quality Summary Statistics\n")
        f.write("=" * 80 + "\n\n")
        
        for pop in populations:
            data = libraries_data[pop]
            reads_per_umi = data['reads_per_umi']
            
            f.write(f"Population: {pop}\n")
            f.write(f"  Number of unique UMI pairs: {data['num_umis']}\n")
            f.write(f"  Total reads: {data['total_reads']}\n")
            
            if reads_per_umi:
                reads_series = pd.Series(reads_per_umi)
                f.write(f"  Reads per UMI - Min: {reads_series.min()}, ")
                f.write(f"Max: {reads_series.max()}, ")
                f.write(f"Mean: {reads_series.mean():.1f}, ")
                f.write(f"Median: {reads_series.median():.1f}\n")
            f.write("\n")
    
    print(f"  Saved: {output_path}")


def check_umi_quality(experiment_name):
    """
    Main function to check UMI quality for an experiment.
    
    Parameters:
    -----------
    experiment_name : str
        Name of the experiment
    """
    print("=" * 60)
    print(f"Checking UMI quality for experiment: {experiment_name}")
    print("=" * 60)
    
    # Find UMI libraries
    print("\nSearching for UMI libraries...")
    try:
        libraries = find_umi_libraries(experiment_name)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return False
    
    if not libraries:
        print("Error: No UMI libraries found")
        return False
    
    print(f"Found {len(libraries)} population(s)")
    
    # Load and analyze dictionaries
    print("\nAnalyzing UMI dictionaries...")
    libraries_data = {}
    
    for population, lib_path in libraries.items():
        print(f"  Loading {population}...")
        try:
            umi_dict = load_umi_dict(lib_path)
            libraries_data[population] = analyze_umi_dict(umi_dict)
            print(f"    UMI pairs: {libraries_data[population]['num_umis']}, "
                  f"Total reads: {libraries_data[population]['total_reads']}")
            print(f"    Reads per UMI - Min: {libraries_data[population]['min_reads_per_umi']}, "
                  f"Max: {libraries_data[population]['max_reads_per_umi']}")
        except Exception as e:
            print(f"    Error loading {population}: {e}")
            continue
    
    if not libraries_data:
        print("Error: No libraries could be loaded")
        return False
    
    # Create output directory
    output_dir = os.path.join('results', experiment_name, 'UMI_quality')
    os.makedirs(output_dir, exist_ok=True)
    print(f"\nOutput directory: {output_dir}")
    
    # Generate visualizations
    print("\nGenerating visualizations...")
    create_umi_count_plot(libraries_data, output_dir)
    create_reads_per_umi_plot(libraries_data, output_dir)
    create_umi_match_status_plot(libraries_data, output_dir, experiment_name)
    create_summary_statistics_plot(libraries_data, output_dir)
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description='Check quality of UMI libraries for each population'
    )
    parser.add_argument(
        'experiment_name',
        help='Name of the experiment (e.g., "example"). Script will analyze UMI libraries from results/{experiment_name}/demultiplexing/'
    )
    
    args = parser.parse_args()
    
    success = check_umi_quality(args.experiment_name)
    
    if success:
        print("\n" + "=" * 60)
        print("UMI quality check complete!")
        print("=" * 60)
    else:
        print("\nUMI quality check failed")


if __name__ == '__main__':
    main()
