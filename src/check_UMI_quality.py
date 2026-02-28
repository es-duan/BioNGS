#!/usr/bin/env python3
"""
Script 5: Check quality of UMI data

This script analyzes the quality of UMI dictionaries created by Script 4 (demultiplex_UMI.py).
It generates visualizations comparing the number of UMIs per population and the distribution 
of reads per UMI.

Input: UMI dictionary pickle files from results/{experiment}/demultiplexing/P{Population}/
Output: Bar plots saved to results/{experiment}/UMI_quality/

Dependencies: Script 4 (demultiplex_UMI.py) must be run first
Usage: python check_UMI_quality.py <experiment_name>
"""

import os
import pickle
import glob
import argparse
from collections import Counter
import altair as alt
import pandas as pd


def find_umi_libraries(experiment_name, umi_dir=None, output_base='results'):
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
    if umi_dir is None:
        search_dir = os.path.join(output_base, experiment_name, 'demultiplexing')
    else:
        search_dir = umi_dir

    pattern = os.path.join(search_dir, 'P*', '*_UMI_dict.pkl')
    libraries = {}
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
    """
    reads_per_umi = []
    
    for umi_pair, sequences in umi_dict.items():
        num_reads = len(sequences['R1'])
        reads_per_umi.append(num_reads)
    
    return {
        'num_umis': len(umi_dict),
        'reads_per_umi': reads_per_umi,
        'total_reads': sum(reads_per_umi) if reads_per_umi else 0
    }


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
    Create a box plot comparing the distribution of reads per UMI across populations.
    
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
    
    # Create box plot
    chart = alt.Chart(df).mark_boxplot(
        color='lightblue',
        opacity=0.7,
        size=50
    ).encode(
        x=alt.X('Population:N', title='Population', sort=populations),
        y=alt.Y('Reads_per_UMI:Q', title='Number of Reads per UMI'),
        tooltip=[
            alt.Tooltip('Population:N', title='Population'),
            alt.Tooltip('min(Reads_per_UMI):Q', title='Min', format='.0f'),
            alt.Tooltip('q1(Reads_per_UMI):Q', title='Q1', format='.0f'),
            alt.Tooltip('median(Reads_per_UMI):Q', title='Median', format='.0f'),
            alt.Tooltip('q3(Reads_per_UMI):Q', title='Q3', format='.0f'),
            alt.Tooltip('max(Reads_per_UMI):Q', title='Max', format='.0f')
        ]
    ).properties(
        width=600,
        height=400,
        title='Distribution of Reads per UMI Pair by Population'
    ).configure_axis(
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


def check_umi_quality(experiment_name, umi_dir=None, outdir=None):
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
        libraries = find_umi_libraries(experiment_name, umi_dir=umi_dir)
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
        except Exception as e:
            print(f"    Error loading {population}: {e}")
            continue
    
    if not libraries_data:
        print("Error: No libraries could be loaded")
        return False
    
    # Create output directory
    if outdir is None:
        output_dir = os.path.join('results', experiment_name, 'UMI_quality')
    else:
        output_dir = outdir
    os.makedirs(output_dir, exist_ok=True)
    print(f"\nOutput directory: {output_dir}")
    
    # Generate visualizations
    print("\nGenerating visualizations...")
    create_umi_count_plot(libraries_data, output_dir)
    create_reads_per_umi_plot(libraries_data, output_dir)
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
    parser.add_argument(
        '--demux_dir',
        default=None,
        help='Path to demultiplexing directory (e.g., results/example/demultiplexing/run_150)'
    )

    parser.add_argument(
        '--outdir',
        default=None,
        help='Output directory for UMI quality plots (e.g., results/example/UMI_quality/run_150)'
    )
    args = parser.parse_args()

    success = check_umi_quality(
        args.experiment_name,
        umi_dir=args.demux_dir,
        outdir=args.outdir
    )

    if success:
        print("\n" + "=" * 60)
        print("UMI quality check complete!")
        print("=" * 60)
    else:
        print("\nUMI quality check failed")


if __name__ == '__main__':
    main()
