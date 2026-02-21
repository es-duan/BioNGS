#!/usr/bin/env python3
"""
Script 5: Check quality of UMI data

This script analyzes the quality of UMI dictionaries created by Script 4 (demultiplex_UMI.py).
It generates visualizations comparing the number of UMIs per population and the distribution 
of reads per UMI.

Input: UMI dictionary pickle files from Outputs/{experiment}/demultiplexing/P{Population}/
Output: Bar plots saved to Outputs/{experiment}/UMI_quality/

Dependencies: Script 4 (demultiplex_UMI.py) must be run first
Usage: python check_UMI_quality.py <experiment_name>
"""

import os
import pickle
import glob
import argparse
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np


def find_umi_libraries(experiment_name, output_base='Outputs'):
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
    
    plt.figure(figsize=(10, 6))
    bars = plt.bar(populations, umi_counts, color='steelblue', edgecolor='black', alpha=0.7)
    
    # Add value labels on top of bars
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}',
                ha='center', va='bottom', fontsize=10)
    
    plt.xlabel('Population', fontsize=12, fontweight='bold')
    plt.ylabel('Number of Unique UMI Pairs', fontsize=12, fontweight='bold')
    plt.title('Number of Unique UMI Pairs per Population', fontsize=14, fontweight='bold')
    plt.grid(axis='y', alpha=0.3, linestyle='--')
    plt.tight_layout()
    
    output_path = os.path.join(output_dir, 'UMI_count_per_population.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_path}")
    plt.close()


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
    reads_distributions = [libraries_data[pop]['reads_per_umi'] for pop in populations]
    
    plt.figure(figsize=(12, 6))
    bp = plt.boxplot(reads_distributions, labels=populations, patch_artist=True)
    
    # Color the boxes
    for patch in bp['boxes']:
        patch.set_facecolor('lightblue')
        patch.set_alpha(0.7)
    
    plt.xlabel('Population', fontsize=12, fontweight='bold')
    plt.ylabel('Number of Reads per UMI', fontsize=12, fontweight='bold')
    plt.title('Distribution of Reads per UMI Pair by Population', fontsize=14, fontweight='bold')
    plt.grid(axis='y', alpha=0.3, linestyle='--')
    plt.tight_layout()
    
    output_path = os.path.join(output_dir, 'reads_per_UMI_distribution.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_path}")
    plt.close()


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
                f.write(f"  Reads per UMI - Min: {min(reads_per_umi)}, ")
                f.write(f"Max: {max(reads_per_umi)}, ")
                f.write(f"Mean: {np.mean(reads_per_umi):.1f}, ")
                f.write(f"Median: {np.median(reads_per_umi):.1f}\n")
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
        except Exception as e:
            print(f"    Error loading {population}: {e}")
            continue
    
    if not libraries_data:
        print("Error: No libraries could be loaded")
        return False
    
    # Create output directory
    output_dir = os.path.join('Outputs', experiment_name, 'UMI_quality')
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
        help='Name of the experiment (e.g., "example"). Script will analyze UMI libraries from Outputs/{experiment_name}/demultiplexing/'
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
