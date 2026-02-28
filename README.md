# BioNGS

BioNGS is a structured Next-Generation Sequencing (NGS) processing pipeline developed for bacterial mutation rate analysis.

Originally built for CHEM E 546, it has evolved into a modular and reproducible demultiplexing + quality control framework supporting adaptive filtering and dual-layer QC reporting.

---

## Overview

This pipeline provides:

- DNA index-based demultiplexing
- Configurable read-length filtering (default: 150 bp)
- Automated discard-rate warning system
- Reproducible multi-threshold runs (run_150, run_130, etc.)
- Dual-level quality control:
  - QC Overview (lab-friendly summary)
  - QC Details (professional FastQC-style report)

The design enforces strict separation of sequencing outputs and QC outputs to maintain clarity and reproducibility.

---

## Project Structure
```
BioNGS/
│
├── src/ # Source code
│    ├── qc/ # QC driver & fancy reports
│    │    ├── entry_qc.py
│    │    ├── qc_driver.py
│    │    ├── oi/ 
│    │    │    ├── manifest.py
│    │    │    └── paths.py
│    │    ├── stages/
│    │    │     ├── stage_after_demux.py
│    │    │     └── stage_raw.py
│    │    ├── entry_qc.py
│    │    └── qc_driver.py
│    │
│    ├── check_index_quality.py
│    ├── check_UMI_quality.py
│    ├── demultiplex_folders.py
│    ├── demultiplex_index.py
│    ├── demultiplex_UMI.py
│    ├── fastqc_unpacker.py
│    └── README.md
│ 
│
├── input_data/ # Raw sequencing data
├── results/ # Generated outputs
│    ├── demultiplexing
│    ├── manifests
│    ├── qc_details
│    └── qc_overview
│ 
├── environment.yml # Conda environment
└── README.md # Project overview
```
Detailed script-level documentation is available:
See [src/README.md](src/README.md) for detailed usage instructions for each script.


---

## Setup

### Recommended: Conda Installation

```bash
conda env create -f environment.yml
conda activate biongs
```

This will install all required dependencies:

- Core

  - Python 3.11

  - BioPython

  - NumPy

  - Pandas

  - Matplotlib

- QC & Visualization

  - FastQC

  - Altair

  - vl-convert-python

- Alignment / Bioinformatics Tools

  - Bowtie2

  - Samtools


### Alternative: Manual Setup
```bash
python -m venv .venv
source .venv/bin/activate

pip install -e .
```
Ensure bowtie2 is installed system-wide:
- macOS: brew install bowtie2

- Linux: sudo apt-get install bowtie2

## Quick Start

Run baseline demultiplexing with default 150 bp filtering:
```bash
python -m src.qc.qc_driver <experiment_name> --gw_name <GW_NAME>
```
If discard rate exceeds threshold (default 30%), the pipeline will:

- Print a warning

- Suggest re-running with a different minimum length

- Optionally prompt for a new threshold

### Each threshold produces an independent run directory:
```
results/<experiment>/
  demultiplexing/run_150/
  demultiplexing/run_130/
  qc_overview/run_150/
  qc_overview/run_130/
  qc_details/...
```

## Quality Control Philosophy
### QC Overview
Designed for rapid laboratory assessment:

- Read distribution

- Length histogram

- Discard metrics

- Machine-readable JSON

### QC Details

Professional inspection via FastQC unpacked reports:

- Raw data

- Post-demultiplex data

This separation ensures:

- Clean sequencing data directories

- Clear QC reporting structure

- Reproducible threshold comparison

## Design Principles

- Reproducibility over convenience

- No overwriting between thresholds

- Clear separation of data and QC

- Machine-readable metrics for downstream automation

- Lab-standard defaults with adaptive flexibility

## Documentation
For detailed usage instructions of each script, see: [src/README.md](src/README.md)
