# BioNGS Pipeline Scripts

## Prerequisites

### Recommended: Use Conda with environment.yml

From the project root directory:
```bash
# Create the conda environment from the environment.yml file
conda env create -f environment.yml

# Activate the environment
conda activate biongs
```

This is the recommended approach as it installs all dependencies including bowtie2 (required for Script 6).

### Alternative: Manual Virtual Environment Setup

If you prefer to use pip:

```bash
# Create virtual environment
python -m venv .venv

# Activate the environment
# On macOS/Linux:
source .venv/bin/activate
# On Windows:
# .venv\Scripts\activate

# Install the package
pip install -e .
```

**Note:** You must also have bowtie2 installed on your system for Script 6:
- macOS: `brew install bowtie2`
- Linux: `sudo apt-get install bowtie2`

This document describes the internal structure and execution logic of the BioNGS pipeline.

The pipeline performs:

1. Demultiplexing by index
2. Minimum read-length filtering (default 150 bp)
3. Adaptive threshold warning
4. Dual-layer quality control reporting

---

#  Pipeline Overview

Execution entry point:
```bash
python -m src.qc.qc_driver <experiment> --gw_name <GW_NAME>
```
## Highpoint flow chart
```

RAW FASTQ
   │
   ├── QC Detail (RAW)  ─────────► qc_details/00_raw
   │
   ├── demultiplex
   │       │
   │       ├── drop reads <150 bp(deforlt)
   │       ├── index Demultiplexing
   │       │
   │       └── output FASTQ ───────► demultiplexing/run_150
   │
   ├── QC Overview ─────────────► qc_overview/run_150
   │
   └── QC Details (After demultiplexing) ────────► qc_details/02_after_demux/run_150

```
At this stage, discard_rate is evaluated.

### If Discard Rate is High
Instead of overwriting results, the pipeline supports branching runs:
```
               ┌────────────── run_150 ───────────────┐
               │   demultiplexing/                    │
RAW FASTQ ─────┼── qc_overview/run_150/               │
               │   qc_details/run_150/                │
               └──────────────────────────────────────┘
                               │
                               ▼
                     keep what pipeline did? -► yes -►keep the run_150 data
                               │
                               ▼
                               No
                               │
                               ▼                               
                     User selects new min_len (i.e. 130)
                               │
                               ▼
               ┌────────────── run_130 ───────────────┐
               │   demultiplexing/                    │
               │   qc_overview/run_130/               │
               │   qc_details/run_130/                │
               └──────────────────────────────────────┘
```

## Key Design Feature

Each threshold creates an independent versioned run:
- run_150
- run_130
- run_120


All stored side-by-side.

No data is overwritten.

---

#  Source Directory Structure
```
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
```

---

# Script Descriptions

---

## demultiplex_folders.py

Purpose:
- Create folder structure for demultiplexing
- Prepare population-specific directories

Input:
- experiment name

Output:
- results/<exp>/demultiplexing/run_<min_len>/

---

## demultiplex_index.py

Purpose:
- Match reads by R1/R2 index
- Apply minimum read-length filtering
- Classify reads into:
  - matched
  - short
  - unmatched

Key parameter: 
- --min_len

Default: 150

Outputs:
- Demultiplexed FASTQ files
- metrics.json (machine-readable summary)

---

## check_index_quality.py

Purpose:
- Generate simplified QC overview plots
- Provide quick lab-level inspection

Outputs:
- Read distribution bar plot
- Read length histogram
- index_quality_summary.txt
- metrics.json (passed through)

Stored under:
```
results/<exp>/qc_overview/run_<min_len>/
```

---

## qc/qc_driver.py

Main pipeline orchestrator.

Responsibilities:

- Execute demultiplexing scripts
- Manage directory structure
- Generate QC overview
- Generate QC details (FastQC unpacked)
- Detect high discard rate
- Support interactive rerun with new threshold

Important arguments:
- --min_len
- --warn_discard_rate
- --interactive

 
---

## qc/entry_qc.py

Purpose:
- Run FastQC
- Unpack FastQC results
- Generate structured QC report
- Export Altair figures as PNG

Used for:

- Raw data QC
- Post-demultiplex QC

---

# QC System

## QC Overview (qc_overview)

Designed for:
- Rapid evaluation
- Threshold comparison
- Lab discussion

Contains:
- Distribution plot
- Length histogram
- metrics.json

---

## QC Details (qc_details)

Designed for:
- Professional inspection
- Deep sequencing diagnostics
- Publication-level validation

Contains:
- Raw FastQC unpacked report
- Post-demux FastQC unpacked report

---

# Adaptive Filtering Logic

1. Run baseline min_len (default 150)
2. Compute discard_rate
3. If discard_rate > threshold (default 30%):
   - Print terminal warning
   - Suggest rerun
   - Optionally prompt user
4. New run_<min_len> generated independently

No runs are overwritten.

---

# metrics.json Structure

Example:

```json
{
  "overall": {
    "total_reads": 165273,
    "matched_reads": 101613,
    "short_reads": 56165,
    "unmatched_reads": 7495,
    "discard_rate": 0.3398
  }
}
```
Used for:

- Threshold warning logic

- Future automation

- Comparative analysis

# Running Individual Components

You may run modules independently:

## Demultiplex only
```
python src/demultiplex_index.py example --min_len 150
```
## QC Overview only
```
python src/demultiplex_index.py example --min_len 150
```
## QC Overview only
```
python src/qc/entry_qc.py --manifest manifest.csv --outdir outdir
```