# BioNGS

BioNGS is a fully structured, stage-based Next-Generation Sequencing (NGS) processing pipeline designed for bacterial mutation rate analysis and UMI-aware processing.

The pipeline enforces strict reproducibility, modular stage separation, machine-readable metrics, and interactive quality-aware preprocessing.

---

## Pipeline Overview

BioNGS follows a deterministic stage-based workflow:

```
Stage 0   Input validation
Stage 1   Raw QC
Stage 2   Demultiplex (inline index)
Stage 2.5 Post-demux QC
Stage 3   UMI extraction (strict + optional relaxed rerun)
Stage 3.5 UMI quality summary
Stage 4   QC after structure stripping
Stage 5   Preprocessing (interactive tail-trim + length filtering)
Stage 6   Post-preprocess QC
Stage 7   UMI collapsing / consensus generation
```

Each stage:

- Writes structured outputs
- Generates machine-readable metrics
- Updates a global run_manifest.json
- Never silently overwrites without explicit confirmation

---

## Project Structure
```
BioNGS/
│
├── src/
│   ├── run_pipeline.py
│   ├── stage0_validate.py
│   ├── stage1_raw_qc.py
│   ├── stage2_demux.py
│   ├── stage3_umi_extract.py
│   ├── stage3_5_umi_quality.py
│   ├── stage5_preprocess.py
│   ├── stage6_post_trim_qc.py
│   └── ...
│
├── input_data/
├── results/
│   ├── <exp>/
│   │    ├── demultiplexing/
│   │    ├── umi_extracted/
│   │    ├── trimmed/
│   │    ├── qc_overview/
│   │    ├── qc_details/
│   │    ├── trim_reports/
│   │    └── run_manifest.json
│
├── environment.yml
└── README.md
```
Detailed script-level documentation is available:
See [src/README.md](src/README.md) for detailed usage instructions for each script.

---

## Key Features
1. Structured Stage Control
  - Stage-by-stage execution
  - Automatic status tracking
  - Manifest-based reproducibility
  - Safe abort logic


2. Inline Index Demultiplexing
  - Index-based population splitting
  - Immediate error detection for broken read pairs
  - Unmatched tracking
  - Terminal warnings:
    - 20% warning
    - 35% strong warning
    - 50% abnormal (optional stop)
3. UMI Extraction with Adaptive Rerun
- UMI extraction supports:
    - Strict mode (mismatch=0, shift=0)
    - Optional relaxed rerun profiles:
      - mismatch=1, shift=0
      - mismatch=1, shift=2

- Pipeline reports transparent fail reasons:
  - missing_prefix_left
  - missing_anchor_right
  - too_short
  - invalid_umi_chars
  - other_parse_fail

User may rerun with relaxed profile.
Rerun overwrites Stage3 outputs safely to preserve downstream compatibility.

4. Interactive Preprocessing (Stage 5)
Stage 5 includes:
- 5.1 Tail Trimming (fastp)
  - Sampled Q evaluation (Q30/Q29/Q28/Q25/Q20)
  - User selects Q threshold
  - Automatic warnings:
    - Q30 < 85% → WARNING
    - Q30 < 80% → STRONG WARNING
    - Q30 < 70% → ABNORMAL (optional stop)

5.2 Length Filtering Preview

Before filtering, pipeline displays:
  - Short-rate table at multiple thresholds (e.g. 180/170/160/150/140/130)
  - air-level short definition:
    - A pair is short if either mate < min_len

  - Warning levels:
    - 30% WARNING
    - 40% STRONG WARNING
    - 50% ABNORMAL (optional stop)

Old preview versions are preserved (v1, v2, ...).


5. Dual-Layer QC System
- QC Overview
  - Lab-friendly:
    - Length distribution
    - Summary statistics
    - Short-rate previews

- QC Details
Professional FastQC-style unpacked reports:
  - Raw
  - Post-demux
  - Post-trim

6. Machine-Readable Metrics
  Every stage writes:
  ```
  results/<exp>/metrics/
  ```

  - Includes:
      - metrics_demux.json
      - metrics_trim.json
      - umi_quality_summary.txt
      - run_manifest.json

  - metrics_trim.json records:
      - q_threshold
      - min_len
      - Q30 before/after
      - short_rate
      - dropped reasons
      - user decisions

## Installation

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

  - Fastp

  - Altair

  - vl-convert-python

- Alignment / Bioinformatics Tools

  - Bowtie2

  - Samtools


## Quick Start

From project root:

```bash
conda activate biongs

python -m src.run_pipeline \
  --exp example \
  --run_tag run1 \
  --r1 input_data/example/R1.fastq.gz \
  --r2 input_data/example/R2.fastq.gz \
  --multiplex_csv input_data/example/example_multiplexing_info.csv
```

### Parameter Meaning

- --exp : experiment name (creates results/<exp>/)

- --run_tag : run version label

- --r1 / --r2 : raw FASTQ files

- --multiplex_csv : index → population mapping table


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

## Design Philosophy

BioNGS prioritizes:

  - Engineering reproducibility

  - Explicit user decisions

  - Transparent failure reporting

  - No silent data loss

  - Strict stage isolation

  - Overwrite safety when rerunning structural stages
  
## Documentation
For detailed usage instructions of each script, see: [src/README.md](src/README.md)
