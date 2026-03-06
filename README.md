# BioNGS

BioNGS is a structured, stage-based Next-Generation Sequencing (NGS) processing pipeline designed for bacterial mutation rate analysis with UMI-aware amplicon sequencing.

The pipeline emphasizes:

- reproducibility

- deterministic stage execution

- transparent quality control

- machine-readable metrics

- interactive preprocessing

BioNGS is particularly suited for **amplicon sequencing libraries containing inline indexes and UMIs**.

```
Library example
INDEX - PRIMER2 - UMI - PRIMER1 - INSERT - PRIMER1 - UMI - PRIMER2 - INDEX
```

---
## Installation

From the project root directory:
```
conda env create -f environment.yml
conda activate biongs

pip install -e .
````
This registers the biongs command-line tool.
Required external tools:

- FastQC

- fastp

- Bowtie2

- Samtools

# Running the Pipeline

After installing the package, the BioNGS pipeline can be executed using the biongs command-line interface.

## Run the pipeline

Example command:
```
biongs run \
  --exp example \
  --run_tag run1 \
  --r1 input_data/example/R1.fastq.gz \
  --r2 input_data/example/R2.fastq.gz \
  --multiplex_csv input_data/example/example_multiplexing_info.csv \
  --umi_primers_csv input_data/example/example_UMI_primers.csv
```
## Parameter description

| Parameter | Description |
|-----------|-------------|
| `--exp` | Experiment name. Results will be written to `results/<exp>/`. |
| `--run_tag` | Run identifier for this execution (e.g. `run1`). |
| `--r1` | Input FASTQ file for read 1. |
| `--r2` | Input FASTQ file for read 2. |
| `--multiplex_csv` | CSV file defining population index mapping. |
| `--umi_primers_csv` | CSV file defining UMI primer structures. |

## Pipeline Workflow

BioNGS executes a deterministic multi-stage workflow:

```
Raw FASTQ
   │
   ▼
Stage0   Input validation
   │
   ▼
Stage1   Raw QC (FastQC + overview)
   │
   ▼
Stage2   Demultiplex  (index retained)
   │
   ▼
Stage2.5 Post-demux QC
   │
   ▼
Stage3   UMI extraction and Flanking sequences stripping
   │
   ▼
Stage3.5 UMI quality summary
   │
   ▼
Stage4   QC after structure stripping
   │
   ▼
Stage5   Preprocessing (tail trimming + length filtering)
   │
   ▼
Stage6   Post-trim QC
   │
   ▼
Stage7   UMI collapsing / consensus generation
```

Each stage:

- produces structured outputs

- writes machine-readable metrics

- updates a global run_manifest.json

- can be safely resumed or rerun

Pipeline behavior is fully deterministic and stage-isolated.

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
│
├── environment.yml
└── README.md
```
Detailed script-level documentation is available:
See [src/README.md](src/README.md) for detailed usage instructions for each script.

---
## Stage Descriptions
### Stage 0 — Input Validation

Validates:

- FASTQ file integrity

- paired-end consistency

- gzip readability

- multiplex CSV structure

- uniqueness of index pairs

If any validation fails, the pipeline stops immediately.


### Stage 1 — Raw Quality Control

Performs two independent QC layers:

### QC Details

Runs FastQC and extracts report summaries.

### QC Overview

Computes statistics directly from FASTQ:

- base-level Q30 rate

- read length distribution

- short read rate (<150bp)

Decision logic:
```
short_rate > 30% → warning
Q30 < 70% → warning
short_rate > 30% AND Q30 < 70% → abnormal (user prompted to stop)
```
This prevents proceeding with severely degraded sequencing data.

### Stage 2 — Demultiplexing

Reads are split by inline index sequences.

Features:

- population-specific FASTQ generation

- unmatched read tracking

- strict FASTQ format validation

- automatic abnormality detection

Metrics recorded:
```
total_pairs
unmatched_pairs
unmatched_rate
per_population_pairs
index_lengths
```
Decision logic:
```
unmatched_rate > 50% → abnormal → user prompt
unmatched_rate > 35% → strong warning
unmatched_rate > 20% → warning
```
This helps detect incorrect index configuration or orientation errors.

### Stage 2.5 — Demultiplex QC

Performs QC for each population.

Outputs:
```
qc_overview/
   per_population/

qc_details/
   fastqc_reports/
```
Also generates:
```
demux_summary.txt
demux_distribution.png
```
showing population read distribution.

### Stage 3 — UMI Extraction

Extracts UMIs embedded in primer structures.

Primer pattern example:
```
R1: CCGCGTGATTACGAGTCG-[10N UMI]-GCAGCAGTGAAAGAGTTCTTCG
R2: GGGTTAGCAAGTGGCAGCCT-[10N UMI]-GGAAGCCGTATTCGTTAGTCTG
```
UMI extraction process:

- detect primer prefix

- extract UMI sequence

- verify suffix anchor

- trim structural sequence

- retain insert sequence

Failure reasons are recorded:
```
missing_prefix
missing_anchor
too_short
invalid_umi_chars
qual_len_mismatch
```
Stage 3 records per-population UMI extraction statistics, including unmatched rate, pair-level failure reasons, and mate-specific failure reasons (R1/R2), and writes a umi_extraction_summary.txt report for each population.

### Stage 3.5 — UMI Family Statistics

Analyzes UMI family structure.

Outputs:
```
UMI_quality_summary.txt
```
Metrics include:

- number of unique UMI families

- mean reads per family

- family size distribution

- fraction of families ≥2 reads

- fraction of families ≥3 reads

These metrics are useful for evaluating sequencing depth and PCR duplication levels.

### Stage 4 — QC after Structure Removal

After removing:

- index

- primers

- UMI sequences

QC is recomputed on clean insert sequences.

Outputs:
```
qc_overview/02_clean_after_umi/
qc_details/02_clean_after_umi/
```
Statistics include:

- read counts

- Q30 rates

- read length distribution.

### Stage 5 — Interactive Preprocessing

Stage5 performs trimming and filtering using fastp.

This stage contains the main interactive decision system.

#### Step 5.1 Candidate Q Evaluation

Candidate trimming thresholds are evaluated:
```
Q30
Q29
Q28
Q25
Q20
```
Evaluation is performed on a sampled subset of reads.

Metrics reported:
```
Q30_after
trimmed_bases%
mean_read_length
```
The user selects a trimming threshold.

#### Step 5.2 Length Filtering Preview

After trimming, the pipeline previews the impact of length filtering.

Preview thresholds:
```
180
170
160
150
140
130
```
Short-read definition:
```
min(read_length_R1, read_length_R2) < min_len
```
Warning system:
```
short_rate > 30% → warning
short_rate > 40% → strong warning
short_rate > 50% → abnormal (user prompt)
```
Users may return to Step 5.1 and rerun trimming with different parameters.

#### Step 5.3 Final Filtering

After user confirmation, final filtering is executed and trimmed FASTQ files are produced.

### Stage 6 — Post-Trim QC

Final QC after preprocessing.

Outputs:
```
qc_overview/01_post_trim/
qc_details/01_post_trim/
```
Reports include:

- per-population read length distributions

- overall aggregated statistics

- FastQC reports.

## Output Directory Structure

Example result structure:
```
results/
   <exp>/
      demultiplexing/
      umi_extracted/
      trimmed/
      qc_overview/
      qc_details/
      metrics/
      run_manifest.json
```
run_manifest.json records the execution status of every pipeline stage.


## Additional options

### List pipeline stages:
```
biongs run --list_stages
```
This command prints the names of all pipeline stages and exits without running the pipeline.

This is useful if you want to know the exact stage names required by the --rerun_from option.

### Rerun the pipeline starting from a specific stage:
```
biongs run --rerun_from stage5_preprocess ...
```
This option forces the pipeline to rerun from the specified stage, even if previous results already exist.

For example:
```
biongs run \
  --exp example \
  --run_tag run1 \
  --r1 R1.fastq.gz \
  --r2 R2.fastq.gz \
  --multiplex_csv multiplex.csv \
  --umi_primers_csv umi_primers.csv \
  --rerun_from stage5_preprocess
```

In this case:

- stages before stage5_preprocess will be skipped (if already completed)

- stage5_preprocess and all following stages will be executed again

This is useful if you want to adjust trimming parameters or rerun preprocessing without repeating the earlier stages.

### Disable automatic resume:
```
biongs run --no_resume ...
```
By default, the pipeline automatically detects completed stages and skips them to save time.

The --no_resume option disables this behavior and forces the pipeline to run every stage again from the beginning, even if previous results exist.

This is helpful when you want a completely fresh run.


## Design Philosophy

BioNGS prioritizes:

- reproducible sequencing analysis

- explicit user decisions

- transparent QC thresholds

- machine-readable pipeline metrics

- strict stage isolation

The pipeline is designed to help researchers detect sequencing issues early while maintaining reproducible preprocessing decisions.