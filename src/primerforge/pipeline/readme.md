# Stage 0 — Input Validation

Stage 0 performs strict validation of the input FASTQ files and multiplex configuration file before any downstream processing.

Implementation file: 
```
src/pipeline/stage0_validation.py
```

## Purpose

Stage 0 ensures that:
  - R1 and R2 FASTQ files are structurally valid.
  - Paired-end record counts match.
  - The multiplex CSV file is complete and logically consistent.
  - All index pairs are unique.

If any rule is violated, the pipeline stops immediately.

## Inputs

Required inputs:

  - Paired-end FASTQ files (R1 and R2; plain text or .gz)
  - Multiplex configuration CSV file

Example directory:
```
input_data/
├── sample_R1.fastq.gz
├── sample_R2.fastq.gz
└── multiplex_info.csv
```

## FASTQ Validation Rules

1. File existence

    Both R1 and R2 must exist.

2. Gzip integrity

    If a file ends with .gz, it must be decompressible.

3. FASTQ structure

    Total line count must be divisible by 4.

    Each FASTQ record must follow the standard 4-line format.

4. Paired-end consistency

    R1 and R2 must contain the same number of records.

If record counts differ, the pipeline raises:
```
ValueError: R1/R2 record counts differ
```

## Multiplex CSV Requirements

The CSV file must contain the following columns:
```
Population
R1_index
R2_index
```

Example:
| Population | R1_index  | R2_index  |
|------------|-----------|-----------|
| P1         | ACGTACGT  | TTGCAACC  |
| P2         | GTCAGTCA  | AACCTTGG  |

## Multiplex Validation Rules

1. At least one data row must exist.

2. Required columns must be present.

Population must match:
```
^[A-Za-z0-9._-]+$
```
Each (R1_index, R2_index) pair must be unique.

Duplicate index pairs will raise an error and stop execution.

## Terminal Output

If validation succeeds, Stage 0 prints a structured summary:
```
========== Stage 0 Summary ==========
Time         : 2026-03-03 11:22:01
Status       : OK
-------------------------------------
FASTQ Pair
  R1          : sample_R1.fastq.gz
  R2          : sample_R2.fastq.gz
  R1 records  : 1250000
  R2 records  : 1250000
-------------------------------------
Multiplex CSV
  Path        : multiplex.csv
  Rows        : 8
  Populations : 2
  Names       : P1, P2
  Index mode  : R1+R2
=====================================
```

If validation fails, the pipeline prints:
```
Stage 0 FAILED.
Error: <error message> (e.g. missing file, mismatched record counts, duplicate index pairs)
```
## Generated Files

If Stage 0 succeeds:

  - results/<exp>/metrics/metrics_stage0_validation.json

    Structured validation result (full JSON output)

  - results/<exp>/run_manifest.json
  
    Updated with Stage 0 status and result

No additional data files are created during Stage 0.