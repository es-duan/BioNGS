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

## Script 1: demultiplex_folders.py

Creates the folder structure for demultiplexed populations.

**Usage:**
```bash
python demultiplex_folders.py <experiment_name>
```

**Example:**
```bash
python demultiplex_folders.py example
```

**Input:**
- Experiment name (script will look for a multiplexing CSV in `input_data/{experiment_name}/`)
- CSV file should contain columns: Time, Population, GW_name, R1_index, R2_index

**Output:**
- Creates output directory at `Outputs/{experiment_name}/demultiplexing/`
- Creates a folder for each population (named P{Population}, e.g., "P1" for Population 1)
- Creates empty R1 and R2 fastq files in each folder

---

## Script 2: demultiplex_index.py

Sorts NGS reads by DNA index into population-specific files. For each population in the CSV, finds fastq files matching that population's GW_name and extracts reads with matching indexes.

**Usage:**
```bash
python demultiplex_index.py <experiment_name>
```

**Example:**
```bash
python demultiplex_index.py example
```

**Input:**
- Experiment name
- Script reads the multiplexing CSV to find all populations and their GW_names
- For each population, finds fastq files in `input_data/{experiment_name}/` matching the GW_name
- Processes each population's fastq files independently

**Output:**
- Populates population-specific fastq files with matched reads in `Outputs/{experiment_name}/demultiplexing/P{Population}/`
- Creates files for short reads (< 150 bp): `{GW_name}_short_reads_R1.fastq` and `{GW_name}_short_reads_R2.fastq`
- Creates files for unmatched reads: `{GW_name}_unmatched_reads_R1.fastq` and `{GW_name}_unmatched_reads_R2.fastq`
- Prints summary statistics for each population showing read distribution

**Quality Checks:**
- Verifies R1 and R2 read pairs match by checking headers
- Filters reads shorter than 150 base pairs
- Separates reads with indexes that don't match any population

---

## Script 4: demultiplex_UMI.py

Detects and extracts UMI (Unique Molecular Identifier) sequences from demultiplexed reads, creating libraries organized by UMI pairs for each population.

**Usage:**
```bash
python demultiplex_UMI.py <experiment_name>
```

**Example:**
```bash
python demultiplex_UMI.py example
```

**Input:**
- Experiment name (script will auto-detect multiplexing and UMI primers CSVs from `input_data/{experiment_name}/`)
- Multiplexing CSV: File with 'multiplexing' in the name containing population information
- UMI Primers CSV: File with 'UMI' or 'primer' in the name containing primer sequences
  - CSV should have columns 'f' (forward primer, matches R1) and 'r' (reverse primer, matches R2)
  - Each primer should contain exactly 10 consecutive N's indicating the UMI position

**Output:**
- Creates a UMI library pickle file for each population: `P{Population}_UMI_library.pkl`
- Library stored in `Outputs/{experiment_name}/demultiplexing/P{Population}/`
- Library structure: Dictionary with (forward_UMI, reverse_UMI) tuples as keys
  - Each value contains:
    - `'R1'`: List of trimmed R1 sequences (primers removed)
    - `'R2'`: List of trimmed R2 sequences (primers removed)
- R1 and R2 sequences are in matched order by read pair

**Processing Details:**
- UMI extraction: Matches primer sequences before and after the 10 N's to find UMI location
- Sequences trimmed to remove primer and UMI regions
- Verifies R1 and R2 read pairs match by checking headers

---

## Complete Pipeline Example

```bash
# Navigate to src directory
cd src

# Step 1: Create folder structure
python demultiplex_folders.py example

# Step 2: Demultiplex reads by index for all populations
python demultiplex_index.py example

# Step 4: Create UMI libraries for each population
python demultiplex_UMI.py example
```
