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

This is the recommended approach as it installs all dependencies including bowtie2 (required for Step 6).

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

**Note:** You must also have bowtie2 installed on your system for Step 6:
- macOS: `brew install bowtie2`
- Linux: `sudo apt-get install bowtie2`

## Step 0: step0_file_confirmation.py

Confirms that uploaded experiment files are present, correctly named, and sufficient for downstream pipeline stages.

**Usage:**
```bash
python -m biongs.step0_file_confirmation <experiment_name>
```

**Example:**
```bash
python -m biongs.step0_file_confirmation example
```

**Input:**
- Experiment name
- Expected folder: `input_data/{experiment_name}/{experiment_name}_fastq/`
- Expected CSV files in `input_data/{experiment_name}/`:
  - `{experiment_name}_multiplexing_info.csv`
  - `{experiment_name}_UMI_primers.csv`

**Output:**
- Terminal report showing:
  - Which required files/folders are present or missing
  - Detected R1/R2 fastq pairs
  - Furthest runnable pipeline stage based on uploaded inputs
  - Possible filename typos for near-matching files
- Saved report file:
  - `results/{experiment_name}/logs/file_confirmation_report.txt`

---

## Step 0a: step0a_fast_qc.py

Runs FastQC on all raw FASTQ inputs and keeps only the HTML quality reports.

**Usage:**
```bash
python -m biongs.step0a_fast_qc <experiment_name>
```

**Example:**
```bash
python -m biongs.step0a_fast_qc example
```

**Input:**
- Experiment name
- FASTQ files in `input_data/{experiment_name}/{experiment_name}_fastq/`

**Output:**
- One HTML report per FASTQ file in `results/{experiment_name}/qc/`
- Named `{fastq_name}_fastqc.html`

**Dependencies:**
- FastQC must be installed. See https://github.com/s-andrews/FastQC/blob/master/INSTALL.md
- Step 0 is recommended to confirm file locations before running

---

## Step 1: step1_demultiplex_folders.py

Creates the folder structure for demultiplexed populations.

**Usage:**
```bash
python step1_demultiplex_folders.py <experiment_name>
```

**Example:**
```bash
python step1_demultiplex_folders.py example
```

**Input:**
- Experiment name (script will look for a multiplexing CSV in `input_data/{experiment_name}/`)
- CSV file should contain columns: Time, Population, GW_name, R1_index, R2_index

**Output:**
- Creates output directory at `results/{experiment_name}/demultiplexing/`
- Creates a folder for each population (named P{Population}, e.g., "P1" for Population 1)
- Creates empty R1 and R2 fastq files in each folder
- Saves terminal output log to `results/{experiment_name}/logs/demultiplex_folders_terminal_output.txt`

---

## Step 2: step2_demultiplex_index.py

Sorts NGS reads by DNA index into population-specific files. For each population in the CSV, finds fastq files matching that population's GW_name and extracts reads with matching indexes.

**Usage:**
```bash
python step2_demultiplex_index.py <experiment_name>
```

**Example:**
```bash
python step2_demultiplex_index.py example
```

**Input:**
- Experiment name
- Script reads the multiplexing CSV to find all populations and their GW_names
- For each population, finds fastq files in `input_data/{experiment_name}/` matching the GW_name
- Processes each population's fastq files independently
- Supports both regular (`*.fastq`, `*.fq`) and gzipped (`*.fastq.gz`, `*.fq.gz`) input files
- Automatically decompresses gzipped files before processing

**Output:**
- Populates population-specific fastq files with matched reads in `results/{experiment_name}/demultiplexing/P{Population}/`
- Creates files for short reads (< 150 bp): `{GW_name}_short_reads_R1.fastq` and `{GW_name}_short_reads_R2.fastq`
- Creates files for unmatched reads: `{GW_name}_unmatched_reads_R1.fastq` and `{GW_name}_unmatched_reads_R2.fastq`
- Prints summary statistics for each population showing read distribution
- Saves terminal output log to `results/{experiment_name}/logs/demultiplex_index_terminal_output.txt`

**Quality Checks:**
- Verifies R1 and R2 read pairs match by checking headers
- Filters reads shorter than 150 base pairs
- Separates reads with indexes that don't match any population

---

## Step 2a: step2a_check_index_quality.py

Analyzes the quality of demultiplexing results from Step 2 and generates visualizations of read distributions and quality metrics.

**Usage:**
```bash
python step2a_check_index_quality.py <experiment_name>
```

**Example:**
```bash
python step2a_check_index_quality.py example
```

**Input:**
- Experiment name
- Reads from the demultiplexing results directory: `results/{experiment_name}/demultiplexing/`
- Analyzes:
  - Input fastq files from `input_data/{experiment_name}/`
  - Population fastq files in each `P{Population}/` folder
  - Short reads files (`{GW_name}_short_reads_R1.fastq` and `_R2.fastq`)
  - Unmatched reads files (`{GW_name}_unmatched_reads_R1.fastq` and `_R2.fastq`)

**Output:**
- Creates output directory at `results/{experiment_name}/index_quality/`
- **Read length histogram**: Single PNG plot (`{experiment_name}_read_length_histogram.png`) with multi-panel layout
  - Separate panels for each input fastq file (columns)
  - Separate rows for R1 and R2 reads
  - Samples up to 10,000 reads from each input file
- **Read distribution bar plot**: Single PNG plot (`{experiment_name}_read_distribution.png`) with multi-panel layout
  - Separate panels for each fastq file showing read distribution across:
    - Each population
    - Short reads category
    - Unmatched reads category
  - Displays read counts and percentages on bars
- **Analysis CSV file**: Tidy R-compatible CSV (`index_quality_results.csv`) with columns:
  - `experiment`, `sample`, `category_type`, `category`, `read_count`, `input_total`, `percentage_of_input`
  - Includes rows for populations, QC categories, and summary totals
  - Easily importable into R with `read.csv()`
- **Summary report**: Text file (`index_quality_summary.txt`) with detailed statistics including:
  - Total input reads
  - Reads per population with percentages
  - Short reads count and percentage
  - Unmatched reads count and percentage
  - Recovery rate (percentage of reads successfully categorized)

**Visualizations:**
- Uses Altair library for high-quality static visualizations
- All plots saved as PNG files with multiple panels for easy comparison (with HTML fallback if PNG export fails)

**Dependencies:**
- Requires Step 2 (step2_demultiplex_index.py) to be run first
- Python packages: altair, pandas, biopython, vl-convert-python (for PNG export)

---

## Step 3: step3_demultiplex_UMI.py

Detects and extracts UMI (Unique Molecular Identifier) sequences from demultiplexed reads, creating libraries organized by UMI pairs for each population.

**Usage:**
```bash
python step3_demultiplex_UMI.py <experiment_name>
```

**Example:**
```bash
python step3_demultiplex_UMI.py example
```

**Input:**
- Experiment name (script will auto-detect multiplexing and UMI primers CSVs from `input_data/{experiment_name}/`)
- Multiplexing CSV: File with 'multiplexing' in the name containing population information
- UMI Primers CSV: File with 'UMI' or 'primer' in the name containing primer sequences
  - CSV should have columns 'f' (forward primer, matches R1) and 'r' (reverse primer, matches R2)
  - Each primer should contain exactly 10 consecutive N's indicating the UMI position

**Output:**
- Creates a UMI library pickle file for each population: `P{Population}_UMI_dict.pkl`
- Library stored in `results/{experiment_name}/demultiplexing/P{Population}/`
- Library structure: Dictionary with (forward_UMI, reverse_UMI) tuples as keys
  - Each value contains:
    - `'R1'`: List of BioPython SeqRecord objects (full untrimmed R1 reads)
    - `'R2'`: List of BioPython SeqRecord objects (full untrimmed R2 reads)
    - `'R1_trim_pos'`: Position where R1 should be trimmed (after primer + UMI)
    - `'R2_trim_pos'`: Position where R2 should be trimmed (after primer + UMI)
- SeqRecord objects preserve original read IDs, sequences, and per-base quality scores
- Trimming is deferred to Step 4 (step4_alignment_prep.py) for flexibility
- R1 and R2 sequences are in matched order by read pair
- Also creates files for unmatched reads: `P{Population}_unmatched_UMI_R1.fastq` and `_R2.fastq`
- Saves terminal output log to `results/{experiment_name}/logs/demultiplex_UMI_terminal_output.txt`

**Processing Details:**
- UMI extraction: Matches primer sequences before and after the 10 N's to find UMI location
- Stores full untrimmed SeqRecord objects along with calculated trim positions
- Trimming is performed downstream in Step 6 to preserve flexibility
- Preserves all read metadata including IDs and quality scores
- Verifies R1 and R2 read pairs match by checking headers

---

## Step 3a: step3a_check_UMI_quality.py

Analyzes the quality of UMI dictionaries created by Step 4 and generates visualizations comparing UMI distributions across populations.

**Usage:**
```bash
python step3a_check_UMI_quality.py <experiment_name>
```

**Example:**
```bash
python step3a_check_UMI_quality.py example
```

**Input:**
- Experiment name
- UMI dictionary pickle files from `results/{experiment_name}/demultiplexing/P{Population}/`
- Reads all files matching pattern `*_UMI_dict.pkl`

**Output:**
- Creates output directory at `results/{experiment_name}/UMI_quality/`
- **UMI count bar plot**: PNG plot comparing the number of unique UMI pairs across populations
- **Reads per UMI violin plot**: PNG plot showing distribution of read counts per UMI for each population
- **Summary statistics**: Displays mean and median reads per UMI for each population
- **Text summary**: `UMI_quality_summary.txt` with detailed statistics

**Quality Metrics:**
- Number of unique UMI pairs per population
- Distribution of reads per UMI (to assess sequencing depth and coverage)
- Populations with very few UMIs or highly uneven read distribution may indicate quality issues

**Visualizations:**
- Uses Altair library for high-quality static visualizations
- All plots saved as PNG files

**Dependencies:**
- Requires Step 3 (step3_demultiplex_UMI.py) to be run first
- Python packages: altair, pandas, pickle, vl-convert-python (for PNG export)

---

## Step 4: step4_alignment_prep.py

Prepares UMI-grouped sequences for downstream alignment by trimming stored SeqRecord objects and organizing them into per-UMI FASTQ files.

**Usage:**
```bash
python step4_alignment_prep.py <experiment_name> [--min-reads-per-umi N]
```

**Example:**
```bash
# Include all UMI pairs with at least 1 read (default)
python step4_alignment_prep.py example

# Only process UMI pairs with at least 5 reads
python step4_alignment_prep.py example --min-reads-per-umi 5

# Only process UMI pairs with at least 10 reads
python step4_alignment_prep.py example --min-reads-per-umi 10
```

**Input:**
- Experiment name
- UMI dictionary pickle files from `results/{experiment_name}/demultiplexing/P{Population}/`
- Reads all files matching pattern `*_UMI_dict.pkl`
- Optional: `--min-reads-per-umi` threshold (default: 2)

**Output:**
- Creates output directory at `results/{experiment_name}/alignment/fastq/`
- For each population, creates a subfolder: `P{Population}/`
- For each UMI pair meeting the minimum read threshold, creates two FASTQ files:
  - `{forward_UMI}_{reverse_UMI}_R1.fastq`: Trimmed R1 reads for that UMI pair
  - `{forward_UMI}_{reverse_UMI}_R2.fastq`: Trimmed R2 reads for that UMI pair
- Sequences are trimmed using trim positions stored in Step 3
- Preserves original read IDs and per-base quality scores from input data
- UMI pair information added to sequence description line
- Saves terminal output log to `results/{experiment_name}/logs/alignment_prep_terminal_output.txt`

**Filtering:**
- Only processes UMI pairs with at least `--min-reads-per-umi` sequences (default: 2)
- Skips UMI pairs below threshold to reduce file count and focus on well-supported variants
- Reports statistics on UMI pairs processed vs. skipped

**Processing Details:**
- Loads SeqRecord objects from UMI dictionaries
- Trims sequences using stored `R1_trim_pos` and `R2_trim_pos` values
- Slicing SeqRecord objects preserves per-base quality annotations
- Backward compatible with older string-based dictionaries (if present)

**Quality Preservation:**
- Original FASTQ quality scores maintained through trimming
- Read IDs preserved from source files
- UMI pair metadata added to description field

**Dependencies:**
- Requires Step 3 (step3_demultiplex_UMI.py) to be run first
- Python packages: biopython, pickle

---

## Complete Pipeline Example

```bash
# Navigate to src directory
cd src

# Step 0: Confirm files are present and correctly named
python -m biongs.step0_file_confirmation example

# Step 0a: Run FastQC quality checks (requires FastQC installed)
python -m biongs.step0a_fast_qc example

# Step 1: Create folder structure
python step1_demultiplex_folders.py example

# Step 2: Demultiplex reads by index for all populations
python step2_demultiplex_index.py example

# Step 2a: Check quality of demultiplexing results
python step2a_check_index_quality.py example

# Step 3: Create UMI libraries for each population
python step3_demultiplex_UMI.py example

# Step 3a: Check quality of UMI data
python step3a_check_UMI_quality.py example

# Step 4: Prepare sequences for alignment (trim and organize by UMI)
python step4_alignment_prep.py example --min-reads-per-umi 5
```
