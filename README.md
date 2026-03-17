# BioNGS
Repo for CHEM E 546 BioNGS Project. 

Pipeline for processing and analyzing Next Generation Sequencing data for bacterial mutation rate analysis.

## Setup

### Install the package and dependencies using Conda (Recommended)

This project requires bowtie2 and other bioinformatics tools. The easiest way to set up the entire environment is using conda:

```bash
# Create the conda environment from the environment.yml file
conda env create -f environment.yml

# Activate the environment
conda activate biongs
```

This will install all required dependencies including:
- Python 3.11
- BioPython
- Bowtie2
- Matplotlib
- NumPy
- Pandas

### Alternative: Manual Setup with pip and system tools

If you prefer pip and already have bowtie2 installed on your system:

```bash
# Create a virtual environment
python -m venv .venv
source .venv/bin/activate  # On macOS/Linux
# .venv\Scripts\activate   # On Windows

# Install Python dependencies
pip install -e .

# Note: You must also install bowtie2 separately:
# macOS: brew install bowtie2
# Linux: sudo apt-get install bowtie2
# Windows: See bowtie2 installation guide
```

## Usage

### Install from pyproject.toml

From the repository root:

```bash
# (optional) activate your environment first
source .venv/bin/activate

# install package in editable mode
pip install -e .
```

This installs the CLI commands defined in `pyproject.toml`.

### Run the pipeline steps (terminal commands)

Replace `example` with your experiment folder name under `input_data/`.

```bash
# Step 0: confirm required files and detect runnable stages
file-confirmation example

# Step 0a: run FastQC on raw reads
fast-qc example

# Step 1: create demultiplexing folder structure
demultiplex-folders example

# Step 2: demultiplex by index (optional threshold flag)
demultiplex-index example --short-read-length 150

# Step 2a: evaluate index/demultiplexing quality
check-index-quality example

# Step 3: build UMI dictionaries
demultiplex-umi example

# Step 3a: evaluate UMI quality
check-umi-quality example

# Step 4: prepare trimmed UMI fastqs for alignment
alignment-prep example --min-reads-per-umi 2
```

### Step 5 alignment command

Step 5 is a shell script and is run directly:

```bash
bash src/biongs/step5_align_reads.sh \
	example \
	input_data/reference_docs/rpoB_index/rpoB_index \
	input_data/reference_docs/rpoBsequence.fasta
```

### Notes

- If a command is not found after install, re-activate your environment and run `pip install -e .` again.
- You can still run modules directly (for example: `python -m biongs.step0_file_confirmation example`).
- See [src/README.md](src/README.md) for detailed per-step behavior and outputs.
