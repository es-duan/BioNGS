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

See [src/README.md](src/README.md) for detailed usage instructions for each script.
