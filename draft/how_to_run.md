# Overview of fastqc_unpacker.py and entry_qc.py

---

# How to run?
run this 3 rows sequentially
1. git clone
2. cd /home/uniscan/BioNGS
3. python draft/entry_qc.py

- you will get a new folder "raw_qc_example" in /home/uniscan/BioNGS/draft
- There are 4 new items in folder "raw_qc_example": R1  R2  fastqc_out  index.html 
- R1 & R2 are QCdashboard reporting quality of R1 and R2 example fastq data
- fastqc_out is the original data from the software fastqc
- index.html is the shortcut to R1 and R2 QCdashboard, with superlink in it just like portal.



# fastqc_unpacker.py — QC Report Generator

## Purpose

fastqc_unpacker.py is responsible for parsing FastQC output files and generating a structured HTML quality control (QC) report.

---

## Functionality

Reads FastQC output files, including:

- `fastqc_data.txt`  
- `Images/*.png`  

Parses QC statistics such as:

- Total reads  
- Read length distribution  
- Sequence quality scores  
- Module pass/warn/fail status  

Computes additional derived metrics:

- Q30 proxy and other summary metrics  

Embeds FastQC images into the report.

Generates a formatted HTML dashboard:

- `raw_qc_dashboard.html`

---

## Input
```
<sample>_fastqc/
├── fastqc_data.txt
└── Images/
```
---

## Output
raw_qc_dashboard.html

---

## Role in the pipeline

fastqc_unpacker.py acts as a **report generation module**.

It does NOT:

- Run FastQC  
- Locate FASTQ files  
- Perform sequencing analysis  

It only processes existing FastQC output and converts it into a structured, human-readable HTML report.

---

# entry_qc.py — QC Pipeline Entry Point

## Purpose

entry_qc.py is the main entry point of the raw QC pipeline. It orchestrates the workflow from FASTQ input to final HTML report generation.

---

## Functionality

- Locates input FASTQ files (R1 and optionally R2)  
- Runs FastQC on each FASTQ file  
- Extracts FastQC output archives (`*_fastqc.zip`)  
- Identifies extracted FastQC result directories  
- Calls fastqc_unpacker.py to generate QC HTML reports  
- Organizes outputs into a structured results directory  
- Optionally generates an index page linking multiple reports  

---

## Input
```
FASTQ files:
sample_R1.fastq
sample_R2.fastq
```

---

## Output
```
results/
├── fastqc_out/
├── R1/raw_qc_dashboard.html
├── R2/raw_qc_dashboard.html
└── index.html
```

---

## Role in the pipeline

entry_qc.py acts as the **pipeline controller (entry point)**.

It is responsible for:

- Managing workflow execution  
- Running analysis tools (FastQC)  
- Connecting analysis outputs to report generation  

---

# Relationship Between the Two Components

## Workflow Diagram

```
FASTQ files
↓
entry_qc.py (pipeline controller)
↓
FastQC (quality analysis tool)
↓
fastqc_data.txt + Images/
↓
fastqc_unpacker.py (report generator)
↓
HTML QC report
```

---

## Software Architecture Roles

| Component | Role |
|----------|------|
| entry_qc.py | Pipeline Controller / Entry Point |
| FastQC | Analysis Engine |
| fastqc_unpacker.py | Report Generator / Rendering Module |


