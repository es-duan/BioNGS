# tentitive software structrue proposal

## Overview
BioNGS Pipeline
├── Input
│   ├── FASTQ files
│   └── Reference sequence
├── Processing
│   ├── Raw data Quality Control (RQC)
│   ├── Split fastq by index（Demultiplexing）
│   ├── UMI extraction and storing
│   ├── Trimming
│   │   ├──UMI and Index
│   │   └──low quality and extreme short reads
│   ├── UMI collapsing / grouping
│   ├── Post-trim QC
│   ├── Alignment
│   └── Variant Detection
├── Inference
│   └── Mutation Rate Calculation
└── Report
    ├── QC Summary
    └── Mutation Results
----
## Function zoom in at different layers
### 1. Input
----
### 2. Processing
#### 2.1. Raw data Quality Control (RQC)

#### 2.2. Split fastq by index（Demultiplexing）

#### 2.3. UMI extraction and storing

#### 2.4. Trimming

#### 2.5. UMI collapsing / grouping

#### 2.6. Post-trim QC

#### 2.7. Alignment

#### 2.8. Variant Detection
----

### 3. Inference
#### 3.1. Mutation Rate Calculation
----

### 4. Report