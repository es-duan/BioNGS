# tentitive software structrue proposal

## Overview
BioNGS Pipeline
├── Input
│   ├── FASTQ files
│   └── Reference sequence
├── Processing
│   ├── Raw data Quality Control (Raw QC)
│   └── Split fastq by index（Demultiplexing）
│       └── unmatch rate report
│   ├── Per-population Quality Control (QC for each populations) 
│   ├── UMI extraction and storing
│   ├── Trimming
│   │   ├──UMI and Index
│   │   └──low quality and extreme short reads
│   ├── UMI collapsing / grouping
│   ├── Post-trim QC (final check)
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
#### 2.1. Raw data Quality Control (Raw QC)
- The first report should be generated using the raw data, before any trimming or demultiplexing is performed. Its purpose is to provide an overall assessment of whether the sequencing dataset is fundamentally compromised. 
  - For example, it should help identify issues such as globally low read quality, severe adapter contamination, abnormal read length distributions, or an excessive proportion of short reads.

- I do not think this stage is appropriate for determining trimming thresholds. Because the samples are mixed, reads originating from different populations may have distinct characteristics. Averaging across mixed populations could mask population-specific problems and lead to misleading trimming decisions.

##### which will be report?
  - Basic Statistics(table) 
    - Total Sequences(reads)
    - Total Bases
    - Sequence length(a range) + median or IQR
    - Q30 ratio
    - Sequences flagged as poor quality

  - pass/warning/fail dignosis from each section of FastQC(table):
    - Per base sequence quality
    - Per tile sequence quality
    - Per sequence quality scores
    - Per base N content
    - Sequence Length Distribution
    - Adapter Content
    - Sequence Duplication Levels
    - Overrepresented sequences

  - plots of the festure above


#### 2.2. Split fastq by index（Demultiplexing）

#### 2.3. per-population Quality Control (QC for each populations) 
- Which tube is wrong.
  - For example: population 12 is poor in average base quality, population 7 has low reads quantity 
- **evidence for trimming** : if there is no diff between each population(stratification)
  - no diff can be: per-base Q same and the read length distribution same

##### which will be report? 
1. whole fastq file after identifing the unmatched reads
  - total read number
  - matched reads vs unmatched reads vs shortreads(reads shorter than median of the raw data which will be set as default) (bar chart)

2. each population fastq file(only for the matched reads)
- pass/warning/fail dignosis from each section of FastQC(big table with all population for comparison)
- read quantity for every population and short reads **(threshold is reads shorter than median of the matched data in each group which will be set as default)**
  - plot all plot below and compair them:
    - Per base sequence quality
    - Per sequence quality scores
    - Per base N content
    - Sequence Length Distribution
    - Sequence Duplication Levels
    - Overrepresented sequences
#### WINDOW: should automatic trimming be stopped to set populaton-specific parameters?

#### 2.4. UMI extraction and storing

#### 2.5.1. automatic Trimming by total
**base on the report above, manual and automatic options can be chosen**
- **question: do you want to drop the low average Q score reads (mean Q score < threshold = 30 as company sat) OR cut bases from 3' to 5' gradually by agrithems**
- cut UMI and Index

#### 2.5.2. manual Trimming by population


#### 2.5.1. automatic filter by total
- drop the shorter reads < threshold 
  - calculate total read length distributions' median
  - discard the read
- drop the reads with N > threshold

#### 2.5.2. manual filter by population
- drop the shorter reads < threshold 
  - for each population
    - calculate read length distributions' median
    - min_len = max( floor(median - 1.5 × IQR) or floor(median), hard_min_len )
      - if population reads < gate, go back to global median reads and report
- drop the reads with N > threshold

#### 2.6. UMI collapsing / grouping

#### 2.7. Post-trim QC (final check)
- the key is : prove the trimming decision was effective and reasonable, which feeds a clean data to the alignment script.
  - the reads quantity of each population  +  2.7 read count/2.3 read count
  - plot all plot below and compair them:
    - Per base sequence quality
    - Per sequence quality scores
    - Per base N content
    - Sequence Length Distribution
    - Sequence Duplication Levels
    - Overrepresented sequences

???how can check whether UMI and INDEX been deleted completely?
#### 2.8. Alignment

#### 2.9. Variant Detection
----

### 3. Inference
#### 3.1. Mutation Rate Calculation
----

### 4. Report