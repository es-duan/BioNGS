# Metrics Explanation

This document explains the key metrics, summary fields, warning thresholds, and report outputs used in the BioNGS pipeline.

The goal of these metrics is to make each stage interpretable, reproducible, and easy to troubleshoot.

## 1. General conventions

### 1.1 Read-level vs pair-level metrics

BioNGS uses both read-level and pair-level metrics.

A read-level metric treats each read independently, for example:

- R1 short rate

- R2 short rate

- base-level Q30 rate

A pair-level metric treats one R1/R2 pair as a single unit, for example:

- total read pairs

- unmatched pairs

- usable pairs

- short pairs in Stage 5 preview

This distinction is important because a metric such as short_rate may be calculated differently in different stages.

### 1.2 Base-level Q30 rate

Unless otherwise stated, Q30 rate means:
```
fraction of bases with Phred quality score greater than or equal to 30.
```
In the code, this is computed by counting all bases whose ASCII-decoded quality score is at least 30, and dividing by the total number of bases.

Formula:
```
Q30 rate = q30_bases / total_bases
```
Interpretation:

- higher values indicate better sequencing quality

- this is a base-level metric, not a read-level pass/fail percentage

### 1.3 Length distribution

Several stages record a length_counter, which is a mapping:
```
read length (bp) → number of reads
```
This is used to generate:

- histogram-style read length plots in raw / intermediate QC

- step-line read length frequency plots in post-trim QC

## 2. Stage 0 — Input validation metrics

Stage 0 checks whether the input data are structurally valid before the pipeline starts. It validates:

- existence of R1 and R2 files

- gzip readability for .gz files

- FASTQ record counts

- paired-end consistency

- multiplex CSV structure

- uniqueness of dual-index pairs (R1_index, R2_index)

## Key validation outputs
`r1_records, r2_records`

Number of FASTQ records in R1 and R2.

Interpretation:

- these should be equal for paired-end input

- if they differ, the pipeline stops before downstream processing 

`n_rows`

Number of rows in the multiplex CSV after parsing.

`unique_populations`

List of distinct population names found in the multiplex CSV.

`n_populations`

Number of distinct populations defined in the multiplex CSV.

`index_key_mode`

Indicates how populations are keyed in the multiplex mapping. In the current implementation this is:
```
R1+R2
```
meaning demultiplexing is based on a dual-index pair rather than a single inline index.

## 3. Stage 1 — Raw QC metrics

Stage 1 generates two QC layers:

- QC details from FastQC

- QC overview computed directly from FASTQ

### Per-read-file summary fields

These are produced separately for R1 and R2.

`total_reads`

Number of reads in that FASTQ file. Computed by iterating through records in the file. 

`q30_rate`

Base-level Q30 rate.

Formula:
```
q30_rate = q30_bases / total_bases
```

`short_rate`

Fraction of reads whose read length is shorter than the short_len threshold.
In Stage 1, `short_len` is fixed at 150 bp.

Formula:
```
short_rate = number of reads with length < 150 / total_reads
```

`length_counter`

Mapping from read length to read count.

This field is used to draw `length_distribution.png`.

## Combined Stage 1 decision metrics

Stage 1 also defines a conservative combined decision rule across R1 and R2.

`short_rate_used`

The worse short-rate value between R1 and R2.

Formula:
```
short_rate_used = max(R1 short_rate, R2 short_rate)
```
`q30_rate_used`

The worse Q30 rate between R1 and R2.

Formula:
```
q30_rate_used = min(R1 q30_rate, R2 q30_rate)
```
This means the pipeline judges raw input quality based on the weaker mate rather than the average.

### Stage 1 warning thresholds

Stage 1 uses the following logic:
```
short_rate > 30% → warning
Q30 rate < 70% → warning
short_rate > 30% AND Q30 rate < 70% → abnormal
```
Interpretation:

- a high short-rate suggests many truncated reads

- a low Q30 rate suggests generally poor sequencing quality

- if both occur together, the pipeline treats the data as abnormally poor and prompts the user to stop

## 4. Stage 2 — Demultiplexing metrics

Stage 2 assigns each read pair to a population using a dual-index mapping:
```
(R1_index, R2_index) → population
```
A pair is assigned only if the R1 read starts with the expected `R1_index` and the R2 read starts with the expected `R2_index`. Unmatched pairs are written to dedicated unmatched FASTQ files. By default, matched reads are not stripped of index sequence in this core demultiplex function because `strip_matched=False`.

## Main Stage 2 metrics

`total_pairs`

Total number of paired reads processed in demultiplexing.

`unmatched_pairs`

Number of pairs that could not be assigned to any population based on the dual-index mapping. These pairs are written unchanged to unmatched output files.

`unmatched_rate`

Fraction of read pairs that could not be assigned to a population.

Formula:
```
unmatched_rate = unmatched_pairs / total_pairs
```
Interpretation:

- high values may indicate wrong index definitions

- high values may also indicate orientation problems or library structure mismatch

- unmatched reads are preserved as evidence for troubleshooting rather than discarded silently

`per_population_pairs`

Dictionary mapping population name to the number of matched read pairs assigned to that population.

`index_lengths`

Summary of expected index lengths in the multiplex mapping. Current fields include:

- r1_unique_lengths

- r2_unique_lengths

- r1_max_len

- r2_max_len 

These values are used internally for prefix matching and are also helpful for validating whether the multiplex CSV is internally consistent.

## Stage 2 warning thresholds
```
unmatched_rate > 20% → warning
unmatched_rate > 35% → strong warning
unmatched_rate > 50% → abnormal, user may stop
```

Interpretation:

- low unmatched rates are expected in a well-configured demultiplex run

- high unmatched rates are usually a configuration or library-structure problem, not just noise

## 5. Stage 2.5 — Post-demultiplex QC metrics

Stage 2.5 evaluates the FASTQ files generated by demultiplexing, including per-population outputs and optionally the unmatched bin. It produces both:

- per-population raw-style QC summaries

- an overall demultiplex distribution summary 

### Per-population QC overview metrics

For each population and mate, Stage 2.5 reuses compute_raw_overview, so the per-population QC summaries include:

- total_reads

- short_rate

- q30_rate

- length_counter

These metrics should be interpreted exactly as in Stage 1, except that now they apply to each demultiplexed population separately.

### Demultiplex summary metrics

`pop_pair_counts`

Population-level counts derived from the number of reads in the demultiplexed R1 files. Since each R1 read corresponds to one paired record, this serves as a pair-count proxy. 

`total_pairs`

Sum of all per-population counts including unmatched.

`unmatched_pairs`

Count of pairs in the unmatched output.

`unmatched_rate`

Fraction of pairs that remained unmatched after demultiplexing.

Formula:
```
unmatched_rate = unmatched_pairs / total_pairs
```
This value should be broadly consistent with Stage 2.

## 6. Stage 3 — UMI extraction metrics

Stage 3 parses UMI-containing primer structures from demultiplexed reads. The current implementation expects the primer CSV to define one forward and one reverse primer pattern of the form:
```
prefix + 10N UMI + suffix
```
The code locates the prefix, extracts a fixed 10-base UMI, verifies the suffix at the expected position, and trims off the structural sequence so that only the insert remains. 

### Core extraction outcomes

For each population, Stage 3 tracks at least the following high-level counts.

`total_pairs`

Number of paired reads examined for UMI extraction in that population.

`extracted_pairs`

Number of pairs for which UMI extraction succeeded.

Interpretation:

- both mates satisfied the required parsing logic

- trimmed output FASTQ files are written for these pairs


`unmatched_pairs`

Number of pairs for which UMI extraction failed.

These failed pairs are written to unmatched UMI FASTQ outputs. 

`unmatched_rate`

Fraction of pairs that failed UMI extraction.

Formula:
```
unmatched_rate = unmatched_pairs / total_pairs
```
This is the broadest failure summary for Stage 3.

### Why Stage 3 separates short pairs from usable pairs

A read may fail UMI extraction not because the pattern is wrong, but because the read is too short to contain the full expected structure. To distinguish these situations, Stage 3 uses a reporting-only threshold:
```
SHORT_MINLEN = 150
```
Pairs are split into two categories:
```
short_pairs_lt150

usable_pairs_ge150 
```

`short_pairs_lt150`

Number of pairs where at least one mate is shorter than 150 bp.

Interpretation:

- these pairs may be structurally incomplete

- failure on these pairs is less informative about true parsing performance

`usable_pairs_ge150`

Number of pairs where both mates are at least 150 bp.

Interpretation:

- these pairs are long enough to be considered structurally interpretable under the current reporting rule

`usable_unmatched_pairs`

Number of failed UMI-extraction pairs among usable_pairs_ge150.

`usable_unmatched_rate`

Failure rate among usable pairs only.

Formula:
```
usable_unmatched_rate = usable_unmatched_pairs / usable_pairs_ge150
```
Interpretation:

- this metric is usually more informative than overall unmatched_rate

- it helps separate genuine parsing failure from short-read failure caused by insufficient read length

This is one of the most useful diagnostic metrics in Stage 3.

### Per-mate and pair-level failure reasons

Stage 3 writes a `umi_extraction_summary.txt` for each population. The summary includes:

- pair-level fail reasons

- R1 fail reasons

- R2 fail reasons

`fail_reason_pair`

Counter of pair-level first failure reasons.

Interpretation:

- summarizes the dominant reason a pair failed extraction

- useful for identifying the main bottleneck

`fail_reason_r1`

Counter of failure reasons observed on R1.

`fail_reason_r2`

Counter of failure reasons observed on R2.

These mate-specific counters help distinguish asymmetric problems, for example when one mate fails mostly at anchor matching while the other mate succeeds more often.

### Common Stage 3 failure reason labels

The exact labels may depend on which extractor path is used, but the uploaded code clearly includes the following possible reason names:

`qual_len_mismatch`

Sequence and quality string lengths do not match.

This indicates malformed FASTQ content and is not a biological failure. 

`too_short`

The read does not extend far enough to cover the expected prefix + UMI + suffix region. 

`missing_prefix_left / missing_prefix_anywhere`

The expected primer prefix could not be found in the required location or anywhere in the read, depending on the extraction mode. 

`missing_anchor_right`

The expected suffix anchor after the UMI was not found at the required position. This often indicates one of the following:

- structural mismatch

- shifted primer placement

- read too short to fully span the region

- wrong orientation or wrong primer definition

This is one of the most informative Stage 3 failure types. 

`invalid_umi_chars`

The extracted 10-base UMI contains characters other than A/C/G/T. 

`other_parse_fail`

A general fallback category for failures that do not fit a more specific rule. This label is referenced in the stage documentation and summary logic. 

### UMI family count data
`umi_dict_path`

Path to a serialized dictionary storing UMI family counts for a given population. This is consumed by Stage 3.5.

Interpretation:

- each UMI key represents one family

- the associated integer value represents the number of reads assigned to that family

## 7. Stage 3.5 — UMI family statistics

Stage 3.5 loads the per-population UMI count dictionaries and computes family-size summary metrics.

### Main Stage 3.5 metrics
`unique_families`

Number of distinct UMI sequences observed.

Interpretation:

- approximates molecular diversity

- larger is not always better; interpretation depends on library complexity and sequencing depth

`mean_reads_per_family`

Average number of reads per UMI family.

Formula:
```
mean_reads_per_family = sum(family sizes) / number of families
```

`median_family_size`

Median family size.

Interpretation:

more robust than the mean when a small number of UMI families are extremely large

`max_family_size`

Largest observed UMI family size.

Large values may indicate over-amplification, duplication, or highly abundant molecules.


### Family-size distribution bins

Stage 3.5 explicitly summarizes the following bins:

- n_size_1, pct_size_1

- n_size_2, pct_size_2

- n_size_3, pct_size_3

- n_size_gt3, pct_size_gt3 

`n_size_1`

Number of UMI families supported by exactly one read.

High values suggest low redundancy and limited support for consensus-based correction.

`n_size_2`

Number of UMI families supported by exactly two reads.

`n_size_3`

Number of UMI families supported by exactly three reads.

`n_size_gt3`

Number of UMI families supported by more than three reads.

`pct_*`

Corresponding percentage among all UMI families.

Formula pattern:
```
pct_size_k = n_size_k / unique_families
```

### Consensus-relevant support metrics
`n_ge2, pct_ge2`

Number and proportion of UMI families supported by at least two reads.

Interpretation:

- important because families with at least two reads are often the minimum for simple duplicate-aware consensus rules

`n_ge3, pct_ge3`

Number and proportion of UMI families supported by at least three reads.

Interpretation:

- useful when more conservative consensus support is desired 


## 8. Stage 4 — Clean insert QC metrics

Stage 4 performs QC after index / primer / UMI structures have been removed, so reads should now primarily contain clean insert sequence. It reuses the raw overview logic on the Stage 3 cleaned FASTQ outputs and also runs FastQC per population.

### Per-population metrics

The code explicitly stores, for each population:

- `r1_total_reads`

- `r2_total_reads`

- `r1_q30_rate`

- `r2_q30_rate`


These are interpreted exactly as in earlier raw-style overview stages.

Since Stage 4 uses `compute_raw_overview(..., short_len=150)`, the underlying summary files also implicitly contain:

- `short_rate`

- `length_counter`

even though the returned per-population dictionary shown in the snippet highlights only read count and Q30 fields.

Interpretation:

- low mean sequence length here may simply reflect expected insert size, depending on the library

- this stage is especially useful for confirming that structural sequence removal behaved as expected

## 9. Stage 5 — Interactive preprocessing metrics

Stage 5 is the most interactive stage. It evaluates multiple trim settings, lets the user choose a quality threshold, then previews length filtering before writing final trimmed outputs. 

9.1 Candidate Q evaluation metrics

Candidate trim quality thresholds:
```
30, 29, 28, 25, 20
```
These are tested on a sampled subset of read pairs. For each candidate, the following are recorded:

`Q30_after`

Base-level Q30 rate after running trim-only fastp on the sample. 

`trimmed_bases_pct`

Fraction of bases trimmed away.

Formula:
```
trimmed_bases_pct = trimmed_bases / total_bases_before
```
Interpretation:

- higher values indicate more aggressive trimming

- high Q30_after with very high trimmed_bases_pct may indicate over-trimming

`mean_len`

Mean read length after trim-only processing.

Interpretation:

- helps assess how much sequence length is preserved at a given quality cutoff

These metrics are extremely useful for selecting a quality threshold that balances quality improvement against read retention. 


### 9.2 Full trim-only metrics

After a user selects q_threshold, the pipeline performs full trim-only processing and records:

- `q30_before`

- `q30_after`

- paths to fastp_html and fastp_json

- trim-only FASTQ output paths per population 


`q30_before`

Base-level Q30 before trim-only processing, as parsed from the fastp JSON summary.

`q30_after`

Base-level Q30 after trim-only processing.

Interpretation:

- the difference between q30_before and q30_after quantifies the quality gain produced by trimming

### 9.3 Length filtering preview metrics

Stage 5 then previews filtering across multiple minimum-length cutoffs. The code uses a pair-level rule based on:
```
Lmin = min(length(R1), length(R2))
```
A pair is considered short for cutoff c if:
```
Lmin < c
```
This is stricter than read-level short-rate and is appropriate for paired-end downstream processing. 

For each cutoff, the preview table records:

`short_rate`

Fraction of pairs that would be filtered out because at least one mate is shorter than the selected cutoff.

Formula:
```
short_rate = short_pairs / total_pairs
```
`kept_rate`

Fraction of pairs that would remain after filtering.

Formula:
```
kept_rate = 1 - short_rate
```

`est_dropped_pairs`

Estimated number of pairs that would be removed under that minimum-length threshold.

`total_pairs`

Total number of pairs examined in the preview table. 

### Stage 5 length filtering warnings

```
short_rate > 30% → warning
short_rate > 40% → strong warning
short_rate > 50% → abnormal
```

Interpretation:

- if too many pairs would be dropped, the chosen trimming + filtering combination may be too aggressive

- the stage allows the user to go back and select a different trimming threshold before finalizing outputs 


## 10. Stage 6 — Post-trim QC metrics

Stage 6 evaluates the final trimmed reads and produces both per-population and overall summaries. It uses a different overview function from raw stages because short-read judgment is no longer the main concern here. Instead, the focus is on the final quality and length distribution of retained reads.

### Per-file post-trim metrics
`total_reads`

Number of retained reads in the post-trim FASTQ file. 

`total_bases`

Total number of retained bases.

`q30_bases`

Number of retained bases with Phred score at least 30.

This field is kept specifically so that population-level summaries can be aggregated exactly rather than approximately. 

`q30_rate`

Final base-level Q30 rate.

Formula:
```
q30_rate = q30_bases / total_bases
```

`mean_len`

Mean read length after all trimming and filtering.

`median_len`

Median read length after all trimming and filtering.

`length_counter`

Read length frequency map used to generate `read_length_frequency.png`. 

### Overall aggregated metrics

Stage 6 also aggregates metrics across all populations by summing:

- `total_reads`

- `total_bases`

- `q30_bases`

- `length_counter`

and then recomputing:

- `q30_rate`

- `mean_len`

- `median_len`


Interpretation:

- overall summaries give a global view of the final dataset

- because they are aggregated from counts rather than averaging per-population percentages, they preserve the correct weighting by population size

## 11. Output report files and how to read them
`raw_overview_summary.txt`

Written by `write_overview_report` in raw-style QC stages.

Contains:

- total reads

- short rate `<150 bp`

- Q30 rate (base-level) 

Use this report to quickly assess:

- whether a file is generally high quality

- whether many reads are unexpectedly short

`length_distribution.png`

Bar chart of read-length frequencies in raw-style QC stages.
Useful for spotting:

- truncated reads

- multimodal length distributions

- strong library heterogeneity 

`post_trim_overview_summary.txt`

Written by `write_post_trim_overview` in Stage 6.

Contains:

- total reads

- total bases

- Q30 rate

- mean read length

- median read length 


This is the most concise final-quality report for trimmed data.

`read_length_frequency.png`

Step-line plot of read-length distribution after trimming.
Useful for evaluating whether trimming and length filtering produced the expected final read-size profile. 

`demux_summary.txt`

Written in Stage 2.5.

Contains:

- per-population pair counts

- total pairs

- unmatched pairs

- unmatched rate 

Use this file to verify that demultiplexing proportions look biologically plausible.

`umi_extraction_summary.txt`

Written in Stage 3 for each population.

Contains:

- total read pairs

- successful extraction count

- failed extraction count

- unmatched rate

- pair-level fail reasons

- R1 fail reasons

- R2 fail reasons 

This is the main troubleshooting report for UMI parsing.

`UMI_quality_summary.txt`

Written in Stage 3.5 for each population.

Contains:

- total pairs

- extracted pairs

- unmatched pairs

- number of unique UMI families

- mean / median / max family size

- family-size category counts and percentages

- number and fraction of families with support ≥2 or ≥3 reads 

This is the main report for evaluating family support and future consensus-read potential.

## 12. Practical interpretation guide
### If Stage 1 Q30 is low

The raw sequencing quality may be poor. Expect trimming to help somewhat, but not completely rescue very poor datasets.

### If Stage 2 unmatched rate is high

Check:

- index definitions in the multiplex CSV

- whether R1/R2 orientation is correct

- whether the expected dual-index design matches the actual library

### If Stage 3 unmatched rate is high but usable unmatched rate is much lower

The main problem is probably short reads, not parsing logic.

If `missing_anchor_right` dominates Stage 3

Check:

- suffix definition in the primer CSV

- whether the reads fully span the expected structure

- whether primer/UMI arrangement matches the actual library

### If Stage 5 `trimmed_bases_pct` is high at strict Q thresholds

The chosen trim setting may be too aggressive.

If Stage 5 preview short-rate is high at your desired minimum length

Consider:

- reducing the minimum length threshold

- choosing a less aggressive trim quality threshold

- confirming whether the library insert is expected to be short

### If Stage 3.5 shows mostly family size 1

Consensus-based correction may be limited because most UMIs do not have redundant read support.