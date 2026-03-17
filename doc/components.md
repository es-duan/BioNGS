# Step 0: check presence of files
- name: step0_file_confirmation.py
- input: name of the experiment (ex. example). the user should have uploaded a folder in the input_data folder names with the experiment and containing various files for analysis
- output: txt file/terminal output with information on what files the user has uploaded, to what step in the pipeline they are able to run, and suggestions if other files are uploaded and potentially named incorrectly
- description: The script uses a command line input of the experiment name. Look in the input_data folder for a folder named {experiment}. Within that folder, there should be a folder named {experiment}_fastq which contains at least 1 pair of fastq files (R1 and R2). Additionally, the input_data/{experiment} folder should contain a csv named {experiment}_multiplexing_info.csv and a csv files named {experiment}_UMI_primers.csv. After files are confirmed present or absent, determine what steps can be run. For example, if only fastq files are provided, only step 0a can be run. If there are additional files that do not match naming patterns, determine if there may be a typo in their names based on how similar they are to the specified patterns. Print out a report of files present and absent, the furthest step down the pipeline that can be run, and the names of any other files that may have been named incorrectly.
- dependencies: no other scripts are necessary

# Step 0a: check the quality of NGS data with FastQC
- name: step0a_fast_qc.py
- inputs: experiment name; fastq files in `input_data/{experiment}/{experiment}_fastq/`; `fastqc` installed and available in PATH
- outputs: HTML FastQC reports in `results/{experiment}/qc/` named `{fastq_basename}_fastqc.html`; non-HTML FastQC artifacts are removed
- dependencies: none (Step 0 recommended)
- description: finds all raw fastq/fq inputs (including gzipped), runs FastQC, then keeps only expected HTML reports
- errors (what errors are raised): `FileNotFoundError` if input fastq directory/files are missing; `RuntimeError` if FastQC is not in PATH or no HTML reports are produced; `subprocess.CalledProcessError` if FastQC execution fails

# Step 1: generate folders for demultiplexed populations
- name: step1_demultiplex_folders.py
- inputs: experiment name; multiplexing CSV in `input_data/{experiment}/` matching `*multiplexing_info*.csv`
- outputs: `results/{experiment}/demultiplexing/P{population}/` with empty `P{population}_R1.fastq` and `P{population}_R2.fastq`; terminal log in `results/{experiment}/logs/demultiplex_folders_terminal_output.txt`
- dependencies: Step 0 recommended
- description: reads population rows from the multiplexing CSV and creates one output folder plus empty paired fastq files per population
- errors (what errors are raised): `FileNotFoundError` if experiment directory or multiplexing CSV is missing

# Step 2: sort NGS reads by DNA index
- name: step2_demultiplex_index.py
- inputs: experiment name; multiplexing CSV; input R1/R2 fastq files (plain or gzipped); optional `--short-read-length` threshold (default 150)
- outputs: populated `results/{experiment}/demultiplexing/P{population}/P{population}_R1.fastq` and `_R2.fastq`; `{GW_name}_short_reads_R1.fastq` and `_R2.fastq`; `{GW_name}_unmatched_reads_R1.fastq` and `_R2.fastq`; terminal log in `results/{experiment}/logs/demultiplex_index_terminal_output.txt`
- dependencies: Step 1
- description: pairs R1/R2 reads by GW_name, validates matching read IDs, extracts index (first 8 bases) and matches it to the multiplexing csv, classifies reads as short, matched-to-population-index, or unmatched, then writes to corresponding output files
- errors (what errors are raised): `FileNotFoundError` if multiplexing CSV is missing; if output directory is missing the script prints an error and returns; CLI argument parsing errors are handled by argparse

# Step 2a: check quality of indexing/NGS data
- name: step2a_check_index_quality.py
- inputs: experiment name; input fastq pairs (plain or gzipped); demultiplexing outputs from Step 2 (population fastqs, short reads, unmatched reads)
- outputs: `results/{experiment}/index_quality/index_quality_results.csv`; `index_quality_summary.txt`; `{experiment}_read_length_histogram.png` (or HTML fallback); `{experiment}_read_distribution.png` (or HTML fallback)
- dependencies: Step 2
- description: counts reads per population/QC category, computes recovery metrics, samples read lengths from raw inputs, and generates QC plots plus tidy CSV and summary text
- errors (what errors are raised): `FileNotFoundError` if multiplexing CSV or demultiplexing output directory is missing; PNG export failures are caught and downgraded to HTML fallback (warning only)

# Step 3: create UMI dictionaries per population
- name: step3_demultiplex_UMI.py
- inputs: experiment name; demultiplexed population fastq files; multiplexing CSV; UMI primers CSV containing `f` and `r` columns with exactly 10 consecutive `N` bases in each primer
- outputs: per-population `P{population}_UMI_dict.pkl` under `results/{experiment}/demultiplexing/P{population}/`; unmatched UMI fastq files; terminal log in `results/{experiment}/logs/demultiplex_UMI_terminal_output.txt`
- dependencies: Step 2
- description: extracts forward/reverse UMIs from read pairs using primer structure, stores matched reads as UMI-pair keyed dictionaries with trim positions for downstream alignment prep
- errors (what errors are raised): `ValueError` if primer format is invalid (not exactly 10 N bases) or required primer columns are missing; `FileNotFoundError` if multiplexing CSV, primer CSV, or required input paths are missing; some missing-path conditions are printed as errors and the script returns

# Step 3a: check quality of UMI data
- name: step3a_check_UMI_quality.py
- inputs: experiment name; UMI dictionaries `*_UMI_dict.pkl` in `results/{experiment}/demultiplexing/P*/`; unmatched UMI fastq files `P*_unmatched_UMI_R1.fastq` in `results/{experiment}/demultiplexing/P*/`
- outputs: `results/{experiment}/UMI_quality/UMI_count_per_population.png`; `reads_per_UMI_distribution.png`; `UMI_match_status.png` (faceted bar chart of matched vs unmatched UMI counts per population); `UMI_quality_summary.txt`; terminal metrics including UMI count, total reads, min/max reads per UMI
- dependencies: Step 3
- description: loads all UMI libraries, computes per-population UMI/read statistics, counts unmatched UMI reads from output fastq files, and generates population-level comparison plots (UMI counts, read distribution, and match status) plus text summary
- errors (what errors are raised): `FileNotFoundError` if demultiplexing directory is missing (caught in main quality check and reported); if no libraries are found/loadable, script prints an error and returns failure (`False`); warnings logged if unmatched fastq files cannot be found or counted; plot export exceptions are caught and printed as warnings

# Step 4: prepare UMI sequences for alignment
- name: step4_alignment_prep.py
- inputs: experiment name; UMI dictionaries from Step 3; optional `--min-reads-per-umi` threshold
- outputs: per-population alignment fastq files in `results/{experiment}/alignment/fastq/P{population}/` named `{forwardUMI}_{reverseUMI}_R1.fastq` and `_R2.fastq`; terminal log in `results/{experiment}/logs/alignment_prep_terminal_output.txt`
- dependencies: Step 3
- description: loads UMI dictionaries, trims reads using stored trim positions, filters by minimum reads-per-UMI, and writes trimmed paired fastq files for downstream alignment
- errors (what errors are raised): `FileNotFoundError` if demultiplexing directory is missing (caught and reported); if no UMI libraries are found, script prints an error and returns failure (`False`)

# Step 5: align reads to reference and call variants
- name: step5_align_reads.sh
- inputs: `experiment`, `reference_index`, and `reference_fasta` CLI args; prepared alignment fastqs in `results/{experiment}/alignment/fastq/P*/`; bowtie2/samtools/bcftools installed
- outputs: per-population SAM/BAM/VCF in `results/{experiment}/alignment/sam/{population}/`, `results/{experiment}/alignment/bam/{population}/`, and `results/{experiment}/alignment/vcf/{population}/`; terminal summary of success/failure counts
- dependencies: Step 4
- description: iterates population folders and UMI pairs, runs Bowtie2 alignment, converts SAM→sorted BAM with SAMtools, then runs BCFtools mpileup/call to produce VCF per UMI pair
- errors (what errors are raised): shell script exits with code `1` for invalid arg count, missing input directories/files, or any failed alignment/conversion/calling task; missing R2 for a UMI pair is logged as a warning and counted as failed
