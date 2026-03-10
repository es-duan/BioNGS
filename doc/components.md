# Script 0: check quality of the NGS file
- name: check_fastq_quality.py
- input: input fastq files (R1.fastq and R2.fastq file)
- output: csvs and histograms with information on how many reads per NGS sample, the length of those reads, and the quality of reads
- description: use fastQC to quality check each input fastq file
- dependencies: only need input fastq files

# Script 1: generate folders to store population demultiplexed files
- name: demultiplex_folders.py
- input: csv with population name and DNA index (forward and reserve) sequences. 
- output: create a folder and empty R1 and R2 fastq files for each population
- description: this is the first script of the pipeline. the user should upload a csv file (see example_multiplexing_info.csv) with information on the populations from the experiment and what indexes were used with each population.

# Script 2: sort NGS reads by DNA index
name: demultiplex_index.py
- input: input fastq files, csv with DNA indexes
- output: populate the empty fastq files with the reads that belong to that population, based on matched index sequences; create an additional fastq file with reads that are too short (under 150 base pairs) and another fastq to store reads that do not match any of the index pairs
- dependencies: script to be run after folders are generated (script 1)
- description: using biopython, loop through each pair of input fastq files (ensure R1 and R2 match by checking the header), detect the forward (first 8 base pairs of R1) and reverse (first 8 bps of R2) indexes, and match those with the populations listed in the multiplexing_info.csv. For each match, open the population fastq file and write in that read. Also, make sure that the input fastq file name matches the GW_name of the population.
    - quality check sequences: to make the script faster, load R1 reads first and check the length of the read. If it is under 150 bps, it is a short read. Place both reads in a short_reads_R1 or R2.fastq. Additionally, if the F and R indexes do not match any populations in the csv, put the R1 and R2 reads in a unmatched_reads_R1 and R2.fastq. There should be one pair of short_read/unmatched_read files per input file

# Script 3: check quality of indexing/NGS data
name: check_index_quality.py
- input: input fastq files, populated demultiplexed population fastq files, short_read fastq files, and unmatched_read fastq files
- output: histogram of read lengths from input fastq files; bar plot with the number of reads in each population fastq, short_read, and unmatched_read files (one per input fastq)
- dependencies: script 2 must be run

# Script 4: create a dictionary of UMI reads for each population
- name: demultiplex_UMI.py
- input: demultiplexed population fastq files, UMI primer sequences
- output: dictionary where keys are UMI pairs and values are a list of a list of  R1 and R2 trimmed sequences that match that UMI. Save these dictionaries in the respective demultiplexing/population folders, named P{population number}_UMI_dict
- dependencies: script 2
- description: Loop through the demultiplexed fastq file for each population. Detect the UMI sequence by aligning to the primer sequence (the UMI is the 10 N bps on the primer) to each read. Forward UMIs are detected from R1, and reverse UMIs from R2. For each unique forward and reverse UMI pair, create a new key in a library. Store the sequences as under their respective UMI pair key. There should be one list of R1 reads and a corresponding list of R2 reads for each UMI. The R1 and R2 lists should be the same length and in order (ids should match).

# Script 5: check quality of UMI data
- name: quality_check_UMI.py
- input: dictionaries of UMI/sequences for each population
- output: bar plot comparing the number of UMIs/population and reads per UMI per population
- dependences: script 4
- description: compare the number of UMI pairs per population (length of dictionary), compare the number of sequences per UMI per population (length of R1 and R2 lists)

# Script 6: prepare UMI sequences for alignment to the genome
- name: alignment_prep.py
- input: UMI dictionary from step 4, UMI primer sequences
- output: in a new folder under results/{experiment}/alignment/fastq, save R1 and R2 files for each UMI pair with more than 2 reads. The names should be {forwardUMI_reverseUMI}_R1.fastq and {forwardUMI_reverseUMI}_R2.fastq. Also include terminal output and csv version of terminal output with quality information that can be easily accessed in dataframe format.
- description: Open the UMI dictionaries from script 4. Filter the libraries to remove UMI keys with only one forward and reverse read value. Trim the sequences to remove the entire primer sequence (detect matches to primer_after and set the last position as the trim position). Report the total number of keys, the number of keys that will move forward with alignment (greater than or equal to the minimum sequence pair), the number of keys removed (less than minimum) + percentage. Then, for the keys with more than or equal to the minimum, save the sequences in the output dir as R1 and R2 fastq files.


# Script 7: use the bowtie aligner and SAMtools to align reads to the genome
- name: align_reads.sh
- input: UMI fastqs, saved from script 6 in results/{experiment}/alignment/fastq, rpoB reference genome index (prepared with bowtie in input_data/reference_docs)
- output: sam files in vcf outputs in results/{experiment}/alignment/sam, bam outputs in results/{experiment}/alignment/bam, vcf outputs in results/{experiment}/alignment/vcf, with population folders, and named by UMI pair (same as fastq)
- dependencies: script 4
- description: use the bowtie_test.sh tool as a template. For each UMI pair in the alignment/fastq file, use bowtie2 to align the reads to the specified reference genome index. Use samtools to convert the outputted sam file to a bam file, then bcftools to convert the bam files to vcf files. Have an option to include command arguments: experiment (ex. example, this is consistent with previous python scripts), index (specify the reference genome index to use for alignment)

# Script 7: 