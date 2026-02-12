# Script 1: generate folders to store population demultiplexed files
- input: csv with population name and DNA index (forward and reserve) sequences
- output: create a folder and empty R1 and R2 fastq files for each population
- description: this is the first script of the pipeline. the user should upload a csv file (see example_multiplexing_info.csv) with information on the populations from the experiment and what indexes were used with each population.

# Script 2: sort NGS reads by DNA index
- input: input fastq files, csv with DNA indexes
- output: populate the empty fastq files with the reads that belong to that population, based on matched index sequences; create an additional fastq file with reads that are too short (under 150 base pairs) and another fastq to store reads that do not match any of the index pairs
- dependencies: script to be run after folders are generated (script 1)
- description: using biopython, loop through each pair of input fastq files (ensure R1 and R2 match by checking the header), detect the forward (first 8 base pairs of R1) and reverse (first 8 bps of R2) indexes, and match those with the populations listed in the multiplexing_info.csv. For each match, open the population fastq file and write in that read. Also, make sure that the input fastq file name matches the GW_name of the population.
    - quality check sequences: to make the script faster, load R1 reads first and check the length of the read. If it is under 150 bps, it is a short read. Place both reads in a short_reads_R1 or R2.fastq. Additionally, if the F and R indexes do not match any populations in the csv, put the R1 and R2 reads in a unmatched_reads_R1 and R2.fastq. There should be one pair of short_read/unmatched_read files per input file

# Script 3: check quality of NGS data
- input: input fastq files, populated demultiplexed population fastq files, short_read fastq files, and unmatched_read fastq files
- output: histogram of read lengths from input fastq files; bar plot with the number of reads in each population fastq, short_read, and unmatched_read files (one per input fastq)
- dependencies: script 2 must be run