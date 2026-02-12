|   Component   |   Description |
|---------------|---------------|
|1. Function to generate folder | A function generates an empty folder |
|2. Function to read large NGS fastq file | Need the appropriate function to read the type of file provided by user (text file?CSV?)
|3. A function to SPLIT text based on indexes | Each bacteria has a primer associated with their DNA to idenitfy them with their population. This is the first node in the sorting process |
|4. Identifying indexes and SORTING into respective folders | Along with identifying indexes to respective populations, this step may also expose errors in DNA scripts and therefore these would be sorted into a separate folder as well as each population being sorted separately | 
|5. Function to further separate populations by UMI | After sorting to appropriate populations via 8 base index, further categorization is needed based on a UMI index to know which sequence to compare against for determining mutation rates |
|6. Next a function to trim sequences (what criteria again?) | After sorting into appropriate UMI folders, sequences need to be trimmed to prep for comparison |
|7. Function to align trimmed genes for comparison and sort into respective folder | Finally, DNA sequences should be prepped for comparison by being aligned |
|8. Calculations for mutation rates | based on comparisons, execute calculations to determine mutation rates |
|9. Function to convert calculations into a colorful convenient visualization | Ideally, a figure for making the data more digestable is an end goal |



# Script 0: check quality of the NGS file
- input: input fastq files
- output: csvs and histograms with information on how many reads per NGS sample, the length of those reads, and the quality of reads
- description: use fastQC to quality check each input fastq file
- dependencies: only need input fastq files

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

# Script 3: check quality of indexing/NGS data
- input: input fastq files, populated demultiplexed population fastq files, short_read fastq files, and unmatched_read fastq files
- output: histogram of read lengths from input fastq files; bar plot with the number of reads in each population fastq, short_read, and unmatched_read files (one per input fastq)
- dependencies: script 2 must be run

# Script 4: create a library of UMI reads for each population
- input: demultiplexed population fastq files, UMI primer sequences
- output: library where keys are UMI pairs and values are a list of a list of  R1 and R2 trimmed sequences that match that UMI
- dependencies: script 2
- description: take the demultiplexed fastq file for each population. Detect the UMI sequence be aligning to the primer sequence (the UMI is the 10 N bps on the primer) to each read. For each unique UMI, create a new key in a library. Then, loop though the the R1 and R2 files, detect UMIs, trim the primer (including index and UMI) from each read and place it in the respective R1 and R2 list for its UMI. Save the outputted list

# Script 5: check quality of UMI data
- input: libraries of UMI/sequences for each population
- output: bar plot comparing the number of UMIs/population and reads per UMI per population
- dependences: script 4
- description: compare the number of UMI pairs per population (length of library), compare the number of sequences per UMI per population (length of R1 and R2 lists)

# Script 6: align sequences to reference gene
- input: libraries of UMI/sequences, reference gene.gb