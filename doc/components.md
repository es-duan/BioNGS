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


# Script: generate folders to store population demultiplexed files
- input: csv with DNA index (forward and reserve) sequences and its corresponding population
- output: create a folder and fastq file for each population

# Script: sort NGS reads by DNA index
- input: NGS file
- output: 