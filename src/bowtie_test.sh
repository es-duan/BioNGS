#mkdir results/example/alignment/sam/
bowtie2 -x input_data/reference_docs/rpoB_index/rpoB_index -1 results/example/alignment/fastq/P1/AAATGCGTGA_GTCGGACCTC_R1.fastq -2 results/example/alignment/fastq/P1/AAATGCGTGA_GTCGGACCTC_R2.fastq -S results/example/alignment/sam/AAATGCGTGA_GTCGGACCTC.sam
#mkdir results/example/alignment/bam/
samtools view -S -b results/example/alignment/sam/AAATGCGTGA_GTCGGACCTC.sam | samtools sort -o results/example/alignment/bam/AAATGCGTGA_GTCGGACCTC.bam
#mkdir results/example/alignment/vcf/
bcftools mpileup -f input_data/reference_docs/rpoBsequence.fasta results/example/alignment/bam/AAATGCGTGA_GTCGGACCTC.bam | bcftools call --ploidy 1 -mv -Ov -o results/example/alignment/vcf/AAATGCGTGA_GTCGGACCTC.vcf