#!/bin/bash

# Step 5: Align reads to reference genome using Bowtie2, SAMtools, and BCFtools
# Usage: bash align_reads.sh <experiment> <reference_index> <reference_fasta>
# Example: bash align_reads.sh example input_data/reference_docs/rpoB_index/rpoB_index input_data/reference_docs/rpoBsequence.fasta

# Check if required arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Error: Incorrect number of arguments."
    echo "Usage: bash align_reads.sh <experiment> <reference_index> <reference_fasta>"
    echo "Example: bash align_reads.sh example input_data/reference_docs/rpoB_index/rpoB_index input_data/reference_docs/rpoBsequence.fasta"
    exit 1
fi

# Parse command line arguments
EXPERIMENT=$1
REFERENCE_INDEX=$2
REFERENCE_FASTA=$3

# Define base directories
FASTQ_DIR="results/${EXPERIMENT}/alignment/fastq"
SAM_DIR="results/${EXPERIMENT}/alignment/sam"
BAM_DIR="results/${EXPERIMENT}/alignment/bam"
VCF_DIR="results/${EXPERIMENT}/alignment/vcf"

# Check if input directory exists
if [ ! -d "$FASTQ_DIR" ]; then
    echo "Error: Input directory $FASTQ_DIR does not exist."
    echo "Please run alignment_prep.py (step 4) first."
    exit 1
fi

# Check if reference index exists
if [ ! -f "${REFERENCE_INDEX}.1.bt2" ]; then
    echo "Error: Reference index ${REFERENCE_INDEX}.1.bt2 does not exist."
    echo "Please provide a valid Bowtie2 index."
    exit 1
fi

# Check if reference fasta exists
if [ ! -f "$REFERENCE_FASTA" ]; then
    echo "Error: Reference fasta file $REFERENCE_FASTA does not exist."
    exit 1
fi

# Create output directories if they don't exist
mkdir -p "$SAM_DIR"
mkdir -p "$BAM_DIR"
mkdir -p "$VCF_DIR"

echo "=========================================="
echo "Starting alignment pipeline for experiment: $EXPERIMENT"
echo "Reference index: $REFERENCE_INDEX"
echo "Reference fasta: $REFERENCE_FASTA"
echo "=========================================="
echo ""

# Initialize counters
total_files=0
successful_alignments=0
failed_alignments=0

# Loop through all population folders
for pop_dir in "$FASTQ_DIR"/*; do
    if [ -d "$pop_dir" ]; then
        population=$(basename "$pop_dir")
        echo "Processing population: $population"
        
        # Create population subdirectories
        mkdir -p "$SAM_DIR/$population"
        mkdir -p "$BAM_DIR/$population"
        mkdir -p "$VCF_DIR/$population"
        
        # Find all R1 fastq files in this population
        for r1_file in "$pop_dir"/*_R1.fastq; do
            if [ -f "$r1_file" ]; then
                # Extract UMI pair name (remove _R1.fastq extension)
                umi_name=$(basename "$r1_file" _R1.fastq)
                r2_file="$pop_dir/${umi_name}_R2.fastq"
                
                # Check if corresponding R2 file exists
                if [ ! -f "$r2_file" ]; then
                    echo "  Warning: Missing R2 file for $umi_name. Skipping."
                    ((failed_alignments++))
                    continue
                fi
                
                ((total_files++))
                echo "  Processing UMI pair: $umi_name"
                
                # Define output file paths
                sam_file="$SAM_DIR/$population/${umi_name}.sam"
                bam_file="$BAM_DIR/$population/${umi_name}.bam"
                vcf_file="$VCF_DIR/$population/${umi_name}.vcf"
                
                # Step 1: Run Bowtie2 alignment
                echo "    Running Bowtie2 alignment..."
                if bowtie2 -x "$REFERENCE_INDEX" \
                           -1 "$r1_file" \
                           -2 "$r2_file" \
                           -S "$sam_file" 2>&1 | grep -q "overall alignment rate"; then
                    
                    # Step 2: Convert SAM to BAM and sort
                    echo "    Converting SAM to sorted BAM..."
                    if samtools view -S -b "$sam_file" | samtools sort -o "$bam_file"; then
                        
                        # Step 3: Generate VCF using bcftools
                        echo "    Generating VCF file..."
                        if bcftools mpileup -f "$REFERENCE_FASTA" "$bam_file" | \
                           bcftools call --ploidy 1 -mv -Ov -o "$vcf_file"; then
                            echo "    ✓ Successfully processed $umi_name"
                            ((successful_alignments++))
                        else
                            echo "    ✗ Failed to generate VCF for $umi_name"
                            ((failed_alignments++))
                        fi
                    else
                        echo "    ✗ Failed to convert SAM to BAM for $umi_name"
                        ((failed_alignments++))
                    fi
                else
                    echo "    ✗ Failed to align $umi_name"
                    ((failed_alignments++))
                fi
                
                echo ""
            fi
        done
    fi
done

# Print summary statistics
echo "=========================================="
echo "Alignment pipeline completed!"
echo "=========================================="
echo "Total UMI pairs processed: $total_files"
echo "Successful alignments: $successful_alignments"
echo "Failed alignments: $failed_alignments"
echo ""
echo "Output locations:"
echo "  SAM files: $SAM_DIR/{population}/"
echo "  BAM files: $BAM_DIR/{population}/"
echo "  VCF files: $VCF_DIR/{population}/"
echo "=========================================="

# Exit with appropriate status code
if [ $failed_alignments -eq 0 ]; then
    exit 0
else
    exit 1
fi
