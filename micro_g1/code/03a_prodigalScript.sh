#!/bin/bash

# Description: Performs PRODIGAL analysis on a list of genome sequences
# Working directory: /projects/intro2gds/I2GDS2024/individual_folders/manisha/
# Input files: data/assembly_data/ncbi_dataset/data/
# Output directory: results/03a_prodigal
# Author: Manisha Chaudhary
# Date: 09/20/2024

# Directories setup
INPUT_DIR="data/assembly_data/ncbi_dataset/data/"
OUTPUT_DIR="results/03_prodigal"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Find all .fna files in the directory and loop through them
find "$INPUT_DIR" -name "*.fna" | while read -r fna_file; do

    # Extract the directory path containing the file (for the sample name)
    sample_dir=$(dirname "$fna_file" | xargs basename)

    # Create the sample subdirectory in OUTPUT_DIR
    sample_output_dir="${OUTPUT_DIR}/${sample_dir}"
    mkdir -p "$sample_output_dir"

    # Extract the filename (without directory path)
    filename=$(basename "$fna_file")

    # Run PRODIGAL for each .fna file and output results to the output directory
    echo "Running PRODIGAL analysis for $filename..."
    prodigal -i "$fna_file" \
             -o "$sample_output_dir/${filename%.fna}_genes.gbk" \
             -a "$sample_output_dir/${filename%.fna}_proteins.faa" \
             -d "$sample_output_dir/${filename%.fna}_nucleotide.ffn" \
             -s "$sample_output_dir/${filename%.fna}_scores.txt"

    # Check if PRODIGAL ran successfully
    if [ $? -eq 0 ]; then
        echo "PRODIGAL analysis completed successfully for $filename."
    else
        echo "PRODIGAL analysis failed for $filename."
    fi
done
