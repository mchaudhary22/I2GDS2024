#!/bin/bash

# Description: Performs QUAST analysis on a list of genome sequences 
# working directory: /projects/intro2gds/I2GDS2024/individual_folders/manisha/
# input files: data/assembly_data/ncbi_dataset/data/
# output directory: results/02_QUAST
# Author: Manisha Chaudhary
# Date: 09/19/2024


# Directories setup
INPUT_DIR="data/assembly_data/ncbi_dataset/data/"
OUTPUT_DIR="results/02_QUAST"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Find all .fna files in the directory and loop through them
find "$INPUT_DIR" -name "*.fna" | while read -r fna_file; do
    # Extract the filename (without directory path)
    filename=$(basename "$fna_file")

    # Run QUAST for each .fna file and output results to a specific folder
    echo "Running QUAST analysis for $filename..."
    quast "$fna_file" -o "$OUTPUT_DIR/${filename%.fna}_quast_output"

    # Check if QUAST ran successfully
    if [ $? -eq 0 ]; then
        echo "QUAST analysis completed successfully for $filename."
    else
        echo "QUAST analysis failed for $filename."
    fi
done

# Notify when all QUAST analyses are done
echo "QUAST analysis completed for all GCA files."
