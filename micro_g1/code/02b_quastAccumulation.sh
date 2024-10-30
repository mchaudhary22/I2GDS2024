#!/bin/bash

# Description: Accumulates all QUAST report.txt files into a single directory
# working directory: /projects/intro2gds/I2GDS2024/individual_folders/manisha/
# input files: results/02_QUAST
# output directory: results/02_QUAST_reports
# Author: Manisha Chaudhary
# Date: 09/20/2024

# Define directories
QUAST_OUTPUT_DIR="results/02_QUAST"
REPORTS_DIR="results/02_QUAST_reports"

# Create the reports directory
mkdir -p "$REPORTS_DIR"

# Find all report.txt files in the QUAST output directory and copy them
find "$QUAST_OUTPUT_DIR" -name "report.txt" | while read -r report_file; do
    # Extract the folder name of the report (assuming folder structure per assembly)
    assembly_name=$(basename "$(dirname "$report_file")")
    
    # Copy the report file to the reports directory with a unique name
    cp "$report_file" "$REPORTS_DIR/${assembly_name}_report.txt"

    # Check if the copy was successful
    if [ $? -eq 0 ]; then
        echo "Copied report.txt for $assembly_name successfully."
    else
        echo "Failed to copy report.txt for $assembly_name."
    fi
done

echo "All report.txt files accumulated in $REPORTS_DIR."
