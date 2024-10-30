#!/usr/bin/env python

"""
Description: Clean and rename sequence IDs in .faa files, appending the filename as a prefix to each sequence ID.
Input directory: results/03_prodigal/all_faa_files/
Output directory: results/03_prodigal/cleaned_faa_files/
Author: Manisha Chaudhary
Date: 10/24/2024
"""

import glob
import os
from Bio import SeqIO

# Define input and output directories
input_dir = "results/03_prodigal/all_faa_files/"
output_dir = "results/03_prodigal/cleaned_faa_files/"

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Process each .faa file in the input directory
for fname in glob.glob(input_dir + "*.faa"):
    # Extract the base name of the file without the extension for use in the sequence ID
    base_name = os.path.basename(fname)[:-4]

    # Open the output file for writing the cleaned sequences
    with open(os.path.join(output_dir, base_name + "_clean.faa"), "w") as output_handle:
        # Parse each sequence in the .faa file
        for seq_record in SeqIO.parse(fname, "fasta"):
            # Modify the sequence ID: append filename prefix to original ID
            seq_record.id = base_name + "_" + seq_record.id.split()[0]
            seq_record.description = ""  # Remove additional description

            # Write the modified sequence to the new file
            SeqIO.write(seq_record, output_handle, "fasta")

    print(f"Completed cleaning {fname}")

print("All .faa files cleaned and saved to the output directory.")
