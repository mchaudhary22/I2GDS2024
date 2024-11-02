# Introduction to Microbiome Pipeline Tools
This pipeline employs three key tools — QUAST, Prodigal, and MMSeq2 — for evaluating genome assembly quality, predicting gene sequences, and clustering gene families. These steps are essential in genome/microbiome studies, especially for accurate assembly assessment, gene identification, and comparative analysis across genomes.

This repo was developed as part of the curriculum for Virginia Tech’s Introduction to Genomic Data Science course. It provides example scripts for performing genomic/microbiome analyses, though some paths may need to be adjusted depending on your installation method or file locations. Please ensure your script paths, data files, and working directory are correctly specified for smooth execution.

For questions or assistance, feel free to reach out to:
Manisha Chaudhary (mchaudhary@vt.edu), Lena Patino Westermann (patinole@vt.edu) or Jiayu Dong (jiayu24@vt.edu)

# Dependencies & Version Information
Ensure the following dependencies are installed to run the pipeline:

- **QUAST**: v5.0.2
- **Prodigal**: V2.6.3
- **MMSeqs2**: 13.45111
- **Additional Python Packages**: biopython 1.84, numpy 1.26.4, pandas 2.2.2, xlrd 2.0.1, openpyxl 3.1.5

# Environment Setup
Load Anaconda3 on ARC to create a conda environment and activate it before running each tool:

    module load Anaconda3
    conda create -n <env name>
    source activate <env name>

After creating and activating the environment, install each tool using conda.

# Data download

Install the **NCBI CLI tool** using Conda

```
conda install -c conda-forge ncbi-datasets-cli
```

Download genomic data by **NCBI CLI tool**
```
datasets download genome accession GCA_041080895.1 --output-dir/path/to/save
```

**Note**: Use the following command to download the provided genome assemblies. The “datasets” CLI tool can be used for downloading the data when the assembly accessions are provided. The command script is in "01_data_download.sh", and the list of genomes is specified in "Genome_list.txt".

<details>
    <summary>01_data_download.sh</summary>

```
#!/bin/bash

# Description: Download genomic data from ncbi website, using the Assembly IDs in Genome_List.txt 
# working directory: /projects/intro2gds/I2GDS2024/individual_folders/manisha
# input files: raw_data/Genome_List.txt
# output directory: data/assembly_data
# Author: Manisha Chaudhary
# Date: 09/16/2024

#SBATCH -J Data_download # job name
#SBATCH --account= ENTER_ACCOUNT
#SBATCH --partition=normal_q
#SBATCH --time=1-00:00:00 # 10 minutes; format: days-hours:minutes:seconds
#SBATCH --mem=128G # based on memory used when testing with interactive jobs
#SBATCH --mail-user= EMAIL #enter desired email address for updates
#SBATCH --mail-type=BEGIN #include to get emailed when job begins
#SBATCH --mail-type=END #include to get emailed when job ends
#SBATCH --mail-type=FAIL #include to get emailed if job fails

# Define the file containing the list of accession codes
ACCESSION_FILE="raw_data/Genome_List.txt"
DOWNLOAD_DIR="data/assembly_data"

mkdir -p "$DOWNLOAD_DIR"

#Use the datasets command to download all genomes listed in the file
datasets download genome accession --inputfile "$ACCESSION_FILE"

#Move downloaded files to the specified directory
#assuming files are named based on the accession codes
	
    for accession in $(cat "$ACCESSION_FILE") ; do
	echo "Moving files for accession: $accession"

#construct the expected filename
filename= "${accession}.zip"
#Move using the constructed filename
    mv "${accession}_*.zip" "$DOWNLOAD_DIR/$filename"
    
 #check if the command was successful
     if [ $? -eq 0 ] ; then
      	echo "Downloaded $accession successfully."
     else
	echo "Failed to download $accession."
     fi
done < "$ACCESSION_FILE"

echo "Download completed."

```

</details>
<p></p>

After downloading the data, unzip the file using the following command:

	unzip filename.zip

# 1. QUAST Analysis (Quality Assessment of Genome Assemblies)

### 1.1 Introduction
[QUAST](https://bioinf.spbau.ru/quast) assesses the quality of genome assemblies by analyzing contigs for metrics like N50, L50, and total contig length. This provides insights into the completeness and accuracy of microbial genome assemblies, a key step before further analysis.

## 1.2 Installation
To install QUAST on ARC:

    conda install -c bioconda quast

After installation, user can verify the version:

    quast --version

## 1.3 Running QUAST
QUAST accepts assembly files and outputs detailed reports on genome quality:

    quast /path/to /downloaded_assembled_data/*.fna -o quast_output

**Note**: If a user needs to run QUAST on multiple .fna files (all starting with "GCA"), you can automate this in the script. After saving the script, make it executable by running. Refer to this bash script as “02a_quastScript.sh”

<details>
    <summary>02a_quastScript.sh</summary>

```
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

```
</details>

## 1.4 Accumulate all the reports of the QUAST results to a single directory using a bash script

Refer this bash script for accumulation of .txt files from QUAST output in separate folder 
"02b_quastAccumulation.sh"

<details>
    <summary>02b_quastAccumulation.sh</summary>

```
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

```

</details>

## 1.5 Input Requirements
- **Input Format**: FASTA files of genome assemblies, with a minimum contig length threshold specified to filter out short contigs.
- **Input Location**: Store assembly files in a dedicated directory, such as assembly_data/

## 1.6 Output Explanation
- Reports:
    - N50 and L50: Metrics for assembly contiguity.
    - Total Contig Length: The total length of all contigs.
    - GC Content: An indicator of the genomic composition.
- Files Generated:
    - report.txt: Summary report of all assemblies.
    - contigs_reports: Folder with detailed contig-based metrics.

**GitHub Link**

For more information, visit the [QUAST GitHub page](https://github.com/ablab/quast).

# 2. Prodigal Analysis (Gene Prediction and Annotation)

## 2.1 Introduction

PRODIGAL predicts protein-coding genes within prokaryotic genomes, identifying open reading frames and translating initiation sites. This is critical for annotating microbial genomes and inferring gene functions.

## 2.2 Installation
To install Prodigal on ARC:

    conda install -c bioconda prodigal

Verify installation with:

    prodigal -v

## 2.3 Running Prodigal
Run Prodigal on each genome assembly to predict genes:

    prodigal -i genome.fasta -a genes.faa -d genes.ffn -o output.gbk -s scores.txt

**Note**: User can create and run the bash script for analyzing your multiple .fna files, make the script executable and run it. Refer this bash script as "03a_prodigalScript.sh"

<details>
    <summary>03a_prodigalScript.sh</summary>

```
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

```
</details>

## 2.4 Input Requirements
- Input Format: Genomic FASTA files from which genes are to be predicted.
- Input Location: Place all input genomes in a directory, e.g., genome_data/.

## 2.5 Output Explanation
- Output Files:
    - .gbk file: Contains gene prediction results in GenBank format.
    - .faa file: FASTA file of predicted protein sequences.
    - .ffn file: FASTA file of predicted nucleotide sequences of genes.
    - .txt file: Score file with translation tables and gene scores.

**GitHub Link**

For more details, refer to the [Prodigal GitHub page](https://github.com/hyattpd/Prodigal).

# 3. MMSeq2 Analysis (Clustering of Protein Sequences)

## 3.1 Introduction
MMSeqs2 performs many-against-many sequence comparisons and clusters sequences based on similarity. This tool is used to create clusters of homologous gene families across genomes, essential for microbial comparative studies.

## 3.2 Installation
To install MMSeqs2 on ARC:

    conda install -c bioconda mmseqs2

Verify installation:

    mmseqs --version

## 3.3 Running MMSeq2
### 3.3.1 Prepare Input Data
You have accumulated all the .faa files from Prodigal analysis under /results/03_prodigal/all_faa_files. These files contain the predicted protein sequences and will be the input for MMseqs2.
- Ensure that your .faa files are in the correct format and named appropriately for each genome. Refer this bash script for accumulation of all .faa files “03b_faaAccumulation.sh”

<details>
    <summary>03b_faaAccumulation.sh</summary>

```
# Description: Collects all .faa files generated by Prodigal analysis
# Working directory: /projects/intro2gds/I2GDS2024/individual_folders/manisha/
# Input directories: results/03_prodigal/sample_subdirectories/
# Output directory: results/03_prodigal/all_faa_files/
# Author: Manisha Chaudhary
# Date: 09/28/2024

# Set the directory where the Prodigal results are stored
INPUT_DIR="results/03_prodigal"
FAA_DIR="$INPUT_DIR/all_faa_files"

# Create a directory for the .faa files if it doesn't exist
mkdir -p "$FAA_DIR"

# Find all .faa files in the subdirectories of INPUT_DIR and move them to FAA_DIR
find "$INPUT_DIR" -type f -name "*.faa" | while read -r faa_file; do
    # Move each .faa file to the all_faa_files directory
    echo "Moving $faa_file to $FAA_DIR"
    cp "$faa_file" "$FAA_DIR/"
done

echo "All .faa files have been moved to $FAA_DIR."

```
</details>

<p></p>
Modify the header of the FASTA files to shorten them and change the ID to just the sequence of the genome. For example, a previous header

"CP074663.1_1 # 1 # 1401 # 1 # ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.525" will now become "GCA_001756025_CP074663.1_1"
- The first 13 characters of this ID signify the genome name and remaining identify the gene ID. Use the following code to change the header document “03c_clean_faa_files.py”

<details>
    <summary>03c_clean_faa_files.py</summary>

```
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
```

</details>

### 3.3.2 Create a Database with MMseqs2
MMseqs2 requires the creation of a sequence database. Use the following command to create the database from your .faa files:

    mmseqs createdb /results/03_prodigal/all_faa_files/*.faa /results/mmseqs2_db

- This command will create a database from all the .faa files.

### 3.3.3 Perform Sequence Clustering with MMseqs2
Now, perform sequence clustering to identify homologous gene families across genomes. Run the following command:

    mmseqs cluster /results/mmseqs2_db /results/mmseqs2_cluster /results/tmp --min-seq-id 0.5 -c 0.8

- --min-seq-id 0.5: Sets a minimum sequence identity threshold of 50% for clustering.
- -c 0.8: Sets a coverage threshold, ensuring that 80% of the sequence is aligned.

### 3.3.4 Extract Cluster Information
Once clustering is complete, extract the cluster information to analyze gene families:

```
mmseqs createtsv /results/mmseqs2_db /results/mmseqs2_db /results/mmseqs2_cluster /results/mmseqs2_cluster.tsv
```

This will generate a TSV file with the clustering results, showing which genes belong to which gene families.

**Note**: Use the following script to run MMSeq2 analysis "04_mmseqScript.sh"

<details>
    <summary>04_mmseqScript.sh</summary>

```
# Description: Performs MMSeq2 analysis on accumulated .faa files obtained from prodial analysis
# Working directory: /projects/intro2gds/I2GDS2024/individual_folders/manisha/
# Input files: results/03_prodigal/cleaned_faa_files/
# Output directory: results/04_mmseq
# Author: Manisha Chaudhary
# Date: 09/28/2024

# Set directories
FAA_DIR="results/03_prodigal/cleaned_faa_files/"
OUTPUT_DIR="results/04_mmseq"
MMSEQS_DB="$OUTPUT_DIR/mmseqs2_db_cleaned"
MMSEQS_CLUSTER="$OUTPUT_DIR/mmseqs2_cluster_cleaned"
MMSEQS_TMP="$OUTPUT_DIR/tmp"
TSV_OUTPUT="$OUTPUT_DIR/mmseqs2_cluster_final.tsv"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check if the faa files directory exists
if [ ! -d "$FAA_DIR" ]; then
    echo "FAA directory not found: $FAA_DIR"
    exit 1
fi

# Create MMseqs2 database
echo "Creating MMseqs2 database..."
mmseqs createdb "$FAA_DIR"/*.faa "$MMSEQS_DB"
if [ $? -ne 0 ]; then
    echo "Error during database creation"
    exit 1
fi

# Perform clustering
echo "Clustering sequences..."
mmseqs cluster "$MMSEQS_DB" "$MMSEQS_CLUSTER" "$MMSEQS_TMP" --min-seq-id 0.5 -c 0.8
if [ $? -ne 0 ]; then
    echo "Error during clustering"
    exit 1
fi

# Extract cluster information into TSV format
echo "Extracting cluster information into TSV file..."
mmseqs createtsv "$MMSEQS_DB" "$MMSEQS_DB" "$MMSEQS_CLUSTER" "$TSV_OUTPUT"
if [ $? -ne 0 ]; then
    echo "Error during TSV creation"
    exit 1
fi

echo "MMseqs2 clustering complete. TSV file created at: $TSV_OUTPUT"
```

 </details>   

### 3.3.5 Generate Presence-Absence Matrix
Now that you have the cluster information, you can compute the genome-gene family presence-absence matrix. This step can be done by parsing the output TSV file and counting the presence of gene families across the genomes.

You can write a Python script to generate this matrix from the mmseqs2_cluster.tsv file. The script will:

- Use the clustering output to create a binary matrix (1 for presence, 0 for absence). Refer the following python script to calculate the MMSeq matrix “04b_mmseq_matrix_construction.py”. Make sure to add the path of your working directory before running the file.

<details>
    <summary>04b_mmseq_matrix_construction.py</summary>

```
#Description: Construct the final MMSeq matrix from the tsv file generated from the software.
#Working directory: /projects/intro2gds/I2GDS2024/individual_folders/manisha/
#Input directory: results/04_mmseq2
#Output directory: results/04_mmseq2
#Author: Manisha Chaudhary
#Date: 10/24/2024

# Import necessary libraries
from collections import defaultdict, OrderedDict
import pandas as pd
import numpy as np
import csv
import os

# Set the path to the .faa files directory
path = "results/03_prodigal/cleaned_faa_files/"

# Initialize an empty list to hold the data from the .tsv file
data_list = []
with open("results/04_mmseq/mmseqs2_cluster_final.tsv", newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t')
    for row in spamreader:
        data_list.append(row)

# Create a defaultdict to store cluster-to-gene mappings
d = defaultdict(list)
for k, v in data_list:
    d[k].append(v)

# Create an ordered dictionary to maintain the order of clusters
d_order = OrderedDict(d)

# Initialize an empty DataFrame with 205 columns (adjust if needed)
ws = pd.DataFrame(index=np.arange(len(d_order)), columns=np.arange(205))
ws.fillna('', inplace=True)

# Create a list of all genome identifiers from .faa files
genome = []
for x in os.listdir(path):
    if x.endswith(".faa"):
        genome.append(x[0:13])  # Assuming the genome identifier is the first 13 characters

# Iterate through the clusters and genes, mapping them to the DataFrame
n = 0
for cluster, gene in d_order.items():
    ws.at[n, 0] = cluster[:]  # Add the cluster name
    ws.at[n, 1] = len(gene)   # Add the number of genes in the cluster
    unique_list = []
    for i in gene:
        meg = i[0:13]
        index = genome.index(meg)
        ws.at[n, index + 4] = str(ws.at[n, index + 4]) + str(i[:]) + ' '  # Map genes to the genome
        if meg not in unique_list:
            unique_list.append(meg)
    ws.at[n, 2] = len(unique_list)  # Add the number of unique genomes
    ws.at[n, 3] = False if ws.at[n, 2] == ws.at[n, 1] else True  # Identify paralogs
    n += 1

# Add a header row to the DataFrame
header_list = ['', '#genes', '#genomes', 'if paralogs']
for assembly in genome:
    header_list.append(assembly)
header_list.append('')
ws.columns = header_list

# Debugging statements
print("Current number of columns in ws:", ws.shape[1])  # Print the current number of columns
print("Length of header_list:", len(header_list))  # Print the length of the header list
print("Header list:", header_list)  # Print the header list to inspect its contents

# Set the columns of the DataFrame
ws.columns = header_list

# Save the final matrix to a CSV file
ws.to_csv("results/04_mmseq/mmseqsFinalMatrix.csv", index=False, header=True)
```

</details>

<p></p>
**GitHub Link**

For further reference, see the [MMSeq2 GitHub page](https://github.com/soedinglab/MMseqs2).

# Other information
For your reference, the steps outlined above correspond to specific sections of the overall workflow.
![common_pipeline](images/CommonPipeline.jpg)
