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
DOWNLOAD_DIR="/projects/intro2gds/I2GDS2024/individual_folders/manisha/data/assembly_data/"

mkdir -p "$DOWNLOAD_DIR"

#Use the datasets command to download all genomes listed in the file
datasets download genome accession --inputfile "$ACCESSION_FILE"

#Move downloaded files to the specified directory
#assuming files are named based on the accession codes
	
    for accession in $(cat "$ACCESSION_FILE") ; do
	echo "Moving files for accession: $accession"

#construct the expected filename
filename="${accession}.zip"
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
