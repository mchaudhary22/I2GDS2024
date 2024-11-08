#!/bin/bash
#SBATCH --job-name=trimmomatic_trim
#SBATCH --cpus-per-task=10              
#SBATCH -A <allocation>
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<user>

# Path to the Trimmomatic JAR file
trimmomatic_jar="/apps/packages/tinkercliffs-rome/trimmomatic/0.39/trimmomatic-0.39.jar"

# Path to the adapter file
adapter_file="/apps/packages/tinkercliffs-rome/trimmomatic/0.39/TruSeq3-PE.fa"

# Define the folder containing raw FASTQ files 
raw_folder="/home/amichael19/rawdata/PEreads/" 

# Define the folder containing trimmed FASTQ files 
trim_folder="/home/amichael19/results/trimmomatic/Petrim"

# Check if the folder exists
if [ ! -d "$raw_folder" ]; then
  echo "Error: The specified raw folder does not exist: $raw_folder"
  exit 1
fi

# Loop through all paired-end forward reads (_1.fq.gz files) in the raw folder
for forward_file in "$raw_folder"/*_1.fq.gz; do

  # Check if the forward file exists
  if [ ! -f "$forward_file" ]; then
    echo "Error: Forward file not found: $forward_file"
    continue
  fi

  # Get the base name for the current file (e.g., D0C1)
  base_name=$(basename "$forward_file" "_1.fq.gz")

  # Define the reverse file name
  reverse_file="${raw_folder}/${base_name}_2.fq.gz"

  # Check if the reverse file exists
  if [ ! -f "$reverse_file" ]; then
    echo "Error: Reverse file not found: $reverse_file"
    continue
  fi

  # Define output file names (output will be saved in the same folder as raw files)
  forward_paired_out="${trim_folder}/${base_name}_1.trim.fq.gz"
  forward_unpaired_out="${trim_folder}/${base_name}_1un.trim.fq.gz"
  reverse_paired_out="${trim_folder}/${base_name}_2.trim.fq.gz"
  reverse_unpaired_out="${trim_folder}/${base_name}_2un.trim.fq.gz"

  # Run Trimmomatic with java -jar
  java -jar "$trimmomatic_jar" PE \
    "$forward_file" "$reverse_file" \
    "$forward_paired_out" "$forward_unpaired_out" \
    "$reverse_paired_out" "$reverse_unpaired_out" \
    ILLUMINACLIP:"$adapter_file":4:30:10 MINLEN:30 HEADCROP:10

  # Check if the Trimmomatic command was successful
  if [ $? -eq 0 ]; then
    echo "Successfully processed $base_name."
  else
    echo "Error processing $base_name. Check the input files and parameters."
  fi

done

echo "All files have been processed and saved in $trim_folder."