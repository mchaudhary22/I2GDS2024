#!/bin/bash
# By Ying-Xian and Clayton with ChatGPT
# 10-31-24
# Use checkm2 to assess the completeness and contamination of the draft genome 

#SBATCH -J yx_checkm2 # modify the job name accordingly
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --time=1-00:00:00
#SBATCH --mem=64G
#SBATCH --mail-user=yingxian@vt.edu # change to your personal email
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=30 # modify if needed

# Activate the environment
source activate checkm2

# Define paths
INPUT_DIR="/projects/intro2gds/I2GDS2024/micro_g2/results/spades_out_yingxian/contigs/" # need to be modified accordingly
OUTPUT_DIR="/projects/intro2gds/I2GDS2024/micro_g2/results/checkm2_out_yingxian/" # need to be modified accordingly

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run CheckM2
echo "Starting CheckM2 analysis..."
checkm2 predict \
    --threads $SLURM_CPUS_PER_TASK \
    --input "$INPUT_DIR" \
    -x fasta \
    --output-directory "$OUTPUT_DIR"

# Check if CheckM2 completed successfully
if [ $? -eq 0 ]; then
    echo "CheckM2 prediction completed successfully. Results saved in $OUTPUT_DIR."
else
    echo "Error: CheckM2 prediction failed. Please check the logs for more information."
    exit 1
fi
