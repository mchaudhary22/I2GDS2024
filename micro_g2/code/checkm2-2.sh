#!/bin/bash
# By Ying-Xian with ChatGPT
# 10-31-24
# Use checkm2 to assess the completeness of the draft genome 

#SBATCH -J yx_checkm2
#SBATCH --account=introtogds # use 'leaph' for lab's allocation
#SBATCH --partition=normal_q
#SBATCH --time=1-00:00:00
#SBATCH --mem=64G
#SBATCH --mail-user=yingxian@vt.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=30

# Activate the environment
source activate checkm2

# Define paths
INPUT_DIR="/projects/leaph/shared/project_data/listeria_US/Lm_clinical_fasta/"
OUTPUT_DIR="/projects/leaph/yingxian/checkm2_for_all/output_lm_clinical_ncbi/"

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
