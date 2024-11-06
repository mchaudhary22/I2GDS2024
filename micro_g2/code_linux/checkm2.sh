#!/bin/bash
# By Ying-Xian and Clayton with ChatGPT
# 11-05-24
# Use CheckM2 to assess the completeness and contamination of genome assemblies

#SBATCH -J yx_checkm2 # modify the job name accordingly
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --time=1-00:00:00
#SBATCH --mem=64G
#SBATCH --mail-user=cjmarkham@vt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=30 #modify if needed

# Activate the environment
source activate checkm2

# Define paths
INPUT_DIR="/projects/intro2gds/I2GDS2024/micro_g2/results/spades_out_yingxian/contigs/"
OUTPUT_DIR="/projects/intro2gds/I2GDS2024/micro_g2/results/checkm2_out_cjmarkham/"
DIAMOND_DIR="/projects/intro2gds/I2GDS2024/micro_g2/software"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Install checkm2 DIAMOND database if it doesn't exist
if [ -f $DIAMOND_DIR/CheckM2_database/uniref100.KO.1.dmnd ]; then
    echo "DIAMOND database exists."
else
    echo "Downloading DIAMOND database."
    checkm2 database --download --path $DIAMOND_DIR
fi

# Run CheckM2
echo "Starting CheckM2 analysis..."
checkm2 predict \
    --threads $SLURM_CPUS_PER_TASK \
    --input $INPUT_DIR \
    -x .fasta \
    --force \
    --output-directory $OUTPUT_DIR \
    --database_path $DIAMOND_DIR/CheckM2_database/uniref100.KO.1.dmnd

# Check if CheckM2 completed successfully
if [ $? -eq 0 ]; then
    echo "CheckM2 prediction completed successfully. Results saved in $OUTPUT_DIR."
else
    echo "Error: CheckM2 prediction failed. Please check the logs for more information."
    exit 1
fi
