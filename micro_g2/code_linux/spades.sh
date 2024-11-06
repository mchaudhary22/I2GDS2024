#!/bin/bash
# Modified by Ying-Xian and Saehah for single-end SPAdes assembly
# 10-20-24

#SBATCH -J yx_spades_single_end # job name
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --time=2-00:00:00 # time: days-hours:minutes:seconds
#SBATCH --mem=128G # based on memory used when testing with interactive jobs
#SBATCH --mail-user=yingxian@vt.edu #enter desired email address for updates
#SBATCH --mail-type=BEGIN #include to get emailed when job begins
#SBATCH --mail-type=END #include to get emailed when job ends
#SBATCH --mail-type=FAIL #include to get emailed if job fails

# Define paths
data_dir="/projects/intro2gds/I2GDS2024/micro_g2/results/fastp_out_yingxian/" # Updated to your input_path
output_dir="/projects/intro2gds/I2GDS2024/micro_g2/results/spades_out_yingxian/" # Updated to your output_path

# Create output directories if they don't exist
mkdir -p $output_dir/contigs
mkdir -p $output_dir/scaffolds

# Activate the Anaconda environment
module load Anaconda3
source activate SPAdes_v4.0.0_env

# Change to the data directory
cd $data_dir

# Run SPAdes for each SE file
for f in *_QC.fastq.gz
do
    # Check if _contigs.fasta already exists
    if [ -f ${f%_QC.fastq.gz}/contigs.fasta ]; then
        echo "_contigs.fasta exists for ${f}, skipping SPAdes..."
    else
        # Assemble using SPAdes with --isolate for single-end data
        echo "Assembling ${f} with SPAdes"
        spades.py --isolate -s $f -o ${f%_QC.fastq.gz} -t 3 -m 60

        # Check the log file for issues
        echo "Check the log file for ${f%_QC.fastq.gz} for any issues."
    fi

    # Move the final contigs and scaffolds to the output directory
    if [ -f ${f%_QC.fastq.gz}/contigs.fasta ]; then
        mv ${f%_QC.fastq.gz}/contigs.fasta $output_dir/contigs/${f%_QC.fastq.gz}_contigs.fasta
    else
        echo "contigs.fasta not found for ${f}, skipping move step..."
    fi

    if [ -f ${f%_QC.fastq.gz}/scaffolds.fasta ]; then
        mv ${f%_QC.fastq.gz}/scaffolds.fasta $output_dir/scaffolds/${f%_QC.fastq.gz}_scaffolds.fasta
    else
        echo "scaffolds.fasta not found for ${f}, skipping move step..."
    fi

done

