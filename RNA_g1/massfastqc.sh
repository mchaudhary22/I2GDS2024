#!/bin/bash
# Mass FastQC
#SBATCH --job-name=Batch_FastQC
#SBATCH --cpus-per-task=6
#SBATCH -A <allocation>
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<user>

#Move into directory w reads
cd /projects/intro2gds/I2GDS2024/individual_folders/jaret/data/trimmedreads/trimmomatic/PEtrim

#run fastqc
/home/amichael19/software/FastQC92424/FastQC/fastqc *.fq.gz -o /projects/intro2gds/I2GDS2024/individual_folders/jaret/data/QualityControl/fastqc/trimmomatic