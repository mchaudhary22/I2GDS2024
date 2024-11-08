#!/bin/bash
# STAR Genome Indexing
#SBATCH --job-name=STAR-hybridgenomeindex 
#SBATCH --cpus-per-task=6
#SBATCH --ntasks=1             
#SBATCH -A <allocation>
#SBATCH --time=62:00:00
#SBATCH -p normal_q
#SBATCH --output=STARslurmlog.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<user>

echo "Starting..."
date
time

cd /home/amichael19/software/STAR10424/STAR-2.7.11b/source

# Define variables

GENOME_DIR=/home/amichael19/results/STAR/trimmomaticindex
FASTA_DIR=/home/amichael19/rawdata/genomicfastas/hybridgenome.fa
GTF_DIR=/home/amichael19/rawdata/gtfs/hybridgenome.gtf

# Run STAR to create genome index
./STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir $GENOME_DIR \
--genomeFastaFiles $FASTA_DIR \
--sjdbGTFfile $GTF_DIR \
--sjdbOverhang 99

echo "Finished!"
date
time

exit;