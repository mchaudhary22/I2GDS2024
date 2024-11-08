#!/bin/bash
# STAR Read Mapping - Test
#SBATCH --job-name=STAR-readmapping-test
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --ntasks=1             
#SBATCH -A <allocation> 
#SBATCH --time=48:00:00
#SBATCH -p normal_q
#SBATCH --output=STARreadmaplog.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<user>

#Don't forget to put in your allocation and username!!

#run STAR using the full path (e.g. home/amichael19/software/STAR-2.7.11b/source/STAR) with
#set thread number (--runThreadN),
#the location of the previously created genome index (--genomeDir),
#read your RNA-Seq files in (--readFilesIn) first with the forward reads separated by a space,
#using compressed gzip files (readFilesCommand),
#and output as a BAM which is sorted (--outSAMtype)

path/to/STAR-2.7.11b/source/STAR runMode alignReads \
--runThreadN 6 \
--genomeDir genomeindex \
--readFilesCommand zcat \
--readFilesIn D0C1_1.trim.fq.gz D0C1_2.trim.fq.gz \
--outSAMtype BAM SortedByCoordinate