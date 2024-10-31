#!/bin/bash
# RNAseq pipeline
#SBATCH --job-name=RNAseq-pipeline
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --ntasks=1             
#SBATCH -A <allocation>
#SBATCH --time=48:00:00
#SBATCH -p normal_q
#SBATCH --output=RNAseqpipetest.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<email>

#!!Set your allocation and email before submitting!!

echo "Pipeline Starting!"
date
time

#Download files we will use
#This is the forward read
wget 'https://drive.usercontent.google.com/download?id=10lMnWkwufamCqlRiHaN_8gL3u9UmK90f&export=download&authuser=1&confirm=t' -O D0C1_1.fq.gz
#This is the reverse read
wget 'https://drive.usercontent.google.com/download?id=1mWzTnFHuSoB43Ma0cfOSMmrK0-DARubV&export=download&authuser=1&confirm=t' -O D0C1_2.fq.gz
#This is the (arabidopsis) gtf file, used in aligning
wget 'https://drive.usercontent.google.com/download?id=160nLEzOYqO-fQ8_Pgbe5k0QuGeRZ3s-Z&export=download&authuser=1&confirm=t' -O arabidopsisgenome.gtf
#This is the (arabidopsis) genomic fasta file, used in indexing
wget 'https://drive.usercontent.google.com/download?id=1ZoM6vfRoSWphNLPwxcnT1fsXLFdbFprd&export=download&authuser=1&confirm=t' -O arabidopsisgenome.fa 

echo "files successfully downloaded!"

#load the req modules
module load FastQC
module load Trimmomatic
module load Miniconda3

echo "modules loaded!"

#create conda environments we will use 
conda create -n STAR -c bioconda star
conda create -n subread -c bioconda subread

echo "environments created!"

#fastqc raw reads
fastqc *.fq.gz -o .

echo "fastQC completed on raw reads!"

#trim via trimmomatic
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
-trimlog trimlog.txt \
D0C1_1.fq.gz D0C1_2.fq.gz \
D0C1_1.trim.fq.gz D0C1_1un.trim.fq.gz \
D0C1_2.trim.fq.gz D0C1_2un.trim.fq.gz \
ILLUMINACLIP:/apps/packages/tinkercliffs-rome/trimmomatic/0.39/TruSeq3-PE.fa:2:30:10 MINLEN:30 HEADCROP:10 

echo "trimming completed!"

#rerun fastqc on trims
fastqc *.trim.fq.gz -o .

echo "fastQC completed on trimmed reads!"

#invoke STAR env and run STAR indexing
source activate STAR

mkdir genomeindex/

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir genomeindex \
--genomeFastaFiles arabidopsisgenome.fa \
--sjdbGTFfile arabidopsisgenome.gtf

echo "STAR completed genome indexing!"

#run read mapping

STAR --runThreadN 6 \
--genomeDir genomeindex \
--readFilesIn D0C1_1.trim.fq.gz D0C1_2.trim.fq.gz \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \

echo "STAR completed read mapping!"

#deactivate STAR env and activate subread
conda deactivate 
source activate subread

featureCounts -p -a arabidopsisgenome.gtf -o testcount.txt Aligned.sortedByCoord.out.bam

echo "Feature counts created successfully!"
echo "Pipeline finished successfully!"
date
time