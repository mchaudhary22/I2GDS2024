#!/bin/bash
#SBATCH --account=chipseq
#SBATCH --partition=normal_q
#SBATCH --time=5:00:00

# This script prepares a SINGLE INPUT FILE for use in the base_precorr.sh code
# Before running:
# 1. Index genome of interest
# 2. Obtain blacklist of genes from genome of interest
# 3. fetchChromSizes for genome of interest (*genome.bed file)
# 4. bedtools makewindows to separate genome into windows of 100 bp (*genome_100.bed file)
# 5. Generate Promoters with promoter.sh (*Promoter.bed file) (for more info read promoter.sh script)


cd /projects/lu_lab/edith/practice/rawdata/input                  # changes to the directory with the input data



(bowtie2 -x /home/edithchen/data/bw2_index/mm10/mm10  -U *.fq -S input.sam) 2> alignSumm.log  


# Identifies location of indexed genome and aligns the input.fq file, the aligned reads go to input.sam 
# It outputs any standard errors from the subshell (in parentheses) to a file called alignSumm.log


samtools view -bq 10 ./input.sam > ./input.bam          # converts the input.sam file into a bam file called input.bam, skipping alignments with MAPQ smaller than 10.

samtools sort input.bam -o input_sort.bam               # sorts input.bam by leftmost coordinate and outputs result to the file input_sort.bam

samtools index input_sort.bam                           # indexes the file input_sort.bam

samtools rmdup -s input_sort.bam input_unique.bam       # removes potential pcr duplicates for single-end reads and output to the file input_unique.bam

bedtools bamtobed -i input_unique.bam > input_pre.bed   # convert input_unique.bam to BED format file input_pre.bed

bedtools subtract -a ./input_pre.bed -b /home/edithchen/data/blacklist/mm10/mm10.blacklist.bed > input.bed   

# subtracts features from file "B" (the genome blacklist) that intersect features from file "A" (input_pre.bed file) and output to input.bed


bedtools slop -i input.bed -g /home/edithchen/genomesetup/fetchromsize/mm10/mm10genome.bed -b 100 > input_extend.bed 

# increases the size of each feature in input.bed by 100bp, increasing the size of the entry by the samenumber of base pairs in each direction and outputs to input_extend.bed


bedtools coverage -counts -b input_extend.bed -a /home/edithchen/data/promoter/mm10/mm10Promotor.bed > input_extend_promotor_win.bed 

# compares the features in "A" (the promotors in the genome of interest) to the file "B" (input_extend.bed) in search of overlaps and outputs overlaps to input_extend_promotor_win.bed


bedtools coverage -counts -b input_extend.bed -a /home/edithchen/genomesetup/fetchromsize/mm10/mm10genome_100.bed > input_extend_geno_win_100.bed 

#compares the features in "A" (the genome separated into windows of 100 bps) to the file "B" (input_extend.bed) in search of overlaps and outputs overlaps to input_extend_geno_win_100.bed
