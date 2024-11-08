<img src=https://github.com/user-attachments/assets/939601c0-5314-47e8-8fbc-b7464a9e4c8e width=30% height=30%>

# MOWChIPseq for Histone Modification Profiling

## Introduction
This repository explains a general pipeline for data analysis with MOWChIP-Seq samples in a linux HPC environment. 
It was developed as part of curriculum for Virginia Tech's ALS 5224 Intro to Genomic Data Science course. 

In order to identify significant histone modications in samples and observe patterns across patients with and without treatment, we obtained and prepared samples with MOWChIP-seq. However, test samples may be used from previously published data for Dr. Chang Lu's lab. 

Contact: Edith Chen (edithchen@vt.edu)


## Installing/Running 

#!/bin/bash
#SBATCH --account=chipseq
#SBATCH --partition=normal_q
#SBATCH --time=5:00:00
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=10


# go to the working directory
path=/projects/lu_lab/jdneice/TBIproject
cd $path

cd $path/rawdata

# trim the fastq data
trim_galore *.fastq
cd $path/results


############################### EXTENSION SETUP #################################
# Set the extend name as variables

FILES=$path/rawdata/*.fq
FQ=.fq
SAM=.sam
LOG=.log
BAM=.bam
SORT=_sort.bam
UNI=_unique.bam
PREBED=_pre.bed
BED=.bed
EXT=_extend.bed
EXPWIN=_extend_promotor_win.bed
EXPCOR=_extend_promotor_nor.bed
G100=_geno_win_100.bed
SORTG100=_geno_win_100_sort.bed
BedG=.bedGraph
NEBw=.bw


################################# FOLDER SETUP ##################################
# Create folders

mkdir Aligned_BAM
mkdir Aligned_SAM
mkdir SORT_BAM
mkdir UNIQUE_BAM
mkdir PRE_BED
mkdir BED
mkdir EXTEND_BED
mkdir EXT_PROMOTOR_WIN
mkdir EXT_GENO_WIN
mkdir Correlation
mkdir MACS
mkdir GenoWin100Sort
mkdir BedGraph
mkdir Nor_Ext_BW


# wc command allows you to count the number of lines, wordx, characters, and bytes of each given file or standard input and print the result
# -l--print the number of lines

input_length=$(wc -l < /projects/lu_lab/jdneice/TBIproject/input/Mouse/input.bed )


# Get the basename of the $FILES
for fn in $FILES
do
echo `basename "$fn"`
f=`basename "${fn%.*}"`
echo $f

done
## Align the reads with bowtie2
bowtie2 -p 16 -x /home/jdneice/Data/bw2_Index/mm10/mm10 -U $path/rawdata/$f$FQ -S $PWD/Aligned_SAM/$f$SAM 2>$PWD/Aligned_SAM/$f$LOG


## Converts SAM files to BAM files -b--output in the BAM format -q INT--skip alignments with MAPQ smaller than INT[0], here it is 10. 
# MAPQ: MAPping Quality. It equals −10 log10 Pr{mapping position is wrong}, 10 means Pr is 0.1
samtools view -bq 10 $PWD/Aligned_SAM/$f$SAM > $PWD/Aligned_BAM/$f$BAM

## Sort alignments -o--write the output to a file
samtools sort $PWD/Aligned_BAM/$f$BAM -o $PWD/SORT_BAM/$f$SORT

## Indexes SAM/BAM/CRAM files
samtools index $PWD/SORT_BAM/$f$SORT

## Remove potential PCR duplicates: if multiple read pairs have identical external coordinates, only retain the pair with highest mapping quality. -s remove duplicates for sing-end reads
samtools rmdup -s $PWD/SORT_BAM/$f$SORT $PWD/UNIQUE_BAM/$f$UNI


## Converts BAM file to BED file
bedtools bamtobed -i $PWD/UNIQUE_BAM/$f$UNI > $PWD/PRE_BED/$f$PREBED

## Subtract searches for features in B that overlap A by at least the number of base pairs given by the -f option
#-a--remove entire feature if any overlap
bedtools subtract -a $PWD/PRE_BED/$f$PREBED -b /home/jdneice/Data/blacklist/mm10/mm10.blacklist.bed > $PWD/BED/$f$BED

# Converts BED file to BAM file
bedToBam -i $PWD/BED/$f$BED -g /home/jdneice/Data/Genome/mm10/mm10genome.bed > $PWD/BED/$f$BAM

## Increase the size of each feature in a feature file by a user-defined number of bases. -b--Increase the BED/GFF/VCF entry by the same number base pairs in each direction. Integer.
bedtools slop -i $PWD/BED/$f$BED -g /home/jdneice/Data/Genome/mm10/mm10genome.bed -b 100 > $PWD/EXTEND_BED/$f$EXT


## Counts the coverage of data to reference genome
# Computes both the depth and breadth of coverage of features in file B on the features in file A.
# -counts--only report the countes of overlaps, don't compute fraction
# -b-- BAM/BED/GFF/VCF file (B)
# -a-- BAM/BED/GFF/VCF file (A)
bedtools coverage -counts -b $PWD/EXTEND_BED/$f$EXT -a /home/jdneice/Data/Promoter/mm10/mm10Promotor.bed > $PWD/EXT_PROMOTOR_WIN/$f$EXPWIN
bedtools coverage -counts -b $PWD/EXTEND_BED/$f$EXT -a /home/jdneice/Data/Genome/mm10/mm10genome_100.bed > $PWD/EXT_GENO_WIN/$f$G100


# Sort a file, -k1,1--sort only at positon 1, 
# -k2,2g--general-numeric-sort only at position 2, -k2,2n--numeric-sort only at position 2
# General numeric sort handles numbers in exponential notation, numeric sort is just a regular alphabetic sort that knows 10 comes after 9.
# -o--write the output to a file
# -u--output only the first of an equal run
sort -k1,1 -k2,2g -u -o $PWD/GenoWin100Sort/$f$SORTG100 $PWD/EXT_GENO_WIN/$f$G100


# Counts the length of BED file
ChIP_length=$(wc -l < $PWD/BED/$f$BED)

#Paste file1 file2--merge file1 and file2
# | pipe symbol
#OFS--output field separator, put the string "\t"--and works as 'tab'
#-v--asign the variable value
# $4/'$ChIP_length'*1000000--normalized value of chip reads, $4--the counts of each section of chip 
# $8/'$input_length'*1000000--normalized value of input reads, $8--the counts of each section of input 

paste $PWD/EXT_PROMOTOR_WIN/$f$EXPWIN /projects/lu_lab/jdneice/chip_practice/input/Mouse/input_extend_promotor_win.bed | awk -v OFS="\t" '{print $4/'$ChIP_length'*1000000-$8/'$input_length'*1000000}' > $PWD/Correlation/$f$EXPCOR
paste $PWD/EXT_GENO_WIN/$f$G100 /projects/lu_lab/jdneice/chip_practice/input/Mouse/input_extend_genome_win_100.bed | awk -v OFS="\t" '{print $1,$2,$3,$4/'$ChIP_length'*1000000-$8/'$input_length'*1000000}' > $PWD/BedGraph/$f$BedG


#################################### MACS2 ######################################
# Peak calling by MACS2

macs2 callpeak -t $PWD/BED/$f$BED -c /projects/lu_lab/jdneice/chip_practice/input/Mouse/input.bed -f BED -g hs -n $f -q 0.05 --outdir $PWD/MACS


# Sort the bed file (package: ucsc-bedsort)
bedSort $PWD/BedGraph/$f$BedG $PWD/BedGraph/$f$BedG

# Generate bedgraph file
bedGraphToBigWig $PWD/BedGraph/$f$BedG /home/jdneice/Data/Genome/mm10/mm10genome.bed $PWD/Nor_Ext_BW/$f$NEBw

done



#rm -r $PWD/Aligned_BAM
#/work/cascades/bz10/VCU_human/scripts .

