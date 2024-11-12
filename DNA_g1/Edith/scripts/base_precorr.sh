#!/bin/bash

Group=GROUPGOESHERE
subgroup=SUBgroupGOESHERE
Histone=HISTONEGOESHERE

cd /projects/lu_lab/edith/MIA/ChIP_seq/$Group/$Histone/$subgroup

cd Raw_Data/
gunzip *.gz
trim_galore *.fastq

cd ..

FILES=$PWD/Raw_Data/*.fq
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
EXG4000=_extend_geno_4000.bed
EXPCOR=_extend_promotor_nor.bed
EXGCOR=_extend_geno_4000_nor.bed
G100=_geno_win_100.bed
SORTG100=_geno_win_100_sort.bed
BedG=.bedGraph
NEBw=.bw

mkdir Aligned_BAM
mkdir Aligned_SAM
mkdir SORT_BAM
mkdir UNIQUE_BAM
mkdir BED
mkdir EXTEND_BED
mkdir EXT_PROMOTOR_WIN
mkdir EXT_GENO_WIN
mkdir Correlation
mkdir MACS
mkdir GenoWin100Sort
mkdir BedGraph
mkdir Nor_Ext_BW
mkdir Enhancer_bed
mkdir MACS_enhancer

input_length=$(wc -l < /projects/lu_lab/Gaoshan/MIA/ChIP_seq/Input/$Group/input.bed )

for fn in $FILES
do
echo `basename "$fn"`
f=`basename "${fn%.*}"`
echo $f


# Read into bowtie2
bowtie2 -p 16 -x /home/gaoshanli/Data/mm10/Sequence/Bowtie2Index/genome -U $PWD/Raw_Data/$f$FQ -S $PWD/Aligned_SAM/$f$SAM 2>$PWD/Aligned_SAM/$f$LOG

# Read into samtools
samtools view -bq 10 $PWD/Aligned_SAM/$f$SAM > $PWD/Aligned_BAM/$f$BAM

samtools sort $PWD/Aligned_BAM/$f$BAM -o $PWD/SORT_BAM/$f$SORT

samtools index $PWD/SORT_BAM/$f$SORT

samtools rmdup -s $PWD/SORT_BAM/$f$SORT $PWD/UNIQUE_BAM/$f$UNI

bedtools bamtobed -i $PWD/UNIQUE_BAM/$f$UNI > $PWD/BED/$f$PREBED

bedtools subtract -a $PWD/BED/$f$PREBED -b /home/gaoshanli/Data/mm10/mm10.blacklist.bed > $PWD/BED/$f$BED

bedtools subtract -a $PWD/BED/$f$BED -b /home/gaoshanli/Data/Promotor/mm10/mm10Promotor.bed > $PWD/Enhancer_bed/$f$BED

bedToBam -i $PWD/Enhancer_bed/$f$BED -g /home/gaoshanli/Data/bed_winbed/mm10/mm10genome.bed > $PWD/Enhancer_bed/$f$BAM

bedtools slop -i $PWD/Enhancer_bed/$f$BED -g /home/gaoshanli/Data/bed_winbed/mm10/mm10genome.bed -b 100 > $PWD/EXTEND_BED/$f$EXT

bedtools coverage -counts -b $PWD/EXTEND_BED/$f$EXT -a /home/gaoshanli/Data/Promotor/mm10/mm10Promotor.bed > $PWD/EXT_PROMOTOR_WIN/$f$EXPWIN

bedtools coverage -counts -b $PWD/EXTEND_BED/$f$EXT -a /home/gaoshanli/Data/bed_winbed/mm10/mm10genome_win_4000.bed > $PWD/EXT_GENO_WIN/$f$EXG4000

bedtools coverage -counts -b $PWD/EXTEND_BED/$f$EXT -a /home/gaoshanli/Data/bed_winbed/mm10/mm10genome_win_100.bed > $PWD/EXT_GENO_WIN/$f$G100

sort -k1,1 -k2,2g -u -o $PWD/GenoWin100Sort/$f$SORTG100 $PWD/EXT_GENO_WIN/$f$G100

ChIP_length=$(wc -l < $PWD/Enhancer_bed/$f$BED)

paste $PWD/EXT_PROMOTOR_WIN/$f$EXPWIN /projects/lu_lab/Gaoshan/MIA/ChIP_seq/Input/$Group/input_extend_promotor_win.bed | awk -v OFS="\t" '{print $4/'$ChIP_length'*1000000-$8/'$input_length'*1000000}' > $PWD/Correlation/$f$EXPCOR

paste $PWD/EXT_GENO_WIN/$f$EXG4000 /projects/lu_lab/Gaoshan/MIA/ChIP_seq/Input/$Group/input_extend_genome_4000win.bed | awk -v OFS="\t" '{print $4/'$ChIP_length'*1000000-$8/'$input_length'*1000000}' > $PWD/Correlation/$f$EXGCOR

macs2 callpeak -t $PWD/Enhancer_bed/$f$BED -c /projects/lu_lab/Gaoshan/MIA/ChIP_seq/Input/$Group/input.bed -f BED -g mm -n $f -p 0.01 --outdir $PWD/MACS_enhancer

paste $PWD/EXT_GENO_WIN/$f$G100 /projects/lu_lab/Gaoshan/MIA/ChIP_seq/Input/$Group/input_extend_genome_100win.bed | awk -v OFS="\t" '{print $1,$2,$3,$4/'$ChIP_length'*1000000-$8/'$input_length'*1000000}' > $PWD/BedGraph/$f$BedG

bedSort $PWD/BedGraph/$f$BedG $PWD/BedGraph/$f$BedG

bedGraphToBigWig $PWD/BedGraph/$f$BedG /home/gaoshanli/Data/bed_winbed/mm10/mm10genome.bed $PWD/Nor_Ext_BW/$f$NEBw

done

/projects/lu_lab/Gaoshan/MIA/ChIP_seq/Summary_enhancer.sh .
