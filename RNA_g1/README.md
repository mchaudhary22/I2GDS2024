# RNA-Seq-Pipeline (Jaret & Lili)

<img src=https://github.com/user-attachments/assets/939601c0-5314-47e8-8fbc-b7464a9e4c8e width=30% height=30%>


## Introduction

This page is a work in progress!

This repo explains a basic pipeline for RNA-Seq analysis in a linux HPC environment. It was developed as part of curriculum for Virginia Tech's Intro to Genomic Data Science course. This pipeline runs in Linux and relies on [FASTQC](https://github.com/s-andrews/FastQC), [Trimmomatic](https://github.com/timflutre/trimmomatic), [STAR](https://github.com/alexdobin/STAR), and [Featurecounts](https://subread.sourceforge.net/featureCounts.html). This example pipeline uses paired-end FastQ reads, but it could be altered for use with paired end data (see below). Throughout the code snippets in this repo we have attempted to make things as clear as possible. However, some paths may need to be changed depending on your method of installation or location of your files so be sure to be mindful of locations of your software and to read the code snippets carefully. 

Contact: Jaret Arnold (amichael19@vt.edu) or Lili Zebluim (liliz@vt.edu)

To do (before finalized):
- [x] Upload new pipeline image 
- [ ] Test All Code Blocks
- [ ] Fix i dont understand image
- [x] Add info on downloading the file
- [ ] Read through/edit blurbs and code snippets
- [x] Add References 

## Downloading Test files
To download the files used in this test workflow, run the following commands in your linux environment. 
```bash
#This is the forward read
wget 'https://drive.usercontent.google.com/download?id=10lMnWkwufamCqlRiHaN_8gL3u9UmK90f&export=download&authuser=1&confirm=t' -O D0C1_1.fq.gz
#This is the reverse read
wget 'https://drive.usercontent.google.com/download?id=1mWzTnFHuSoB43Ma0cfOSMmrK0-DARubV&export=download&authuser=1&confirm=t' -O D0C1_2.fq.gz
#This is the (arabidopsis) gtf file, used in aligning
wget 'https://drive.usercontent.google.com/download?id=160nLEzOYqO-fQ8_Pgbe5k0QuGeRZ3s-Z&export=download&authuser=1&confirm=t' -O arabidopsisgenome.gtf
#This is the (arabidopsis) genomic fasta file, used in indexing
wget 'https://drive.usercontent.google.com/download?id=1ZoM6vfRoSWphNLPwxcnT1fsXLFdbFprd&export=download&authuser=1&confirm=t' -O arabidopsisgenome.fa 
```

## FastQC
FastQC will be used to assess the quality of the raw reads and generate an html report detailing sequence quality, adapter contamination, GC content, etc. If FastQC is available on your computing environment, installation may be as easy as invoking module load. If not, try installing via the download. 

#### Installation via module load:
```bash
#load the fastqc module
module load FastQC
#test if install worked
fastqc --version

#TESTED
```

<details>
<summary> 
  
#### Installation via download:
</summary>
  
```bash
#download fastqc files
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
#unzip the file
unzip fastqc_v0.12.1
#add fastqc to your path (if this doesn't work try running using the full path)
export PATH:$PATH/to/fastqc
#refresh your paths in bashrc
source ~/.bashrc
#test if fastqc is working
fastqc --version

#UNTESTED
```
</details>

#### Running FastQC:

```bash
#move to location of reads (if neccesary)
#cd /path/to/reads   
#run fastqc on any files ending in .fq.gz and outputs in current file (-o)
fastqc *.fq.gz -o .

#TESTED
```

## Trimmomatic
Trimmomatic is used to remove adapter sequence contamination and low quality reads. Use your fastqc report to inform you how to best trim your reads. In the case of the demo file, trimming bases at the front (head) will improve the quality of our reads so we will use the HEADCROP option. If trimmomatic is available on your computing environment, installation may be as easy as invoking module load. If not, try installing via the download. Keep in mind larger files (as those in the example) can take some time to be trimmed, so consider batching your job using slurm/sbatch.
<need to doublecheck the installation code>

#### Installation via module load:
```bash
#load trimmomatic module
module load Trimmomatic
#test trimmomatic 
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar

#TESTED
```

<details>
<summary>
  
#### Installation via download: 
</summary>

```bash
#download trimmomatic files
wget https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip
#unzip the file
unzip Trimmomatic-0.39.zip
#test trimmomatic
java -jar /path/to/Trimmomatic-0.39/trimmomatic-0.39.jar 

#TESTED
```

</details>

#### Trimming paired-end reads:
```bash

#run trimmomatic on paired end (PE) data and output a log file of trims (-trimlog)
#firstly providing the forward read and the reverse read (D0C1_1.fq.gz D0C1_2.fq.gz),
#and then the names of the paired (D0C1_1.trim.fq.gz) and unpaired (D0C1_1un.trim.fq.gz) trimmed reads
#and the names for the paired/unpaired trimmed reads for the reverse read
#clipping the TruSeq3-PE adaptor (ILLUMINACLIP), the head of the read (HEADCROP)
#and dropping any reads which aren't minimum length 30 (MINLEN)

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
-trimlog trimlog.txt \
D0C1_1.fq.gz D0C1_2.fq.gz \
D0C1_1.trim.fq.gz D0C1_1un.trim.fq.gz \
D0C1_2.trim.fq.gz D0C1_2un.trim.fq.gz \
ILLUMINACLIP:/apps/packages/tinkercliffs-rome/trimmomatic/0.39/TruSeq3-PE.fa:2:30:10 MINLEN:30 HEADCROP:10 

#TESTED
```

> [!NOTE]
> Adapter selection may vary depending on the method of sequencing and therefore may need to be changed depending upon your data. Simply change TruSeq3-SE to the applicable adapter file provided by trimmomatic. 

Explanation of common trimming parameters: 

<img src="https://github.com/user-attachments/assets/013ce5ba-483c-4e31-aede-cf344c9a3eb9" width=75% height=75%>

<details>
  
<summary> 

  #### Trimming single-end reads: 
</summary>

```bash
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE \
-trimlog trimlog.txt \
example.fastq \
example.trim.fastq \
ILLUMINACLIP:TruSeq3-SE:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:20 \
MINLEN:36

#TESTED
```

</details>

After trimming, it is advisable to generate a second FastQC report to assess the success of trimming. For example:

```bash
fastqc *.trim.fq.gz -o .

```

## STAR
STAR is the program used to index and align reads to a reference genome. It is advisable to schedule STAR and other alignment processes on ARC as they can often require large memory and time requirements, particularly with more reads (see [example slurm scripts](#slurm-job-examples)). For instance, in our trials the test readmapping took ~50 min using 64GB of memory and 6 threads.
<need to check installation instructions>

#### Installation via download:
```bash
#download STAR
wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz
#untar the files
tar -xzf 2.7.11b.tar.gz
#move into the source files
cd STAR-2.7.11b/source
#compile STAR
make STAR

#TESTED
```

#### Genome indexing:
```bash
#create a dictory for the output genome index
mkdir /path/to/genomeindex

#Run STAR with 6 Threads (--runThreadN) using the test genomic fasta (--genomeFastaFiles) and gtf file (--sjdbGTFfile)
#to output in the folder previously created

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir genomeindex. \
--genomeFastaFiles arabidopsisgenome.fa \
--sjdbGTFfile arabidopsisgenome.gtf

#TESTED
```
<details>
<summary> Output file descriptions </summary>
<img src=https://github.com/user-attachments/assets/f2f91c6d-dc1a-470e-9087-fa28a96b18db width=100% height=100%>
</details>

> [!WARNING]
> Both STAR genome indexing and read mapping can be computationally intensive and require time. If working on ARC these should be submitted using slurm to efficiently schedule them. For Readmapping we have provided some example settings, but these can vary depending on your job and data size. To create a slurm batch job, simply nano <desiredfilename>.sh and input the your sbatch settings at the top of the file below the shebang. You may have to convert from dos to unix linebreaks if on windows using dos2unix <filename>. See the [example slurm scripts](#slurm-job-examples). 

#### Genome Read Mapping:
```bash
#!/bin/bash
# STAR Read Mapping - Test
#SBATCH --job-name=STAR-readmapping-test
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --ntasks=1             
#SBATCH -A <allocation> 
#SBATCH --time=48:00:00
#SBATCH -p normal_q
#SBATCH --output=STARslurmlogtest.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<user>

#Don't forget to put in your allocation and username!!

#run STAR using the full path (e.g. home/amichael19/software/STAR-2.7.11b/source/STAR) with
#set thread number (--runThreadN),
#the location of the previously created genome index (--genomeDir),
#read your RNA-Seq files in (--readFilesIn) first with the forward reads separated by a space,
#using compressed gzip files (readFilesCommand),
#and output as a BAM which is sorted (--outSAMtype)

path/to/STAR-2.7.11b/source/STAR --runThreadN 6 \
--genomeDir genomeindex/ \
--readFilesIn D0C1_1.fq.gz D0C1_2.fq.gz \ 
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate


#TESTED
```
<details>
<summary> Output file descriptions </summary>
<img src=https://github.com/user-attachments/assets/28f4b1b2-5b5b-4f8a-b03d-07b3ab967de3 width=100% height=100%>
</details>


## FeatureCounts
Feature counts is used to produce a matrix of genes and the count of their transcripts and it is contained within the subread package. Installation is easiest via conda (e.g. Miniconda3). 
<need to check installation instructions>

#### Installation via conda:
```bash
#load Miniconda3 to create an environment for Featurecounts
module load Miniconda3
#Create environment named subread (-n) on bioconda channel (-c) using the subread package 
conda create -n subread -c bioconda subread
#ARC current version to activate the previously created environment
source activate subread
#test if featurecounts work by invoking it
featureCounts


#TESTED
```

#### Running Feature counts:
```bash

#run featurecounts on paired end (p) data with an (-a) annotation file (.gtf/gff) using the previously created .bam (Aligned.sortedByCoord.out.bam) and output (-o) as a file named testcount.txt
featureCounts -p -a arabidopsisgenome.gtf -o testcount.txt Aligned.sortedByCoord.out.bam

#TESTED
```



## Slurm job Examples

<details>
<summary>FASTQC slurm job </summary>

```bash
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
```
</details>

<details>
<summary>Trimmomatic slurm job </summary>
  
```bash
#!/bin/bash
#SBATCH --job-name=trimmomatic_trim
#SBATCH --cpus-per-task=10              
#SBATCH -A <allocation>
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<user>

# Path to the Trimmomatic JAR file
trimmomatic_jar="/apps/packages/tinkercliffs-rome/trimmomatic/0.39/trimmomatic-0.39.jar"

# Path to the adapter file
adapter_file="/apps/packages/tinkercliffs-rome/trimmomatic/0.39/TruSeq3-PE.fa"

# Define the folder containing raw FASTQ files 
raw_folder="/home/amichael19/rawdata/PEreads/" 

# Define the folder containing trimmed FASTQ files 
trim_folder="/home/amichael19/results/trimmomatic/Petrim"

# Check if the folder exists
if [ ! -d "$raw_folder" ]; then
  echo "Error: The specified raw folder does not exist: $raw_folder"
  exit 1
fi

# Loop through all paired-end forward reads (_1.fq.gz files) in the raw folder
for forward_file in "$raw_folder"/*_1.fq.gz; do

  # Check if the forward file exists
  if [ ! -f "$forward_file" ]; then
    echo "Error: Forward file not found: $forward_file"
    continue
  fi

  # Get the base name for the current file (e.g., D0C1)
  base_name=$(basename "$forward_file" "_1.fq.gz")

  # Define the reverse file name
  reverse_file="${raw_folder}/${base_name}_2.fq.gz"

  # Check if the reverse file exists
  if [ ! -f "$reverse_file" ]; then
    echo "Error: Reverse file not found: $reverse_file"
    continue
  fi

  # Define output file names (output will be saved in the same folder as raw files)
  forward_paired_out="${trim_folder}/${base_name}_1.trim.fq.gz"
  forward_unpaired_out="${trim_folder}/${base_name}_1un.trim.fq.gz"
  reverse_paired_out="${trim_folder}/${base_name}_2.trim.fq.gz"
  reverse_unpaired_out="${trim_folder}/${base_name}_2un.trim.fq.gz"

  # Run Trimmomatic with java -jar
  java -jar "$trimmomatic_jar" PE \
    "$forward_file" "$reverse_file" \
    "$forward_paired_out" "$forward_unpaired_out" \
    "$reverse_paired_out" "$reverse_unpaired_out" \
    ILLUMINACLIP:"$adapter_file":4:30:10 MINLEN:30 HEADCROP:10

  # Check if the Trimmomatic command was successful
  if [ $? -eq 0 ]; then
    echo "Successfully processed $base_name."
  else
    echo "Error processing $base_name. Check the input files and parameters."
  fi

done

echo "All files have been processed and saved in $trim_folder."
```
</details>

<details>
<summary>STAR Genome Indexing slurm job </summary>

```
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
```
</details>


<details>
<summary>STAR Read Alignment slurm job </summary>
  
```bash
#!/bin/bash
# STAR Read Mapping - D3
#SBATCH --job-name=STAR-hybridreadmapping-Day3
#SBATCH --cpus-per-task=10
#SBATCH --mem=96G
#SBATCH --ntasks=1             
#SBATCH -A <allocation>
#SBATCH --time=48:00:00
#SBATCH -p normal_q
#SBATCH --output=STARslurmlogD3.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<user>

echo "Starting..."
date
time

cd /home/amichael19/software/STAR10424/STAR-2.7.11b/source

# Define variables

GENOME_DIR=/projects/intro2gds/I2GDS2024/individual_folders/jaret/data/STAR/PEgenomeindex
OUTPUT_DIR=/projects/intro2gds/I2GDS2024/individual_folders/jaret/data/bams/trimmobam/D3
TRIM_DIR=/projects/intro2gds/I2GDS2024/individual_folders/jaret/data/trimmedreads/trimmomatic/PEtrim
PREFIX=D3trimmo

# Run STAR to align paired-end RNA-seq data
./STAR --runThreadN 10 \
--genomeDir $GENOME_DIR \
--readFilesIn ${TRIM_DIR}/D3C1_1.fq.gz,${TRIM_DIR}/D3C2_1.fq.gz,${TRIM_DIR}/D3C4_1.fq.gz,${TRIM_DIR}/D3P1_1.fq.gz,${TRIM_DIR}/D3P2_1.fq.gz,${TRIM_DIR}/D3P3_1.fq.gz,${TRIM_DIR}/D3P4_1.fq.gz ${TRIM_DIR}/D3C1_2.fq.gz,${TRIM_DIR}/D3C2_2.fq.gz,${TRIM_DIR}/D3C4_2.fq.gz,${TRIM_DIR}/D3P1_2.fq.gz,${TRIM_DIR}/D3P2_2.fq.gz,${TRIM_DIR}/D3P3_2.fq.gz,${TRIM_DIR}/D3P4_2.fq.gz \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ${OUTPUT_DIR}/${PREFIX}_ 

echo "Finished!"
date
time

exit;
```
</details>

<details>
<summary>Feature Counts slurm job</summary>
  
```
#!/bin/bash
# Featurecounts
#SBATCH --job-name=Featurecounts
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH -A <allocation>
#SBATCH --time=48:00:00
#SBATCH -p normal_q
#SBATCH --output=featurecountslog.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<user>

# Declare Variables
ANNO_DIR=/home/amichael19/rawdata/gtfs/hybridgenome.gtf
OUTPUT_DIR=/projects/intro2gds/I2GDS2024/individual_folders/jaret/data/featurecounts/raw/hybridfeaturecounts.txt

# Run featurecounts
featureCounts -a $ANNO_DIR -o $OUTPUT_DIR \
/projects/intro2gds/I2GDS2024/individual_folders/jaret/data/bams/rawbam/D0/D0raw_Aligned.sortedByCoord.out.bam \
/projects/intro2gds/I2GDS2024/individual_folders/jaret/data/bams/rawbam/D1/D1raw_Aligned.sortedByCoord.out.bam \
/projects/intro2gds/I2GDS2024/individual_folders/jaret/data/bams/rawbam/D2/D2raw_Aligned.sortedByCoord.out.bam \
/projects/intro2gds/I2GDS2024/individual_folders/jaret/data/bams/rawbam/D3/D3raw_Aligned.sortedByCoord.out.bam \
/projects/intro2gds/I2GDS2024/individual_folders/jaret/data/bams/rawbam/D5/D5raw_Aligned.sortedByCoord.out.bam \
/projects/intro2gds/I2GDS2024/individual_folders/jaret/data/bams/rawbam/D7/D7raw_Aligned.sortedByCoord.out.bam \
/projects/intro2gds/I2GDS2024/individual_folders/jaret/data/bams/rawbam/D10/D10raw_Aligned.sortedByCoord.out.bam
```
</details>

## References
Data Download and workflow
https://www.youtube.com/watch?v=lG11JjovJHE

Installing FastQC
https://www.youtube.com/watch?v=5nth7o_-f0Q 

Installing Trimmomatic
https://www.youtube.com/watch?v=sG6b1aGEdCQ

Installing STAR
https://www.youtube.com/watch?v=RhEinNds1uc&t=0s

Making a Genomic Index in STAR
https://www.youtube.com/watch?v=Ju6PtQD-H34&t=320s

Installing Feature Counts 
https://www.youtube.com/watch?v=2VhNCYe8nQw

https://subread.sourceforge.net/featureCounts.html

https://github.com/alexdobin/STAR

https://github.com/timflutre/trimmomatic

http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

https://github.com/s-andrews/FastQC

https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

<add citations and any other refs here>
