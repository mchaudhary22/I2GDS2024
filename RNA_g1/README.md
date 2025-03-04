# RNA-Seq-Pipeline 💻⚙️(Jaret & Lili)

![RNAseq workflow ](https://github.com/user-attachments/assets/4e5f4768-09be-4302-809c-eff8fbda234f)

<a name="top"></a>

## Introduction

This page is a work in progress!
This repo explains a basic pipeline for RNA-Seq analysis. This pipeline relies heavily on FASTQC, Trimmomatic, STAR, and Featurecounts. For installation instructions see below... (add more?)  
!!To download files see RNA_G2!! 

## FastQC
<Blurb about fastqc>
<need to doublecheck the installation code>

Installation via module load:
```bash
module load FastQC

#UNTESTED
```

Installation on the cluster:
```bash
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
unzip fastqc_v0.12.1
export PATH:$PATH/to/fastqc
source ~/.bashrc
fastqc --version

#UNTESTED
```

Running FastQC:

```bash
cd /path/to/reads
fastqc *.fastq

#UNTESTED
```

## Trimmomatic
<Blurb about trimmomatic>
<need to doublecheck the installation code>
  
Installation via module load:
```bash
module load Trimmomatic
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar

#UNTESTED
```

Installation on cluster:

```bash
wget https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar 

#UNTESTED
```

Trimming:
```bash
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE \
-trimlog trimlog.txt \
demo.fastq \
demo.trim.fastq \
ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \

#UNTESTED
```



## STAR
<Blurb about STAR>
<need to check installation instructions>

Installation via download:
```bash
wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz
tar -xzf 2.7.11b.tar.gz
cd STAR-2.7.11b
cd STAR/source
make STAR

#UNTESTED
```


## FeatureCounts
<Blurb about Featurecounts>
<need to check installation instructions>

Installation via conda:
```bash
module load Miniconda3
conda create -n subread -c bioconda subread
source activate subread
featureCounts

#UNTESTED
```

Running Feature counts:
```bash
featureCounts -a <annotation file> -o <path/to/outputfile.txt> <bamfile>

#UNTESTED
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

[Back to top](#top)
