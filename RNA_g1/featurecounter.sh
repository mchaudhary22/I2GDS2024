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