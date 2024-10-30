#!/bin/bash
# Modified by Ying-Xian and Saehah for single-end SPAdes assembly
# 10-20-24

#SBATCH -J yx_spades_single_end # job name
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --time=2-00:00:00 # time: days-hours:minutes:seconds
#SBATCH --mem=128G # based on memory used when testing with interactive jobs
#SBATCH --mail-user=yingxian@vt.edu #enter desired email address for updates
#SBATCH --mail-type=BEGIN #include to get emailed when job begins
#SBATCH --mail-type=END #include to get emailed when job ends
#SBATCH --mail-type=FAIL #include to get emailed if job fails

# Define paths
data_dir="/projects/intro2gds/I2GDS2024/micro_g2/results/fastp_out_yingxian/"
output_dir="/projects/intro2gds/I2GDS2024/micro_g2/results/spades_out_yingxian/"
spades_path="/projects/leaph/.pyenv/versions/mambaforge/bin/spades.py"  # Updated path to spades.py
min_length=500  # Minimum length for contigs
min_coverage=2.0  # Minimum coverage for contigs

# Create output directories if they don't exist
mkdir -p $output_dir/contigs
mkdir -p $output_dir/scaffolds

# Change to the data directory
cd $data_dir

# Run SPAdes for each SE file
for f in *_QC.fastq.gz
do
    # Check if _contigs.fasta already exists
    if [ -f ${f%_QC.fastq.gz}/contigs.fasta ]; then
        echo "_contigs.fasta exists for ${f}, skipping SPAdes and proceeding to contig filtering..."
    else
        # Assemble using SPAdes with --isolate for single-end data
        echo "Assembling ${f} with SPAdes"
        python $spades_path --isolate -s $f -o ${f%_QC.fastq.gz} -t 3 -m 60

        # Check the log file for issues
        echo "Check the log file for ${f%_QC.fastq.gz} for any issues."
    fi

    # Proceed to filtering contigs based on length using seqtk
    cd ${f%_QC.fastq.gz}

    if [ -f contigs.fasta ]; then
        # Delete the old filtered file if it exists
        if [ -f $output_dir/contigs/${f%_QC.fastq.gz}_contigs.fasta ]; then
            rm $output_dir/contigs/${f%_QC.fastq.gz}_contigs.fasta
        fi

        # Use seqtk to filter out contigs shorter than min_length
        seqtk seq -L $min_length contigs.fasta > ${f%_QC.fastq.gz}_filtered_contigs.fasta

        # Filter based on coverage using awk
        awk -v min_cov=$min_coverage '
            /^>/ {
                split($0, a, "_");
                cov_pass = (a[6] >= min_cov);
            }
            {if (cov_pass) print}
        ' ${f%_QC.fastq.gz}_filtered_contigs.fasta > ${f%_QC.fastq.gz}_contigs.fasta

        # Move the final filtered contigs to the output directory
        mv ${f%_QC.fastq.gz}_contigs.fasta $output_dir/contigs/
    else
        echo "contigs.fasta not found for ${f}, skipping filtering step..."
    fi

    # Filter scaffolds based on length and coverage
    if [ -f scaffolds.fasta ]; then
        echo "scaffolds.fasta exists, proceeding to scaffold filtering..."

        # Delete the old filtered file if it exists
        if [ -f $output_dir/scaffolds/${f%_QC.fastq.gz}_scaffolds.fasta ]; then
            rm $output_dir/scaffolds/${f%_QC.fastq.gz}_scaffolds.fasta
        fi

        # Use seqtk to filter out scaffolds shorter than min_length
        seqtk seq -L $min_length scaffolds.fasta > ${f%_QC.fastq.gz}_filtered_scaffolds.fasta

        # Now filter based on coverage using awk
        awk -v min_cov=$min_coverage '
            /^>/ {
                split($0, a, "_");
                cov_pass = (a[6] >= min_cov);
            }
            {if (cov_pass) print}
        ' ${f%_QC.fastq.gz}_filtered_scaffolds.fasta > ${f%_QC.fastq.gz}_scaffolds.fasta

        # Rename the final output file as sample name + 'scaffolds.fasta'
        mv ${f%_QC.fastq.gz}_scaffolds.fasta $output_dir/scaffolds/${f%_QC.fastq.gz}_scaffolds.fasta
    else
        echo "scaffolds.fasta not found, skipping scaffold filtering step for ${f}..."
    fi

    cd ..

done

