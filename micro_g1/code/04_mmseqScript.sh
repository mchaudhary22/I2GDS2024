# Description: Performs MMSeq2 analysis on accumulated .faa files obtained from prodial analysis
# Working directory: /projects/intro2gds/I2GDS2024/individual_folders/manisha/
# Input files: results/03_prodigal/cleaned_faa_files/
# Output directory: results/04_mmseq
# Author: Manisha Chaudhary
# Date: 09/28/2024

# Set directories
FAA_DIR="results/03_prodigal/cleaned_faa_files/"
OUTPUT_DIR="results/04_mmseq"
MMSEQS_DB="$OUTPUT_DIR/mmseqs2_db_cleaned"
MMSEQS_CLUSTER="$OUTPUT_DIR/mmseqs2_cluster_cleaned"
MMSEQS_TMP="$OUTPUT_DIR/tmp"
TSV_OUTPUT="$OUTPUT_DIR/mmseqs2_cluster_final.tsv"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check if the faa files directory exists
if [ ! -d "$FAA_DIR" ]; then
    echo "FAA directory not found: $FAA_DIR"
    exit 1
fi

# Create MMseqs2 database
echo "Creating MMseqs2 database..."
mmseqs createdb "$FAA_DIR"/*.faa "$MMSEQS_DB"
if [ $? -ne 0 ]; then
    echo "Error during database creation"
    exit 1
fi

# Perform clustering
echo "Clustering sequences..."
mmseqs cluster "$MMSEQS_DB" "$MMSEQS_CLUSTER" "$MMSEQS_TMP" --min-seq-id 0.5 -c 0.8
if [ $? -ne 0 ]; then
    echo "Error during clustering"
    exit 1
fi

# Extract cluster information into TSV format
echo "Extracting cluster information into TSV file..."
mmseqs createtsv "$MMSEQS_DB" "$MMSEQS_DB" "$MMSEQS_CLUSTER" "$TSV_OUTPUT"
if [ $? -ne 0 ]; then
    echo "Error during TSV creation"
    exit 1
fi

echo "MMseqs2 clustering complete. TSV file created at: $TSV_OUTPUT"
