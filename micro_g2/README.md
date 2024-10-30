## Pipeline for assembling short read sequence data
This pipeline integrates three essential tools — fastp, SPAdes, and CheckM2 — to ensure comprehensive processing and assessment of sequencing data. It begins by evaluating and trimming sequencing quality and adapters using fastp, followed by read assembly with SPAdes, and concludes with an analysis of assembly quality, specifically contamination and completeness, using CheckM2.

## Dependencies & Version Information
Ensure the following dependencies are installed to run the pipeline:
- fastp v.0.23.4
- SPAdes v.4.0.0
- checkm2 v.1.0.2

## Data download
Download genomic data from [NCBI Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) using the [SRA Toolkit v.3.1.1](https://github.com/ncbi/sra-tools).
- The SRA Toolkit is a collection of tools and libraries that allow users to access and manipulate data from SRA (a repository of molecular data from high-throughput sequencing platforms).

Useful tutorial videos for installing SRA Toolkit and downloading sequencing data from NCBI SRA
- [Download Fastq or SRA files - Whole Genome Sequencing Analysis. Step 1](https://www.youtube.com/watch?v=dZGf8D2WO44)
- [Downloading sequencing data on ubuntu/linux - SRA toolki](https://www.youtube.com/watch?v=E1n-Z2HDAD0)

Use the following command to download SRA Toolkit:
```
curl --output sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz
```
See [here](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) for more detailed instructions on how to install SRA Toolkit.\
Use `fastq-dump` <SRR-of-interest> to retrieve the file that you want. For example: `fastq-dump SRR21285231`

## Step 1: Adapter trimming and quality filtering using [fastp](https://doi.org/10.1093/bioinformatics/bty560)
The code is located in the following path: `code/fastp.py`.
- `fastp.py` : Python script to handle raw FASTQ data and run fastp for adapter identification and quality control.
  - This script will look for the files in the specified directory, automatically identify adapters, trim adapters, and run quality control.
  - Specify the path to locate `fastp` and change the input and output directories in `fastp.py` before use.
  - Install `fastp` by following the instructions at the [GitHub](https://github.com/OpenGene/fastp). The method that was used in this tutorial was:
    
    ```
    wget http://opengene.org/fastp/fastp
    chmod a+x ./fastp
    ```
  - **Input**: FASTQ files (either single-end or paired-end)
    - for single-end data, specify read1 input by `-i` or `--in1`, and specify read1 output by `-o` or `--out1`.
    - for paired-end data, specify read2 input by `-I` or `--in2`, and specify read2 output by `-O` or `--out2`.
    - `fastp.py` in `code/` is written for single-end data, need to modify if your data is paired-end
  - **Expected output**: `_QC.fastq.gz` gzip-compressed file
  - Run the following to execute the code.
    
    ```
    python3 fastp.py
    ```
## Step 2: Assembly using [SPAdes](https://doi.org/10.1089/cmb.2012.0021)
The code is located in the following path: `code/spades.sh`.
- `spades.sh` : Linux bash script to assemble quality-filtered and trimmed reads into contigs.
  - **note to Ying-Xian to edit this line** This script will look for the files in the specified directory, automatically identify adapters, trim adapters, and run QC.
  - **note to Ying-Xian to edit this line** Specify the path to locate `fastp` and change the input and output directories in `fastp.py` before use.
  - Install `spades` by following the instructions at the [GitHub](https://github.com/ablab/spades). An alternative method that was used in this tutorial was installing SPAdes into a Anaconda virtual environment using the following code:
    ```
    module load Anaconda3
    conda create --name SPAdes_v4.0.0_env
    source activate SPAdes_v4.0.0_env
    conda install bioconda::spades
    ```
  - SPAdes v4.0.0 was used in this tutorial. To ensure the same version rather than the most recent version of SPAdes is installed, specify `conda install bioconda::spades=4.0.0`.



**Note to Ying-Xian**: You will want to make sure that this code is present at the top of your script for running SPAdes in the virtual environment:
```
    module load Anaconda3
    source activate SPAdes_v4.0.0_env
    ```
