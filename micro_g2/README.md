## Pipeline for assembling short read sequence data
This pipeline integrates three essential tools — fastp, SPAdes, and CheckM — to ensure comprehensive processing and assessment of sequencing data. It begins by evaluating and trimming sequencing quality and adapters using fastp, followed by read assembly with SPAdes, and concludes with an analysis of assembly quality, specifically contamination and completeness, using CheckM.

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

Use `fastq-dump` <SRR-of-interest> to retrieve the file that you want. For example: `fastq-dump SRR21285231`

## Step 1. To perform adapter trimming and quality filtering using [fastp](https://doi.org/10.1093/bioinformatics/bty560)
The code is located in the following path: `code/`.
- `fastp.py` : Python script to handle raw FASTQ data and run fastp for adapter identification and quality control.
  - This script will look for the files in the specified directory, automatically identify adapters, trim adapters, and run QC.
  - Need to specify the path to locate `fastp` and change the input and output directories in `fastp.py` before use.
  - Need to download `fastp` at the [GitHub](https://github.com/OpenGene/fastp)
    
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
