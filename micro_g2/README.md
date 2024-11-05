## Pipeline for assembling short read sequence data

Contact: Ying-Xian Goh (yingxian@vt.edu), Clayton Markham (cjmarkham@vt.edu), or Saehah Yi (shyi@vt.edu)

This pipeline integrates three essential tools — fastp, SPAdes, and CheckM2 — to ensure comprehensive processing and assessment of sequencing data. It begins by evaluating and trimming sequencing quality and adapters using fastp, followed by read assembly with SPAdes, and concludes with an analysis of assembly quality, specifically contamination and completeness, using CheckM2. All the code needed for this pipeline can be found at the `code_linux/` directory.

This repo was developed as part of the curriculum for Virginia Tech’s Introduction to Genomic Data Science course. Please make sure the script paths, data files, and working directory are correctly specified for smooth execution.

The image below, modified from [Del Angel et al., 2018[(https://f1000research.com/articles/7-148/v1) summarizes our 3-step pipeline:
![image](https://github.com/user-attachments/assets/8ccd585f-a9e8-4cca-9413-2e154410c557)

For another whole genome assembling pipeline example, see [micro_g1](https://github.com/LiLabAtVT/I2GDS2024/tree/main/micro_g1)!

## Dependencies & version information
Ensure the following dependencies are installed to run the pipeline:
- fastp v.0.23.4
- SPAdes v.4.0.0
- CheckM2 v.1.0.2

## Data download
Download genomic data from [NCBI Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) using the [SRA Toolkit v.3.1.1](https://github.com/ncbi/sra-tools).
- The SRA Toolkit is a collection of tools and libraries that allow users to access and manipulate data from SRA (a repository of molecular data from high-throughput sequencing platforms).

Useful tutorial videos for installing SRA Toolkit and downloading sequencing data from NCBI SRA
- [Download Fastq or SRA files - Whole Genome Sequencing Analysis. Step 1](https://www.youtube.com/watch?v=dZGf8D2WO44)
- [Downloading sequencing data on ubuntu/linux - SRA toolkit](https://www.youtube.com/watch?v=E1n-Z2HDAD0)

Use the following command to download SRA Toolkit:
```
curl --output sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz
```
See [here](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) for more detailed instructions on how to install SRA Toolkit.\
Use `fastq-dump` <SRR-of-interest> to retrieve the file that you want. For example: `fastq-dump SRR21285231`

## Step 1: Adapter trimming and quality filtering using [fastp](https://doi.org/10.1093/bioinformatics/bty560)
- Install **fastp** by following the instructions at the [GitHub](https://github.com/OpenGene/fastp). The method that was used in this tutorial was:
  ```
  wget http://opengene.org/fastp/fastp     #to download a file from the specified URL
  chmod a+x ./fastp     #to make the downloaded file executable
  ```
- `fastp.py`: Python script to handle raw FASTQ data and run **fastp** for adapter identification and quality control.
    - This script will look for the files in the specified directory, automatically identify adapters, trim adapters, and run quality control.
    - Specify/Modify the path to locate **fastp**` and change the input and output directories in `fastp.py` before use.
  - **Input**: FASTQ files (either single-end or paired-end)
    - for single-end data, specify read1 input by `-i` or `--in1`, and specify read1 output by `-o` or `--out1`.
    - for paired-end data, specify read2 input by `-I` or `--in2`, and specify read2 output by `-O` or `--out2`.
    - `fastp.py` in `code_linux/` is written for single-end data, you will need to modify the code your data is paired-end.
  - **Expected output**: `_QC.fastq.gz` gzip-compressed file
  - Run the following to execute the code.
    ```
    python3 fastp.py
    ```
> [!NOTE]
> The files required for Step 1 should be downloaded from NCBI SRA, or it can be directly copied from our directory using the following command:
> ```
> cp -r /projects/intro2gds/I2GDS2024/micro_g2/data/1_fastp_input /path/to/destination     #change '/path/to/destination' to your location
> ```

## Step 2: De novo assemble the genome using [SPAdes](https://doi.org/10.1002/cpbi.102)
- Install **SPAdes** by following the instructions at the [GitHub](https://github.com/ablab/spades). An alternative method that was used in this tutorial was installing **SPAdes** into a Anaconda virtual environment using the following code:
  ```
  module load Anaconda3     #to load the Anaconda3 module
  conda create --name SPAdes_v4.0.0_env     #to create a new Conda environment
  source activate SPAdes_v4.0.0_env     #to activate the newly created environment
  conda install bioconda::spades=4.0.0     #to install SPAdes v4 from the Bioconda channel into the active environment
  ```
- `spades.sh`: Bash script to performs a single-end (SE) genome assembly using SPAdes, then filters contigs and scaffolds based on length and coverage.
  - You need to specify the path to `spades` installed and adjust the input and output directories in the script. You can find the path for your SPAdes installed using:
    ```
    which spades.py
    ```
  - The `--isolate` option in SPAdes was used because the genome being sequenced in this tutorial is from a single, isolated organism (e.g., a bacterial isolate), rather than a mixed sample or metagenome.
  - **Input**: Single-end FASTQ files in your defined `data_dir`, where each file is named in the format `*_QC.fastq.gz` (output of fastp in previous step).
  - **Expected output**: It will create 2 directories in your defined `output_dir/`:
    - `contigs/`: Directory contains filtered contigs files for each sample.
    - `scaffolds/`: Directory contains filtered scaffolds files for each sample.
  - After modifications, run the following to execute the code.
    ```
    sbatch spades.sh
    ```
> [!NOTE]
> The files required for Step 2 (`*_QC.fastq.gz`) can be can be directly copied from our directory using the following command:
> ```
> cp -r /projects/intro2gds/I2GDS2024/micro_g2/data/2_spades_input /path/to/destination    #change '/path/to/destination' to your location
> ```

> [!NOTE]
> We understand that this is a lengthy step, so we've attached an example SLURM output file, `slurm-2758888.out`, in the `code_linux/` directory in case you need a reference for how the log file should look when the job runs successfully.   

## Step 3: Assess the quality of genome assemblies using [CheckM2](https://doi.org/10.1038/s41592-023-01940-w)
- Install **CheckM2** by following the instructions at the [GitHub](https://github.com/chklovski/CheckM2).
- For simplicity, you can just download CheckM2 from GitHub and run it directly without installing.
  - Retrieve the files:
    ```
    git clone --recursive https://github.com/chklovski/checkm2.git && cd checkm2
    ```
  - Create an appropriate Conda environment with prerequisites using the checkm2.yml file:
    ```
    conda env create -n checkm2 -f checkm2.yml
    conda activate checkm2
    ```
  - Finally, check the help menu for CheckM2:
    ```
    bin/checkm2 -h
    ```
- `checkm2.sh`: Bash script to assess assembly quality using **CheckM2**.
  - Please be reminded to change the path for the input and output directories before running the code. 
  - **Input**: Assembled genomes in your defined `data_dir`, where each file is named in the format `*_contigs.fasta` (output of SPAdes in previous step).
  - **Expected output**: It will create 2 folders and 2 files in your defined `output_dir/`:
    - `diamond_output/`: Directory contains DIAMOND results.
    - `protein_files/`: Directory contains .faa files.
    - `checkm2.log`: Log file summarizing the CheckM2 run.
    - `quality_report.tsv`: .tsv file summarizing the CheckM2 output. This is the **main output file** for running this program.
  - After modifications, run the following to execute the code.
    ```
    sbatch checkm2.sh
    ```
> [!NOTE]
> The files required for Step 3 can be can be directly copied from our directory using the following command:
> ```
> cp -r /projects/intro2gds/I2GDS2024/micro_g2/data/3_checkm2_input /path/to/destination     #change '/path/to/destination' to your location
> ```

> [!NOTE]
> Example `quality_report.tsv` output file were also provided in the `code_linux/` directory for reference.
