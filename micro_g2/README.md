## Pipeline for assembling short read sequence data

Contact: Ying-Xian Goh (yingxian@vt.edu), Clayton Markham (cjmarkham@vt.edu), or Saehah Yi (shyi@vt.edu)

This pipeline integrates three essential tools — fastp, SPAdes, and CheckM2 — to ensure comprehensive processing and assessment of short read sequencing data. It begins by evaluating and trimming sequencing quality and adapters using fastp, followed by read assembly with SPAdes, and concludes with an analysis of assembly quality, specifically contamination and completeness, using CheckM2. All the code needed for this pipeline can be found at the `code_linux/` directory.

> [!NOTE]
> This repo was developed as part of the curriculum for Virginia Tech’s Introduction to Genomic Data Science course. Please make sure the script paths, data files, and working directory are correctly specified for smooth execution.

The image below, modified from [Del Angel et al. (2018)](https://f1000research.com/articles/7-148/v1) summarizes our 3-step pipeline:

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

For this tutorial, we used three of the four draft genomes from [this publication](https://journals.asm.org/doi/10.1128/mra.00867-22), which can be downloaded from the SRA using the accession numbers `SRP394588`, `SRP394590`, and `SRP394593` or by using the already downloaded data presented in Step 1 below.

## Step 1: Adapter trimming and quality filtering using [fastp](https://doi.org/10.1093/bioinformatics/bty560)
- Install **fastp** by following the instructions at the [GitHub](https://github.com/OpenGene/fastp). The method that was used in this tutorial was:
  ```
  wget http://opengene.org/fastp/fastp     #to download a file from the specified URL
  chmod a+x ./fastp     #to make the downloaded file executable
  ```
- `fastp.py`: Python script to handle raw FASTQ data and run **fastp** for adapter identification and quality control.
    - This script will look for the files in the specified directory, automatically identify adapters, trim adapters, and run quality control.
    - Specify/modify the path to locate **fastp** and change the input and output directories in `fastp.py` before use.
      
> [!NOTE]
> `fastp.py` here is a Bash script that we wrote to automate fastp. It is different from the fastp software that you downloaded in the previous step using`wget http://opengene.org/fastp/fastp`
    
  - **Input**: FASTQ files (either single-end or paired-end)
    - For single-end data, specify read1 input by `-i` or `--in1`, and specify read1 output by `-o` or `--out1`.
    - For paired-end data, also specify read2 input by `-I` or `--in2`, and specify read2 output by `-O` or `--out2`.
    - Since `fastp.py` in `code_linux/` is written for single-end data, you will need to modify the code your data is paired-end.
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

## Step 2: _De novo_ assembly of genomes using [SPAdes](https://doi.org/10.1002/cpbi.102)
- Install **SPAdes** by following the instructions at the [GitHub](https://github.com/ablab/spades). An alternative method that was used in this tutorial was installing **SPAdes** into a Anaconda virtual environment using the following code:
  ```
  module load Anaconda3     #to load the Anaconda3 module
  conda create --name SPAdes_v4.0.0_env     #to create a new Conda environment
  source activate SPAdes_v4.0.0_env     #to activate the newly created environment
  conda install bioconda::spades=4.0.0     #to install SPAdes v4 from the Bioconda channel into the active environment
  ```
- `spades.sh`: Bash script to perform a single-end (SE) genome assembly using SPAdes, then filter contigs and scaffolds based on length and coverage.
  - If you installed SPAdes into the ARC instead of using a virtual environment, you need to specify the path to `spades` installed and adjust the input and output directories in the script. You can find the path for your SPAdes installed using:
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
> The files required for Step 2 (`*_QC.fastq.gz`) can be directly copied from our directory using the following command:
> ```
> cp -r /projects/intro2gds/I2GDS2024/micro_g2/data/2_spades_input /path/to/destination    #change '/path/to/destination' to your location
> ```

> [!NOTE]
> We understand that this is a lengthy step, so we've attached an example SLURM output file, `slurm-2758888.out`, in the `code_linux/` directory in case you need a reference for how the log file should look when the job runs successfully.   

## Step 3: Assess the quality of genome assemblies using [CheckM2](https://doi.org/10.1038/s41592-023-01940-w)
- Install **CheckM2** by following the instructions at the [GitHub](https://github.com/chklovski/CheckM2).
- The method that was used for this tutorial was installing **CheckM2** in a virtual environment using different code than the GitHub page which seems to run quicker on the VT ARC. This code is as follows:
  ```
  module load Anaconda3     #to load the Anaconda3 module
  conda create --name checkm2     #to create a new Conda environment
  source activate checkm2     #to activate the newly created environment
  conda install bioconda::checkm2     #to install checkm2 v1.0.2 from the Bioconda channel into the active environment
  ```
- You can open the help menu for CheckM2 to make sure that it correctly installed:
  ```
  checkm2 -h
  ```
- You can also do a test run to make sure CheckM2 was installed and is working correctly. See the [CheckM2 Github](https://github.com/chklovski/CheckM2) for the expected results:
  ```
  checkm2 testrun
  ```
- `checkm2.sh`: Bash script to assess assembly quality using **CheckM2**.
  - Please be reminded to change the path for the input, output, and diamond directories before running the code. 
  - **Input**: These are your assembled genomes from SPAdes. You can use either the `*_contigs.fasta` or `*_scaffolds.fasta` files, but you should assess the quality of whichever file you will continue to use in your analysis pipeline. Contigs are contiguous sequences of DNA that have been pieced together during the assembly process. Scaffolds are one order higher of genome assembly in which contigs are positioned relative to each other with gaps in between based on the known information about that genome. For genomes with high completeness, it should not matter much file is used.  In this tutorial, we have selected the contigs for further analysis.
  - The CheckM2 DIAMOND database will need to be installed for rapid annotation of the assemblies. This is done automatically through execution of `checkm2.sh`.
  - **Expected output**: It will create 2 folders and 2 files in your defined `output_dir/`:
    - `diamond_output/`: Directory contains DIAMOND results.
    - `protein_files/`: Directory contains .faa files.
    - `checkm2.log`: Log file summarizing the CheckM2 run.
    - `quality_report.tsv`: .tsv file summarizing the CheckM2 output. This is the **main output file** for running this program.
  - After the appropriate modifications, run the following to execute the code.
    ```
    sbatch checkm2.sh
    ```
> [!NOTE]
> The files required for Step 3 can be can be directly copied from our directory using the following command:
> ```
> cp -r /projects/intro2gds/I2GDS2024/micro_g2/data/3_checkm2_input /path/to/destination     #change '/path/to/destination' to your location
> ```

> [!NOTE]
> The example `quality_report.tsv` output file was also provided in the `code_linux/` directory for your reference.

-----

## Creating scatter plot and heatmap using RStudio
Contact: Ying-Xian Goh (yingxian@vt.edu), Saehah Yi (shyi@vt.edu), or Clayton Markham (cjmarkham@vt.edu)

This section explains the packages used and how the scatterplot and heatmap were generated for our data. All the data and code (written in R notebook) needed for these visualizations can be found in the `r_studio/data` and `r_studio/code` directories, respectively.

## Packages information
Ensure the following packages are installed and loaded:
- For creating scatter plot:
  -  `phyloseq`: For handling and analyzing microbiome data, specifically designed for ecological and genomic data formats. We need this for reading `phyloseq` object (.rds files in `r_studio/data`)
  - `ggplot2`: For creating complex and customizable visualizations using the Grammar of Graphics.
  - `vegan`: For ecological data analysis, offering tools for diversity analysis, ordination, and statistical analysis. We need this for beta diversity ordination analysis.
- For creating heatmap:
  - `pheatmap`: For creating heatmaps with additional customization options for annotating rows and columns, adjusting color schemes, clustering, and scaling. It simplifies the process of generating visually informative heatmaps from data matrices.
  - `reshape2`: For transforming data between wide and long formats, which is often necessary for preparing data for visualization or analysis. It includes functions like melt (to reshape data from wide to long format) and dcast (to go from long to wide format), making it easier to manipulate data frames.

> [!NOTE]
> You can install a package using code like the example below. Here’s how to install `phyloseq`:
> ```
> if (!require("BiocManager", quietly = TRUE))
>     install.packages("BiocManager")
>
> BiocManager::install("phyloseq")
> ```

## Scatter plot
- Input: `boo_wels_o_phyloseq.rds`, `boo_wels_o_taxa.rds`, `boo_wels_o_sample.names.rds`, and `meta_data_2.csv`.
<details>
    <summary>Step 1: Loading packages and setting working directory</summary>

```
# Load the packages
library(phyloseq)
library(ggplot2)
library(vegan)

# Setting working directory
setwd("/Users/Windows10/rstudio/i2gds/project_microplastic/") # change this to your working directory - the place where all your files are located

# Loading data/ Calling the object
boo_wels_o_phyobj <- readRDS("boo_wels_o_phyloseq.rds")
taxa <- readRDS("boo_wels_o_taxa.rds")
sample.names <- readRDS("boo_wels_o_sample.names.rds")

sample_data <- read.csv("meta_data_2.csv")
sample_data$Species <- ifelse(sample_data$Species == 'boo', 'L. booriae', 
                              ifelse(sample_data$Species == 'wels', 'L. welshimeri', 
                                     sample_data$Species))
```

</details>
<p></p>

<details>
    <summary>Step 2: Creating an ordination plot comparing the beta diversity of 2 groups</summary>

```
# Identify empty rows
empty_rows <- which(rowSums(otu_table(boo_wels_o_phyobj)) == 0)

# Filter out empty rows
boo_wels_o_phyobj_bray_pcoa <- prune_samples(sample_names(boo_wels_o_phyobj)[-empty_rows], boo_wels_o_phyobj)

# Perform ordination on the filtered data
ordination1 <- ordinate(boo_wels_o_phyobj_bray_pcoa, method = "PCoA", distance = "bray")

# Create the plot
plot1 <- plot_ordination(boo_wels_o_phyobj_bray_pcoa, ordination1, color="Species", title="Bray PCoA") +
  geom_point(size=2) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Add axis lines
    panel.border = element_rect(color = "black", fill = NA)
  )

# Add 95% CI around the data points
plot1 <- plot1 + stat_ellipse(level = 0.95)

# Calculate Bray-Curtis dissimilarity
braydist <- phyloseq::distance(boo_wels_o_phyobj, method = "bray")

# Perform PERMANOVA
permanova_result <- adonis2(braydist ~ Species, data = sample_data, permutations = 999)

# Extract p-value from PERMANOVA result
p_value <- permanova_result$`Pr(>F)`
p_value <- p_value[!is.na(p_value)]

# Add p-value annotation to the plot
plot1 <- plot1 +
  annotate(
    "text",
    x = 0.15,
    y = -0.5,
    label = paste("PERMANOVA P =", p_value),
    size = 4.5  # Increase the font size
  )

# Safe the plot
ggsave("bray_pcoa_plot.png", plot = plot1, device = "png", width = 8, height = 6)

# Show the plot
plot1
```

</details>
<p></p>

The output plot should look like the example below:
![image](https://github.com/user-attachments/assets/e4017ae7-79c1-4be2-9d5d-4e40a38321a0)

## Heatmap
- Input: `arg_pre_abs.csv`
<details>
    <summary>Loading packages and creating heatmapy</summary>

```
# Load the library
library(ggplot2)
library(reshape2)
library(pheatmap)

# Load the data
data <- read.csv("arg_pre_abs.csv", row.names = 1)

# Convert data to matrix
data_matrix <- as.matrix(data)

# Generate the heatmap
plot_2 <- pheatmap(
  data_matrix,
  color = colorRampPalette(c("lightgrey", "navy"))(50), # Color gradient
  cluster_rows = FALSE,                              # Add row dendrogram
  cluster_cols = FALSE,                              # Add column dendrogram
  main = "Heatmap of Antimicrobial Resistance Genes", # Title
  fontsize_row = 8,                                 # Row font size
  fontsize_col = 8,                                 # Column font size
  border_color = NA                                 # Remove lines inside the heatmap
)

# Safe the plot
ggsave("heatmap.png", plot = plot_2, device = "png", width = 8, height = 6)
```

</details>
<p></p>

The output plot should look like the example below:
![image](https://github.com/user-attachments/assets/b05e9af1-116b-4cb3-b7ed-5e7926820e09)
