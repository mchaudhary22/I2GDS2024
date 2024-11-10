<img src= width=30% height=30%>
**Figure Source:** [Nature Protocols, 2019](https://www.nature.com/articles/s41596-019-0223-x)

# MOWChIPseq for Histone Modification Profiling

## Introduction
This repository explains a general pipeline for data analysis with MOWChIP-Seq samples in a linux HPC environment. 
It was developed as part of curriculum for Virginia Tech's ALS 5224 Intro to Genomic Data Science course. 

In order to identify significant histone modications in samples and observe patterns across patients with and without treatment, we obtained and prepared samples with MOWChIP-seq. However, samples may be used from previously published data for Dr. Chang Lu's lab for running and testing the data analysis process. 

Contact: Edith Chen (edithchen@vt.edu)


## Create conda environment
Python is required to create a conda environment


Replace <my-env> with the name of your environment.

```bash
conda create —name <env> python=3.7.12 --file <this file>

```


