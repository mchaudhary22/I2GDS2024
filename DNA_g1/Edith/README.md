<img src= width=30% height=30%>
**Figure Source:** [Nature Protocols, 2019](https://www.nature.com/articles/s41596-019-0223-x)

# MOWChIPseq for Histone Modification Profiling

## Introduction
This repository explains a general pipeline for data analysis with MOWChIP-Seq samples in a linux HPC environment. 
It was developed as part of the curriculum for Virginia Tech's ALS 5224 Intro to Genomic Data Science course. 

Epigenetic alterations, such as histone modifications, can help us understand how gene expression is regulated, how environmental factors can impact health, develop and evaluate  epigenetic therapeutic agents.

In order to identify significant histone modications in samples and observe patterns across patients with and without treatment, we obtained and prepared samples with MOWChIP-seq. However, samples may be used from previously published data can be used for running and testing the data analysis process. 

Contact: Edith Chen (edithchen@vt.edu)

# 1 Setup/ Preparation

## Create conda environment
A conda environment is a directory that contains a specific set of software packages, including libraries, dependencies, and Python versions.
This enables easy access to the specific versions of software packages needed for a project. For ChIP-seq analysis, I utilised the file ```spec_

To create an environment:

```bash
conda create --name <my-env>
```
Replace <my-env> with the name of your environment.


For downloading 

```bash
conda create —name <env> python=3.7.12 --file <this file>

```
Replace <this file> with the name of your file.



