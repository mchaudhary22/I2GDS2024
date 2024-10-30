#Description: Construct the final MMSeq matrix from the tsv file generated from the software.
#Working directory: /projects/intro2gds/I2GDS2024/individual_folders/manisha/
#Input directory: results/04_mmseq2
#Output directory: results/04_mmseq2
#Author: Manisha Chaudhary
#Date: 10/24/2024

# Import necessary libraries
from collections import defaultdict, OrderedDict
import pandas as pd
import numpy as np
import csv
import os

# Set the path to the .faa files directory
path = "results/03_prodigal/cleaned_faa_files/"

# Initialize an empty list to hold the data from the .tsv file
data_list = []
with open("results/04_mmseq/mmseqs2_cluster_final.tsv", newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t')
    for row in spamreader:
        data_list.append(row)

# Create a defaultdict to store cluster-to-gene mappings
d = defaultdict(list)
for k, v in data_list:
    d[k].append(v)

# Create an ordered dictionary to maintain the order of clusters
d_order = OrderedDict(d)

# Initialize an empty DataFrame with 205 columns (adjust if needed)
ws = pd.DataFrame(index=np.arange(len(d_order)), columns=np.arange(205))
ws.fillna('', inplace=True)

# Create a list of all genome identifiers from .faa files
genome = []
for x in os.listdir(path):
    if x.endswith(".faa"):
        genome.append(x[0:13])  # Assuming the genome identifier is the first 13 characters

# Iterate through the clusters and genes, mapping them to the DataFrame
n = 0
for cluster, gene in d_order.items():
    ws.at[n, 0] = cluster[:]  # Add the cluster name
    ws.at[n, 1] = len(gene)   # Add the number of genes in the cluster
    unique_list = []
    for i in gene:
        meg = i[0:13]
        index = genome.index(meg)
        ws.at[n, index + 4] = str(ws.at[n, index + 4]) + str(i[:]) + ' '  # Map genes to the genome
        if meg not in unique_list:
            unique_list.append(meg)
    ws.at[n, 2] = len(unique_list)  # Add the number of unique genomes
    ws.at[n, 3] = False if ws.at[n, 2] == ws.at[n, 1] else True  # Identify paralogs
    n += 1

# Add a header row to the DataFrame
header_list = ['', '#genes', '#genomes', 'if paralogs']
for assembly in genome:
    header_list.append(assembly)
header_list.append('')
ws.columns = header_list

# Debugging statements
print("Current number of columns in ws:", ws.shape[1])  # Print the current number of columns
print("Length of header_list:", len(header_list))  # Print the length of the header list
print("Header list:", header_list)  # Print the header list to inspect its contents

# Set the columns of the DataFrame
ws.columns = header_list

# Save the final matrix to a CSV file
ws.to_csv("results/04_mmseq/mmseqsFinalMatrix.csv", index=False, header=True)
