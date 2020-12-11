# VENUS：a Viral genome Evolution Network Analysis System

## Introduction

Comprehensive analyses of viral genomes can provide a global picture of SARS-CoV-2 transmission and help to predict the oncoming trends of the pandemic. However, the rapid accumulation of SARS-CoV-2 genomes presents an unprecedented data size and complexity that has exceeded the capacity of existing methods in constructing evolution network through virus genotyping. The VENAS seeks to apply reliable computational algorithms to build an integrative genomic analysis system that enables researchers to trace viral mutations along the transmission routes using the daily updated SARS-CoV-2 genomes.

VENAS can construct the network from an alignment file containing 10k sequences in about 10 minutes, including:
1. Use a neighbor-joining based evolution network instead of the traditional phylogenetic tree to combine viral genomic alterations with transmission events, which avoids the bifurcation limitations of the tree structure. 
1. Use parsimonious information sites (PISs) instead of all variants to construct the evolution network, reducing the interference caused by sequencing errors.
1. Use Hamming distances adjusted by the minor allele frequency to construct the viral genome evolution network, enhancing network construction's robustness.
1. Use a community detection method to transform the evolution network into a two-dimensional isomorphic topological space, and apply a network disassortativity trimming algorithm to extract the backbone network of the topological space, which can be used to trace the evolution of viruses and detect core mutations among distinct strains.
1. Taking advantage of the multi-core and multi-threaded features in high-performance computers and developed a highly parallel network construction pipeline, which could handle a massive amount of viral genome in a limited time.
1. Tracing the viral evolution network through different transmission stages, revealing several transmission branches and mutation "hot-spots" of SARS-CoV-2.

## Pre-requisites
VENAS requires python 3 with PyPy, argparse, pandas, numpy (<http://www.numpy.org/>), networkx, CDlib, matplotlib, biopython (<http://biopython.org/wiki/Main_Page>), click, tqdm libraries installed.
If you want to provide a fasta file as input file, VENAS also needs the MAFFT (<https://mafft.cbrc.jp/alignment/software/>) in the executable path. Then you can use the “multi_mafft.py” to perform a multi-threaded multiple sequence alignment.

## Installation
Cloning the repository via the following commands `$ git clone https://github.com/qianjiaqiang/VENAS.git`

## Basic Usage
This section presents some basic usages of VENAS. We assume here that all the scripts are in the system path.

### Part 1: Effective parsimony-informative site (ePIS) finding and Minor allele frequency calculating

1.	All the sequences in a fasta file should be multiple-aligned to obtain a consensus alignment in ma format.
2.	The site was a parsimony-informative which contained at least two types of nucleotides, and at least two of them occur with a minimum frequency of two.
3.	In the default parameters, the PIS was effective if the number of unambiguous bases ≥80% of the total genomes.
4.	In the default parameters, all sequences with any ePIS containing bases other than "A,T,C,G" will be removed and output to the “rm_non-ATCG_genomes.ma”.
*Note:* The i parameter is the directory where the input file is located. The f parameter is the reference genome sequence id in the ma file.

```#!bash

python -u parsimony-informative.py -i example_data -m variation_graph_taxonid_2697049_outgroupid_none.ma -b none -r 0 -f "OEAV139851"```

### Part 2: Viral genome evolution network construction

```#!bash

python -u haplotype_network.py example_data```

### Part3: Topological classification and major path recognition

*Note:* Only the first two columns are needed in the output “net_all.txt” file of the Part 2 step, which can be handled as described below.
```#!bash
awk -F'\t' '{print $1","$2}' example_data/net_all.txt > net.csv
sed -i '1i\Source,Target' net.csv```
Example input net_all.csv:
```
Source,Target
1,57
5,6
1,210
10,23
10,59
10,69
3,191
10,91```

If you have already processed the net.csv file, you are ready for Part 3.

```#!bash

python main_path.py```

The result net.csv and nodeTable.csv files are in the current working directory. You can visualize the result viral genome evolution network using a general relationship graph or force-directed graph tools, such as the web-based Apache Echarts (<https://echarts.apache.org/>), d3.js (<https://d3js.org/>), or the application-based Gephi (recommend).
The net.csv contains all the edges in the evolution network. The edgeTable.csv file contains the major transmission paths of the evolution network. The nodes with value 1000 in the nodeTable.csv file represent the main nodes on the major transmission paths.

