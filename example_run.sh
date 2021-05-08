#!/bin/bash
PY=../pypy3.7-v7.3.4-linux64/bin/pypy

# Part 1
echo
echo "Part 1"
echo
$PY -u parsimony-informative.py -i example_data -m variation_graph_taxonid_2697049_outgroupid_none.ma -b none -r 0 -f "OEAV139851"

# Part 2
echo
echo "Part 2"
echo
$PY -u haplotype_network.py example_data

# Filter net_all.txt
awk -F'\t' '{print $1","$2}' example_data/net_all.txt > net.csv
sed -i '1i\Source,Target' net.csv

# Part 3
echo
echo "Part 3"
echo
$PY main_path_example.py
