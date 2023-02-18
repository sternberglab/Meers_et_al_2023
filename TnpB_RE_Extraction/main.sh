#!/bin/sh

##############################################################
#
# Script Name   :   main.sh
# Description   :   TnpB RE Boundary Detection For wRNA Building
# Author        :   Hoang C. Le
# Institution   :   Columbia University, Sternberg Laboratory
# Last updated  :   2023-02-16
# Requirements  :   All python packages in requirements.txt
#
##############################################################
read -p 'Output File Name Prefix:' file_prefix
read -p 'TnpB Bed File: ' tnpB_bed_file
read -p 'TnpB Protein Contigs Dataset Fasta:' tnpB_contig_fasta
read -p 'TnpB Intial Right End Expansion:' tnpB_initial_RE_expansion
read -p 'TnpB For Clustering Right End Expansion:' tnpB_final_RE_expansion
read -p 'How Many Clusters To Analyze:' cluster_limit


python3 $PWD/RE_extraction_TnpB.py tnpB_bed_file 0 tnpB_initial_RE_expansion file_prefix tnpB_contig_fasta tnpB_final_RE_expansion cluster_limit
