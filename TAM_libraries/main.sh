#!/bin/sh

##############################################################
#
# Script Name   :   main.sh
# Description   :   Integration site library analysis
# Author        :   Matt W.G. Walker
# Institution   :   Columbia University, Sternberg Laboratory
# Last updated  :   February 15, 2023
# Requirements  :   All python packages in requirements.txt
#
##############################################################

python3 $PWD/generate_deg_seqs.py

python3 $PWD/filter_and_count_reads.py

python3 $PWD/analyze_libraries.py