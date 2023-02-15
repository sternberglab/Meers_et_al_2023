#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------
# Author: Matt W.G. Walker
# Email: matt.walker@columbia.edu
# Institution: Columbia University, Sternberg Laboratory
# Last updated: 2023-02-15
#----------------------------------
"""Generate all possible degenerate sequences"""
#----------------------------------
import csv
import itertools

def main():
    generate_N6s()

def generate_N6s():
    Ns = 'ACTG'
    output_dir_deg = "0_info/N6s.csv"

    with open(output_dir_deg, 'w') as out:
        write = csv.writer(out)
        for output in itertools.product(Ns, repeat=6):
            deglist = [''.join(output)]
            write.writerow(deglist)

main()