#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------
# Author: Matt W.G. Walker
# Email: matt.walker@columbia.edu
# Institution: Columbia University, Sternberg and Landweber Laboratories
# Last updated: 2022-08-22
#----------------------------------
"""Determine the fold-change compared to input"""
#----------------------------------
import pandas as pd
import csv
import os
import math
from decimal import *

def main():
    sum_counts()
    abundance()
    log2_fold_change()
    summary()

def sum_counts(): 
    libraries = ["GstIscB", "GstTnpB3"]
    os.makedirs("3_sum_counts", exist_ok=True)
    for library in libraries:
        output_dir = "3_sum_counts/" + library + "_counts.csv"
        input = pd.read_csv("2_sample_counts/input/counts.csv", index_col=0, usecols = ["N6", "N6_count"], squeeze=True).to_dict()
        output = pd.read_csv("2_sample_counts/output_" + library + "/counts.csv", index_col=0, usecols = ["N6", "N6_count"], squeeze=True).to_dict()
        with open(output_dir, 'w') as out:
            writer = csv.writer(out)
            writer.writerow(["N6","input", "output"])
            for n6 in input:
                writer.writerow([n6, input[n6], output[n6]])

def abundance():
    libraries = ["GstIscB", "GstTnpB3"]
    os.makedirs("4_abundance", exist_ok=True)
    for library in libraries:
        input_dir = "3_sum_counts/" + library + "_counts.csv"
        output_dir = "4_abundance/" + library + "_abundance.csv"
        df = pd.read_csv(input_dir)
        input_sum = df["input"].sum()
        output_sum = df["output"].sum()
        with open(input_dir) as n6_counts, open(output_dir, 'w') as out:
            reader = csv.reader(n6_counts)
            writer = csv.writer(out)
            writer.writerow(["N6", "input", "output"])
            next(reader)
            for row in reader:
                writer.writerow([row[0], (Decimal(row[1])/input_sum*100),(Decimal(row[2])/output_sum*100)])

def log2_fold_change():
    libraries = ["GstIscB", "GstTnpB3"]
    os.makedirs("5_fold-change", exist_ok=True)
    for library in libraries:
        input_dir = "4_abundance/" + library + "_abundance.csv"
        output_dir = "5_fold-change/" + library + "_fold-change.csv"
        with open(input_dir) as input_abundance, open(output_dir, 'w') as out:
            reader = csv.reader(input_abundance)
            writer = csv.writer(out)
            next(reader)
            writer.writerow(["N6", "log2_fold-change"])
            for rows in reader:
                input = Decimal(rows[1])
                output = Decimal(rows[2])
                row = [rows[0]]
                if input == 0 and output == 0:
                    row.append("NA_uncovered_in_input_and_output")
                elif input == 0:
                    row.append("NA_uncovered_in_input_only")
                elif output == 0:
                    row.append("NA_dropout")
                else:
                    try:
                        row.append(math.log2(output/input))
                    except:
                        row.append("NA") # zeros are converted to NAs
                writer.writerow(row)
        print("output directory = ", output_dir)

def summary():
    libraries = ["GstIscB", "GstTnpB3"]
    os.makedirs("6_summary", exist_ok=True)
    for library in libraries:
        output_dir = "6_summary/" + library + ".csv"
        counts_input = pd.read_csv("3_sum_counts/" + library + "_counts.csv", index_col=0, usecols = ["N6", "input"], squeeze=True).to_dict()
        counts_output = pd.read_csv("3_sum_counts/" + library + "_counts.csv", index_col=0, usecols = ["N6", "output"], squeeze=True).to_dict()
        # abundance_input = pd.read_csv("4_abundance/" + library + "_abundance.csv", index_col=0, usecols = ["N6",  "input"], squeeze=True).to_dict()
        # abundance_output = pd.read_csv("4_abundance/" + library + "_abundance.csv", index_col=0, usecols = ["N6",  "output"], squeeze=True).to_dict()
        fc = pd.read_csv("5_fold-change/" + library + "_fold-change.csv", index_col=0, usecols = ["N6", "log2_fold-change"], squeeze=True).to_dict()
        with open(output_dir, 'w') as out:
            writer = csv.writer(out)
            writer.writerow(["N6","counts_input", "counts_output","log2(fold-change)"]) #"abundance_input", "abundance_output",
            for n6 in counts_input:
                writer.writerow([n6, counts_input[n6], counts_output[n6], fc[n6]]) #abundance_input[n6], abundance_output[n6], 

main()