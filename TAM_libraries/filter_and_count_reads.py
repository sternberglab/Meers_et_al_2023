#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------
# Author: Matt W.G. Walker
# Email: matt.walker@columbia.edu
# Institution: Columbia University, Sternberg Laboratory
# Last updated: 09-22-2022
#----------------------------------
"""Trim paired-end input and output reads"""
#----------------------------------
from Bio import SeqIO
import pandas as pd
import numpy as np
import csv
from decimal import *
import os

def main():
  info = pd.read_excel("0_info/sample_info.xlsx", index_col=0).to_dict('index')

  datasets = [
    "220921_TAM_input",
    "220921_TAM_output_GstTnpB3",
    "220921_TAM_output_GstIscB"
  ]

  for dataset in datasets: 
    os.makedirs(info[dataset]['output_dir'], exist_ok=True)
    print("analyzing:", dataset)

    def analyze():
      create_info_csv()
      count_reads()
      trim_pre_N6()
      count_N6s()
      stats()
      abundance()
    
    def create_info_csv():
      input_dir = info[dataset]['fastq']
      output_dir_filtering = info[dataset]['output_dir'] + "info.csv"
      df = pd.DataFrame([input_dir], columns = ["sample_name"])
      df.to_csv(output_dir_filtering, index=False)
      print("Created info csv at ", output_dir_filtering)

    def count_reads():
      input_dir = "1_seq_data/" + info[dataset]['fastq'] + ".fastq"
      count = 0
      for rec in SeqIO.parse(input_dir, "fastq"):
        count += 1
      add_to_info_csv("total_reads", count)
      print("Total reads = ", count)

    def trim_pre_N6():
      input_dir = "1_seq_data/" + info[dataset]['fastq'] + ".fastq"
      output_dir_trimmed = info[dataset]['output_dir'] + "trimmed.fastq"
      PRE_N6_SEQ = "CTGGTTTAAGCTTACCGGTGCGGAGCCGTTGCAACAAGGGCGCTTGTGTGTAGGGAAG"
      trimmed = (
          rec[len(PRE_N6_SEQ):] for rec in SeqIO.parse(input_dir, "fastq") if rec.seq.startswith(PRE_N6_SEQ)
      )
      count_trimmed = SeqIO.write(trimmed, output_dir_trimmed, "fastq")
      print("Reads with pre_N6_seq = ", count_trimmed)
      add_to_info_csv("reads_with_pre_N6_seq", count_trimmed)

    def count_N6s():
      input_dir = info[dataset]['output_dir'] + "trimmed.fastq"
      output_dir = info[dataset]['output_dir'] + "counts_temp.csv"

      # create dictionary to use as counter for degenerate sequences
      N6_counter_dict = dict()
      N6_file = open("0_info/N6s.csv")
      for N6s in N6_file:
        N6s = N6s.strip('\n')
        (key) = N6s
        N6_counter_dict[key] = 0

      # Count matches per degenerate seq
      N6_parser = SeqIO.parse(input_dir, "fastq") 
      for rec in N6_parser:
        N6_seq = rec.seq[:6]
        if N6_seq in N6_counter_dict:
          N6_counter_dict[N6_seq] += 1
      sum_N6_counts = sum(N6_counter_dict.values())
      add_to_info_csv("reads_with_N6", sum_N6_counts)
      print("Reads with N6 = ", sum_N6_counts)

      # Save N6 counter as csv
      s = pd.Series(N6_counter_dict, name="N6_count")
      s.index.name = "N6"
      s.to_csv(output_dir)

    def stats():
      input_dir = info[dataset]['output_dir'] + "counts_temp.csv"
      counts = pd.read_csv(input_dir)

      # MEAN, STD, MIN, MAX
      mean = np.mean(counts.loc[:,"N6_count"],0)
      std = np.std(counts.loc[:,"N6_count"],0)
      min = np.percentile(counts.loc[:,"N6_count"],0)
      max = np.percentile(counts.loc[:,"N6_count"],100)

      # 90:10
      ninety = np.percentile(counts.loc[:,"N6_count"],90)
      ten = np.percentile(counts.loc[:,"N6_count"],10)
      sum_above_ninety = sum(e for e in counts.loc[:,"N6_count"] if e >= ninety)
      sum_below_ten = sum(e for e in counts.loc[:,"N6_count"] if e <= ten)
      ninety_ten = sum_above_ninety/sum_below_ten

      add_to_info_csv("mean", mean)
      add_to_info_csv("min", min)
      add_to_info_csv("max", max)
      add_to_info_csv("std", std)
      add_to_info_csv("90th %", ninety)
      add_to_info_csv("10th %", ten)
      add_to_info_csv("ninety-ten ratio", ninety_ten)

    def abundance():
      input_dir_counts = info[dataset]['output_dir'] + "counts_temp.csv"
      input_dir_info = info[dataset]['output_dir'] + "info.csv"
      output_dir_counts = info[dataset]['output_dir'] + "counts.csv"
      filtering_info = pd.read_csv(input_dir_info, index_col = "sample_name")
      reads_with_N6s = filtering_info["reads_with_N6"][0]

      with open(input_dir_counts) as N6_counts, open(output_dir_counts, 'w') as out:
          reader = csv.reader(N6_counts)
          writer = csv.writer(out)
          next(reader)
          writer.writerow(["N6", "N6_count", "N6_abundance"])
          for row in reader:
              writer.writerow([row[0], row[1], (Decimal(row[1])/reads_with_N6s*100)])

    def remove_large_files():
      trimmed_dir = info[dataset]['output_dir'] + "trimmed.fastq"
      counts_temp_dir = info[dataset]['output_dir'] + "counts_temp.csv"
      os.remove(trimmed_dir)
      os.remove(counts_temp_dir)

    def add_to_info_csv(colname, value):
      output_dir_info = info[dataset]['output_dir'] + "info.csv"
      df = pd.read_csv(output_dir_info)
      df.insert(len(df.columns), colname, value, allow_duplicates=True)
      df.to_csv(output_dir_info, index=False)


    analyze()

main()