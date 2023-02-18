#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#----------------------------------
# Author: Hoang C. Le
# Email: lehoang@wharton.upenn.edu
# Institution: Columbia University, Sternberg Laboratories
# Last updated: 2023-02-16
#----------------------------------
"""Pipeline To Find RE Boundaries For TnpB IS Elements: Part 1"""
#----------------------------------
import pandas as pd
import sys
import os
import re

#----------------------------------
"""1. Get Putative Coordinates"""
#----------------------------------
def get_putative_coordinates(tnpB_coords, expansion_in, expansion_out, file_prefix):
    putative_RE_coords_bed = file_prefix + "_putative_RE_coords_bed" + ".bed"
    final_bed = tnpB_coords

    #Expand
    for i in final_bed.index:
        start = final_bed[1][i]
        stop = final_bed[2][i]

        if final_bed[5][i] == "+":
            final_bed[1][i] = max(stop - expansion_in,0)
            final_bed[2][i] = stop + expansion_out
        else:
            final_bed[1][i] = max(start - expansion_out,0)
            final_bed[2][i] = start + expansion_in

    final_bed.to_csv(putative_RE_coords_bed, sep='\t', header=False, index=False)
    return putative_RE_coords_bed

#----------------------------------
"""2. Get Sequences and Cluster"""
#----------------------------------
def convert_clstr_to_txt(file, file_new):
    file1 = open(file, 'r')
    file2 = open(file_new, 'w')
    Lines = file1.readlines()

    # since I am using a while loop I need the following index
    LineIndex = 0
    while LineIndex < len(Lines):
        # check if >Cluster lines are there
        FF = re.search(r'^(>Cluster.*)', Lines[LineIndex])
        # if the lines re found, then prep it for first column
        if FF != None:
            #print(FF.group())
            firstcol = FF.group().replace('>','')
            file2.write(firstcol+',')
            print(firstcol, end =",")
            LineIndex = LineIndex + 1
            while LineIndex < len(Lines):
                KK = re.search(r'^(.*\*)', Lines[LineIndex])
                if KK != None:

                    seccol = KK.group()
                    seccol = re.sub(r'.*>','',seccol)
                    seccol = re.sub(r'\.\.\. \*','',seccol)
                    print(seccol)
                    file2.write(seccol+'\n')
                    break
                LineIndex = LineIndex + 1
        LineIndex = LineIndex + 1

    file1.close()
    file2.close()

def convert_txt_to_df(file, file_new):
    df = pd.read_csv(file, sep='\t')
    df["clstr"] = df["clstr"].astype(int)
    df["clstr_rep"] = df["clstr_rep"].astype(int)

    ##Create List Where Index Is Cluster Number and Item Is Representative
    reps = [""] * int(df["clstr"].max() + 1)

    for i in df.index:
        if df["clstr_rep"][i] == 1:
            clstr = df["clstr"][i]

            reps[clstr] = df["id"][i]

    df['Group_Rep'] = df['clstr'].apply(lambda x: reps[x])
    print(df)

    df.to_csv(file_new, sep='\t', index=False)

    filtered_list_df = pd.DataFrame({'Group_Reps': reps})
    print(filtered_list_df)
    filtered_list_df.to_csv(file_new + "_ID_list", sep='\t', header=False, index=False)

def get_sequences_and_cluster(tnpB_contig_fasta, file_prefix, putative_RE_coords_bed):
    putative_RE_sequence_fasta = file_prefix + "_putative_RE_sequence" + ".fa"
    os.popen('bedtools getfasta -s -name -fo ' + putative_RE_sequence_fasta + ' -fi ' tnpB_contig_fasta + ' -bed ' + putative_RE_coords_bed).read()

    putative_RE_sequence_fasta_clst_99 = file_prefix + "putative_RE_sequence_fasta_clst_99" + ".fa"
    os.popen('cd-hit -i ' + putative_RE_sequence_fasta + ' -d 0 -c 0.99 -s 0.99 -A 0.99 -g 1 -sc 1 -o ' + putative_RE_sequence_fasta_clst_99).read()
    putative_RE_sequence_fasta_clst_99_clstr = file_prefix + "putative_RE_sequence_fasta_clst_99" + ".clstr"
    putative_RE_sequence_fasta_clst_99_txt = file_prefix + "putative_RE_sequence_fasta_clst_99" + ".clstr.txt"
    putative_RE_sequence_fasta_clst_99_df = file_prefix + "putative_RE_sequence_fasta_clst_99" + ".clstr.df"
    convert_clstr_to_txt(putative_RE_sequence_fasta_clst_99_clstr, putative_RE_sequence_fasta_clst_99_txt)
    convert_txt_to_df(putative_RE_sequence_fasta_clst_99_txt, putative_RE_sequence_fasta_clst_99_df)
    os.popen('rm ' + putative_RE_sequence_fasta_clst_99_txt).read()

    putative_RE_sequence_fasta_clst_99_clst_95 = file_prefix + "putative_RE_sequence_fasta_clst_99_clst_95" + ".fa"
    os.popen('cd-hit -i ' + putative_RE_sequence_fasta_clst_99 + ' -d 0 -c 0.95 -s 0.95 -A 0.95 -g 1 -sc 1 -o ' + putative_RE_sequence_fasta_clst_99_clst_95).read()
    putative_RE_sequence_fasta_clst_99_clst_95_clstr = file_prefix + "putative_RE_sequence_fasta_clst_99_clst_95" + ".clstr"
    putative_RE_sequence_fasta_clst_99_clst_95_txt = file_prefix + "putative_RE_sequence_fasta_clst_99_clst_95" + ".clstr.txt"
    putative_RE_sequence_fasta_clst_99_clst_95_df = file_prefix + "putative_RE_sequence_fasta_clst_99_clst_95" + ".clstr.df"
    convert_clstr_to_txt(putative_RE_sequence_fasta_clst_99_clst_95_clstr, putative_RE_sequence_fasta_clst_99_clst_95_txt)
    convert_txt_to_df(putative_RE_sequence_fasta_clst_99_clst_95_txt, putative_RE_sequence_fasta_clst_99_clst_95_df)
    os.popen('rm ' + putative_RE_sequence_fasta_clst_99_clst_95_txt).read()
    return putative_RE_sequence_fasta_clst_99_clst_95_df

#----------------------------------
"""3. Take Top Clusters, Expand, and Align"""
#----------------------------------
def align_top_clusters(putative_RE_sequence_fasta_clst_99_clst_95_df, tnpB_coords, tnpB_contig_fasta, expansion, cluster_limit):
    os.popen('mkdir ./cluster_alignments')
    os.popen('mkdir ./cluster_alignments/fasta')
    os.popen('mkdir ./cluster_alignments/bed')
    os.popen('mkdir ./cluster_alignments/muscle_alignments')

    clst = pd.read_csv(putative_RE_sequence_fasta_clst_99_clst_95_df, sep='\t')
    full = tnpB_coords
    contig_file = tnpB_contig_fasta
    expansion = int(expansion)
    nt = int(expansion)

    clst = clst[clst['clstr_size'] > 1]
    cluster_list = set(clst['clstr'].tolist())
    cluster_list = sorted(cluster_list)
    cluster_list = range(0,int(cluster_limit))
    print(cluster_list)

    for clster_num in cluster_list:
        df = clst[clst['clstr'] == clster_num]
        df_new = df[['id']]
        df_new['Protein'] = df_new['id'].str.split('::').str[0]
        df_new['Contig'] = df_new['id'].str.split('::').str[1].str.split(':').str[0]
        df_new = df_new[['Protein','Contig']]

        merged = pd.merge(left = df_new, right = full, left_on = ['Contig', 'Protein'], right_on = [0,3], how='left').drop(columns=['Contig', 'Protein'])

        final_bed = merged
        for i in final_bed.index:
            start = final_bed[1][i]
            stop = final_bed[2][i]

            if final_bed[5][i] == "+":
                final_bed[1][i] = stop - 1
                final_bed[2][i] = stop + expansion
            else:
                final_bed[1][i] = max(start - expansion,0)
                final_bed[2][i] = start + 1

        file_name = './cluster_alignments/bed/' + str(nt)  + 'cluster_' + str(clster_num) + '_expanded_coords.bed'
        final_bed = final_bed[abs(final_bed[1] - final_bed[2]) > (0.9*expansion)]
        #print(final_bed)
        final_bed.to_csv(file_name, sep='\t', header=False, index=False)

        fasta_file_name = './cluster_alignments/fasta/' + str(nt)  + 'cluster_' + str(clster_num) + '_RE_' + str(expansion) + 'ex.fa'

        alignment_file_name = './cluster_alignments/muscle_alignments/' +str(nt)  + 'cluster_' + str(clster_num) + '_RE_' + str(expansion) + 'ex_muscle_align.fa'

        print("Aligning Cluster: " + int(clster_num))
        os.popen('bedtools getfasta -s -name -fo ' + fasta_file_name + ' -fi ' + contig_file + ' -bed ' + file_name).read()
        os.popen('muscle -in ' + fasta_file_name + ' -out ' + alignment_file_name).read()

def main():
    tnpB_coords = pd.read_csv(sys.argv[1], sep='\t', header=None)
    expansion_in = int(sys.argv[2])
    expansion_out = int(sys.argv[3])
    file_prefix = sys.argv[4]
    tnpB_contig_fasta = sys.argv[5]
    expansion = sys.argv[6]
    cluster_limit = sys.argv[7]

    putative_RE_coords_bed = get_putative_coordinates(tnpB_coords, expansion_in, expansion_out, file_prefix)
    putative_RE_sequence_fasta_clst_99_clst_95_df = get_sequences_and_cluster(tnpB_contig_fasta, file_prefix, putative_RE_coords_bed)
    align_top_clusters(putative_RE_sequence_fasta_clst_99_clst_95_df, tnpB_coords, tnpB_contig_fasta, expansion, cluster_limit)
