# TAM Library Analysis

These scripts are used to analysis NGS data related to the TAM libraries from Meers et al., 2023

## Installation

Use a package manager, i.e. [pip] (https://pip.pypa.io/en/stable/), to install the required packages listed in requirements.txt, if you do not have them already installed.

Ensure that NGS data are in the directory, "1_seq_data". The scripts require the following fastq files (all of which may be downloaded from the SRA. BioProject Accession #PRJNA925099). Details about fastq files may be found in the SRA and in the sample_info.xlsx file, located in the directory, "0_info".

A3270_220921_R1.fastq
A3269_220921_R1.fastq
A3281_220921_R1.fastq

## Usage

Simply execute "main.sh", which will run the three python scripts in the appropriate order for analysis. 

The critical output file may be found in the directory, "6_summary". The .csv files contained in this directory summarize the input and output counts as well as the log2(fold-change) for each TAM library.