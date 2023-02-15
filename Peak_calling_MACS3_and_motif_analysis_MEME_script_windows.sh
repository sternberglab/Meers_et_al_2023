#########################################################################
## Florian T. Hoffmann | Last updated: August 2, 2022                  ##
## Columbia University, Sternberg Lab                                  ##
## Title: Peak_calling_MACS3_and_motif_analysis_MEME_script_windows.sh ##
#########################################################################

# Required Command line packages for this script (these need to be installed prior to running this script):
    # MACS3 (had to be installed manually)
    # NOTE: MACS3 required the installation of the package build-essential

# RAW DATA:
    # All my sequencing run data can be found here: cd ~/florianhoffmann/BaseSpace.
    # Sequencing data folders are named as follows: date_lab_sequencer, e.g. 121820_Sternberg_MiniSeq.

# BEFORE YOU START:
    # 0. Windows Ubuntu directories are slightly different compared to MacOS Terminal. For example, instead of using cd /BaseSpace, use either cd Basespace or cd ~/BaseSpace.
    # 1. This script requires the ChIP-seq data analysis (i.e. mapping of reads to a reference genome) to be completed already as it requires '.bam' files as its input.
    # 2. Note, this script does not require the presence of the index ('.bai') file for each '.bam' file.
    # 3. Please update the file path in the command (Step 1.5) below before executing this script, and save the script. Use the file path of the sequencing folder (e.g. '121820_Sternberg_MiniSeq') as the working directory.
    # 4. Please update the name of the input file that will be used for peak calling (step 2.2).
    # 5. Execute this bash file by navigating to the directory of this script and typing 'bash Peak_calling_MACS3_windows.sh'.
    # 6. Interrupt the running script at any time by pressing ctrl + c.
    # 7. You can read up on different peak calling options available online here: https://github.com/macs3-project/MACS. 

###########################################################################################################################################################################################################################################################################################################

#######################
### START OF SCRIPT ###
#######################

##########################################################################
### STEP 1: Print a welcome message and prepare the coding environment ###
##########################################################################

    # Step 1.1: Define the script as a bash script.

    #!/bin/bash   
    
    # Step 1.2: Start a timer for this script.
    # At the end of the script, a command will print how long it has taken to run the full scrip.

    START_TIME=$SECONDS
    
    # Step 1.3: clear the Terminal window.

    clear

    # Step 1.4: Print welcome message.

    echo ""
    echo "Script for MACS3 peak calling on .bam files and"
    echo "MEME motif analysis on peaks called by MACS3"
    echo ""

    # Step 1.5: Set the working directory.
    # Update the directory below every time.
    # The 'echo' command will print the set working directory upon execution of this script.

    cd ~/BaseSpace/121820_Sternberg_MiniSeq_test
    echo "Your working directory is set to:" $PWD

######################################################
### STEP 2: Call peaks on '.bam' files using MACS3 ###
######################################################

    # Step 2.1: Create a directory into which the MACS3 output files will be placed.
    
    mkdir -p $PWD/MACS3_peak_calling

    # Step 2.2: Change the directory for peak calling.

    cd $PWD/filtered_reads_aligned_sorted_bam_bai

    # Step 2.3: Use the 'macs3 callpeak' command to determine the number of peaks in '.bam' files.
    # Use aligned (to the reference genome), sorted and filtered (i.e. uniquely mapping) reads as input. These can be found in the ChIP-seq analysis folder named 'filtered_reads_aligned_sorted_bam_bai'.
    # This command requires only the ChIP '.bam' files to call peaks.
    # The ChIP sample '.bam' file should be referenced after '-t'.
    # The input sample (non-immunoprecipitated) should be referenced after '-c'.
    # Since ChIP and input file are provided as '.bam' files, for the format select '-f BAM'.
    # The genome size of E. coli BL21(DE3) needs to be provided. Here, we use an estimate of '-g 4500000'. The full genome size is 4558953.
    # '--nomodel' is used to set a custom fragment size using --extsize.
    # '--extsize 400' is used to treat all reads with a fixed fragment size of 400. In ChIP-seq NGS library preparation, I size select for ~450 bp fragment size.
    # '-n' indicates the name of the output files.
    # 'B' will generate a control lambda file.
    # '-q 0.05' indicates the q-value (minimum FDR). This is the default value.
    # '--outdir' indicates the directory into which the output will be placed.
    # The main output will be an 'excel' ('.xls') file that can be opened in Excel as well as a peaks_summits BED ('.bed') file that can be opened in any text editor such as Visual Studio Code.

    for i in *aligned_sorted_filtered.bam; do macs3 callpeak -t "$i" -c A4422_aligned_sorted_filtered.bam -f BAM -g 4500000 --nomodel --extsize 400 -n "${i/_aligned_sorted_filtered.bam/_macs3}" -B -q 0.05 --outdir $PWD/../MACS3_peak_calling/"${i/_aligned_sorted_filtered.bam/_macs3_peaks}"; done

#######################################################################
### STEP 3: Run sequence motif analysis on MACS3 results using MEME ###
#######################################################################

    # Step 3.1: Change into the directory that contains the reference genome file.

    cd $PWD/../reference_genome
    echo ""
    echo "Your working directory has been updated to:" $PWD
    echo ""    
    
    # Step 3.2: Create an indexed version of the reference genome.
    # Use the 'samtools faidx' command to create the '.fai' file from the original '.fasta' reference genome file.

    samtools faidx *.fasta

    # Step 3.3: Create a truncated version of the '.fasta.fai' index reference genome file.
    # In the bioinformatics literature, this truncated file is referred to as the 'chrom.size' file.
    # The output genome file is a '.bed' file that is used by bedtools instead of the '.fasta' genome file.
    # The 'chrom.size' file will be placed into to the 'reference_genome' folder.
    
    cut -f 1,2 *.fasta.fai > $PWD/bl21_chromosome_size.bed
    
    # Step 3.4: Create a directory into which the MEME output files will be placed.
    
    mkdir -p $PWD/../MEME_motif_analysis

    # Step 3.5: Change the directory to access the MACS3 peak calling output files.

    cd $PWD/../MACS3_peak_calling

    # Step 3.6: Create a directory into which copies of all MACS3 output files named 'macs3_summits.bed' will be placed.
    # Thus, the original 'macs3_summits.bed' files will not be modified and usable for other analyses.

    mkdir -p $PWD/summits_macs3_200bp_bed

    # Step 3.7: Generate a copy of all MACS3 output files named 'macs3_summits.bed'.

    find . -name '*macs3_summits.bed' -exec cp {} $PWD/summits_macs3_200bp_bed \;

    # Step 3.8: Change the directory to rename the copied 'macs3_summits.bed' files.

    cd $PWD/summits_macs3_200bp_bed
     
    # Step 3.9: Increase the size of the peak summit coordinates to 200 bp.
    # For this, the MACS3 output file 'AXXXX_macs3_summits.bed' (or 'BXXX_macs3_summits.bed' in case of BaseSpace IDs from a Chavez lab run) will be modified.
    # Flanks of 100 bp will be added upstream and downstream of the original peak summit coordinates.
    # The total summit interval covered by the new coordinates will be 200 bp after executing this step.

    for i in *summits.bed; do bedtools slop -i "$i" -g $PWD/../../reference_genome/*.bed -l 100 -r 99 > "${i/_summits.bed/_summits_200bp.bed}"; done

    # Step 3.10: Delete all original 'macs3_summits.bed' files.
    # Only the 'macs3_summits_200bp.bed' files will be kept.

    rm *summits.bed

    # Step 3.11: Extract 200 bp sequences from the reference genome using peak coordinates called by MACS3.
    # Here, the reference genome in '.fasta' format will be used to extract genomic sequences.
    # Peak summit coordinates from the 'macs3_summits_200bp.bed' files will be used for genomic sequence extraction.
    # The output will be a '.fasta' file containing the 200 bp peak summit cooridnates and the corresponding 200 bp nucleotide sequence from the reference genome.    

    for i in *summits_200bp.bed; do bedtools getfasta -fi $PWD/../../reference_genome/*.fasta -bed "$i" -fo $PWD/../../MEME_motif_analysis/"${i/_summits_200bp.bed/_summits_200bp.fasta}"; done

    # Step 3.12: Change the directory to run the MEME motif analysis.   
    
    cd $PWD/../../MEME_motif_analysis

    # Step 3.13: Print a message that the MEME-ChIP will run for a while.
    # In particular, when analyzing many samples simulatenously, this will take a long time.

    echo ""
    echo "The MEME-ChIP analysis is starting."
    echo "Please be patient, this will take a while..."
    echo ""

    # Step 3.14: Run the MEME motif analysis on all peaks called by MACS3.
    # MEME-ChIP will create a separate output folder in the folder 'MEME_motif_analysis' for each Illumina BaseSpace ID.

    for i in *.fasta; do meme-chip "$i" -o $PWD/"${i/_macs3_summits_200bp.fasta/_macs3_peaks_meme_motifs}"; done

    # Step 3.15: Delete the now redundant copies of the 'macs3_summits_200bp.fasta' files.
    # While the command 'meme-chip' was running, a copy of each the 'macs3_summits_200bp.fasta' file was created in the output directory.
    # Therefore, the original copy can be deleted.

    rm *summits_200bp.fasta 

#############################################################
### STEP 4: End of peak calling and motif analysis script ###
#############################################################

    # Step 4.1: Stop the timer and print how long the execution of this script has taken.
    
    ELAPSED_TIME=$(($SECONDS - $START_TIME))

    echo ""
    echo "Runtime of script: $(($ELAPSED_TIME/60)) min $(($ELAPSED_TIME%60)) sec"  
    echo ""

    sleep 0.5

    # Step 4.2: Remind the user to save the terminal output after running this script.

    echo ""
    echo "REMINDER: Please save the Terminal output in a file."
    echo "A folder named 'run_log' can be used for this purpose."
    echo ""

    # Step 4.3: Print end note.

    echo ""
    echo ""
    echo "END OF SCRIPT"
    echo ""
    echo ""

#####################
### END OF SCRIPT ###
#####################

###########################################################################################################################################################################################################################################################################################################