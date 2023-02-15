#########################################################################
## Florian T. Hoffmann | Last updated: December 7, 2022                ##
## Columbia University, Sternberg Lab                                  ##
## Title: Combined_ChIPseq_analysis_script_paired-end_reads_Windows.sh ##
#########################################################################

# Required Command line packages for this script (these need to be installed prior to running this script):
    # fastp (available through conda; bioconda channel)
    # Bowtie 2 (downloaded latest version from http://bowtie-bio.sourceforge.net/bowtie2/index.shtml and placed into the conda packages directory because the conda version did not work)
    # Samtools (installed manually since bioconda did not work)
    # MACS2 (available through conda; bioconda channel)
    # deepTools (available through conda; bioconda channel)
    # NOTE 1: conda packages can be found in this directory: cd /Users/florian/opt/anaconda3/pkgs
    # NOTE 2: If the 'conda' command cannot be executed, add conda to the PATH again: export PATH=/Users/florian/miniconda3/bin:$PATH

# RAW DATA:
    # All my sequencing run data can be found here: cd /Users/florian/BaseSpace.
    # Sequencing data folders are named as follows: date_lab_sequencer, e.g. 121820_Sternberg_MiniSeq.

# BEFORE YOU START: 
    # 0. Windows Ubuntu directories are slightly different compared to MacOS Terminal. For example, instead of using cd /BaseSpace, use either cd Basespace or cd ~/BaseSpace.
    # 1. This script is written for processing paired-end reads and needs to be modified (for fastp, Bowtie 2) for processing single-end reads.
    # 2. This script requires lanes to be merged already (use 'merge_lanes.sh' file if you have not merged the '.fastq.gz' files yet).
    # 3. Rename the folder that contains the .fastq files (e.g. change it to '121820_Sternberg_MiniSeq') (optional).
    # 4. Delete all BaseSpace ID folders that are not yours/that were created due to contamination (all folders that you do not want to have analyzed).
    # 5. Please update the file path in the command (Step 1.7) below before executing this script, and save the script.
    # 6. Please update the directory in which Bowtie 2 was installed (Step 1.5).
    # 7. Create a folder in your working directory (must be same as $PWD/..) named 'reference_genome' and place your reference genome (in '.fasta' format) in there.
    # 8. Please place the 'bedGraphToBigWig' Unix shell file into a folder called 'src'.
    # 9. Make 'bedGraphToBigWig' executable by running 'chmod +x bedGraphToBigWig' before executing this script.
    # 9. Execute this bash file by navigating to the directory of this script and typing 'bash ChIPseq_analysis_script'.
    # 10. Interrupt the running script at any time by pressing ctrl + c.

# Required inputs for this script:
    # 'src' folder: must contain the 'bedGraphToBigWig' file; place this shell script into this folder as well.
    # 'reference genome' folder: must contain your reference genome in '.fasta' format.
    # 'fastQ_files' folder: must contain your raw/unprocessed fastQ files (one file for read 1 and another file for read 2 for paired-end data). Please merge sequencing files from all lanes prior to adding the files.

# Outputs of this script:
    # 'fastQ_files_backup' folder: is a duplicate/backup of all original sequencing files.
    # 'fastQ_files_trimmed' folder: adapter-trimmed sequences.
    # 'reference_genome_indexed' folder: contains processed files made from the original reference genome that are required for read mapping by bowtie2.
    # 'reads_aligned_sam' folder: mapped (=aligned) read files in '.sam' format.
    # 'reads_aligned_bam' folder: mapped (=aligned) read files in '.bam' format.
    # 'unfiltered_reads_aligned_sorted_bam_bai' folder: mapped and sorted reads in '.bam' and '.bai' formats.
    # 'filtered_reads_aligned_sorted_bam_bai' folder: mapped, sorted and uniquely mapping (=filtered) reads in '.bam' and '.bai' formats.
    # 'normalization' folder: contains all normalized read files.
    # 'run_log' folder: save the output of this script as a '.txt' file and put it into this folder.

###########################################################################################################################################################################################################################################################################################################

#######################
### START OF SCRIPT ###
#######################

#################################################################################################
### STEP 1: Print a welcome message, prepare the coding environment and command line packages ###
#################################################################################################

    # Step 1.1: Define the script as a bash script.

    #!/bin/bash
    
    # Step 1.2: Start a timer for this script.
    # At the end of the script, a command will print how long it has taken to run the full scrip.

    START_TIME=$SECONDS
    
    # Step 1.3: clear the Terminal window.

    clear

    # Step 1.4: 'Title of script' note.
    # First, print a welcome message stating the title and author of this shell script.
    
    echo ""
    echo ""
    echo ""
    echo "Welcome to the ChIP-seq analysis Unix shell script."
    echo "Read mapping, quality filtering, RPKM normalization and" 
    echo "max. value per 1 kb plotting will be performed."
    echo ""
    echo "This script was written by Florian T. Hoffmann." 
    echo ""
    echo ""
    echo ""
    
    # Step 1.5: Set the working directory.
    # Update the directory below every time.
    # The 'echo' command will print the set working directory upon execution of this script.

    cd ~/BaseSpace/Test
    echo ""
    echo "Your working directory is set to:" $PWD
    echo ""

    # Step 1.6: Create a directory into which a log file of all Terminal outputs will be placed.

    mkdir -p $PWD/run_log

############################################
### STEP 2: Prepare the 'fastq.gz' files ###
############################################

    # Step 2.1: Duplicate the 'fastQ_files' folder and its contents, so there is a back-up.
    # (In case the code messes up the files.)

    cp -rp $PWD/fastQ_files/. $PWD/fastQ_files_backup/

    # Step 2.2: Set the working directory to the 'fastQ_file' folder.
    # Upon execution of this script, the 'echo' command will print that working directory has been changed in this step.

    cd $PWD/fastQ_files
    echo ""
    echo "Your working directory has been updated to:" $PWD
    echo ""

    # Step 2.3: Each '.fastq.gz' file is packaged into a folder. This folder is redundant for the analysis.
    # Thus, this folder (once empty) can be deleted later.
    # Extract the '.fastq.gz' files from that folder and place them into the parent folder.

    find . -name '*.fastq.gz' -exec mv {} $PWD/ \;

    # Step 2.4: The folder is empty now.
    # Delete the empty folder.

    find . -depth -type d -empty -exec rmdir {} \;

    # Step 2.5: The commands below shorten '.fastq.gz' file names of the .fastq files.
    # Only the BaseSpace ID (Starting with 'A4...') and the read ID are kept as the file name.
    # (In the paired-end mode, two reads are generated for each sample, so one has the ID '1' the other '2'.)
    # For example, a file named 'A4414_S19_L001_R1_001.fastq.gz' will be converted into a file named 'A4414_1'.
    # The new shortened file names will be printed using the 'echo' command.

    for i in *.fastq.gz; do mv "$i" "${i/_S?_L001_R1_001/_1}"; done
    for i in *.fastq.gz; do mv "$i" "${i/_S?_L001_R2_001/_2}"; done
    for i in *.fastq.gz; do mv "$i" "${i/_S??_L001_R1_001/_1}"; done
    for i in *.fastq.gz; do mv "$i" "${i/_S??_L001_R2_001/_2}"; done
    for i in *.fastq.gz; do mv "$i" "${i/_S???_L001_R1_001/_1}"; done
    for i in *.fastq.gz; do mv "$i" "${i/_S???_L001_R2_001/_2}"; done

    for i in *.fastq.gz; do echo "$i"; done
    echo ""

#############################################################
### STEP 3: Trim and quality-filter the reads using fastp ###
#############################################################

    # Step 3.1: Create a new folder into which trimmed reads will be saved.
    # The folder will be named 'fastQ_files_trimmed'

    mkdir -p $PWD/../fastQ_files_trimmed

    # Step 3.2: Use fastp to trim and quality filter the reads.
    # Use default parameters.
    # When using paired-end reads, use -i and -o for the first read.
    # And use -I and -O for the second read.
    # When usign single-end reads, only use -i and -o and delete -I and -O.

    for i in *1.fastq.gz; do fastp -i "$i" -I "${i/_1/_2}" -o "${i/_1/_1_trimmed}" -O "${i/_1/_2_trimmed}" -j "${i/_1.fastq.gz/_trimmed}".json -h "${i/_1.fastq.gz/_trimmed}".html; done

    # Step 3.3: Move the trimmed reads generated above (Step 2.2) into the 'fastQ_files_trimmed' folder.
    # Step 3.2 also created summary files: '_trimmed.html' and '_trimmed.json'.
    # This command also moves the summary files into the same folder due to similar naming.

    mv $PWD/*trimmed* $PWD/../fastQ_files_trimmed

#########################################################################################
### STEP 4: Index reference genome and align reads to reference genome using Bowtie 2 ###
#########################################################################################

    # Step 4.1: Create a directory into which the indexed reference genome will be placed.

    mkdir -p $PWD/../reference_genome_indexed

    # Step 4.2: Change into the newly created directory. 
    # Bowtie 2 will deposit the indexed reference genome in the directory that you are currently in.
    # Upon execution of this script, the 'echo' command will print that working directory has been changed in this step.

    cd $PWD/../reference_genome_indexed
    echo ""
    echo "Your working directory has been updated to:" $PWD
    echo ""

    # Step 4.3: Index the reference genome.
    # This command generates 6 output files that end in '.bt2'.
    # The first part (starting with '$BT2_HOME/') of this command locates the file encoding the 'bowtie2-build' command.
    # The second part (starting with '$PWD/') of this command locates the unindexed reference genome. '*.fasta' indicates that any file ending in '.fasta' will be used.
    # The third/last part (strating with 'indexed') of this command determines the names of the indexing output files.
   
    bowtie2-build $PWD/../reference_genome/*.fasta indexed_genome

    # Step 4.4: Create a directory into which the aligned reads will be placed.

    mkdir -p $PWD/../reads_aligned_sam

    # Step 4.5: Change into the 'fastQ_files_trimmed' directory.
    # Bowtie 2 will deposit the aligned reads into this directory in the next step.
    # Upon execution of this script, the 'echo' command will print that working directory has been changed in this step.

    cd $PWD/../fastQ_files_trimmed
    echo ""
    echo "Your working directory has been updated to:" $PWD
    echo ""

    # Step 4.6: Align paired-end reads to the indexed reference genome.
    # This command generates '.sam' output files. 
    # The first part (starting with '*1_trimmed.fastq.gz') of this command marks all read '1' files relevant for the 'for loop'.
    # The second part (starting with 'echo') of this command prints what file is being processed.
    # The third part (starting with '$BT2_HOME/') of this command locates the file encoding the 'bowtie2' command.
    # The fourth part (starting with '$PWD/../reference') of this command locates the folder containing the indexed reference genome files.
    # The fifth part (starting with '-1 $PWD/../fastQ') of this command locates read '1' of the paired-end reads.
    # The sixth part (starting with '-2 $PWD/../fastQ') of this command locates read '2' of the paired-end reads.
    # The seventh/last part (starting with '-S') of this command determines the names of the aligning output files and puts them into the '$PWD/../reads_aligned' directory.
    
    for i in *1_trimmed.fastq.gz; do echo ""; echo ""; echo "The file $i is being aligned..."; echo ""; echo ""; bowtie2 -x $PWD/../reference_genome_indexed/indexed_genome -1 "$i" -2 "${i/_1/_2}" -S $PWD/../reads_aligned_sam/"${i/_1_trimmed.fastq.gz/_aligned}".sam; done

##############################################################################################################################
### STEP 5: Convert '.sam' files into '.bam' files, sort and index files, and eliminate multi-mapping reads using Samtools ###
##############################################################################################################################

    # Step 5.1: Create a directory into which the '.bam' files will be placed.

    mkdir -p $PWD/../reads_aligned_bam

    # Step 5.2: Change into the directory that contains the '.sam' files.
    # Samtools will use '.sam' files in this directory and generate '.bam' files in the next step.
    # Upon execution of this script, the 'echo' command will print that working directory has been changed in this step.

    cd $PWD/../reads_aligned_sam
    echo ""
    echo "Your working directory has been updated to:" $PWD
    echo ""

    # Step 5.3: Convert '.sam' files into '.bam' files using the 'samtools view' command.
    # The '-b' option is used to produce '.bam' output files.
    # This step is required since IGV does not accept '.sam' input files.
    
    for i in *aligned.sam; do samtools view -b "$i" > $PWD/../reads_aligned_bam/"${i/_aligned.sam/_aligned.bam}"; done

    # Step 5.4: Create a directory into which the sorted '.bam' and index '.bai' files will be placed.
    # This directory will later contain the unfiltered (i.e. multi-mapping) '.bam' reads and their corresponding index '.bai' files.

    mkdir -p $PWD/../unfiltered_reads_aligned_sorted_bam_bai
    
    # Step 5.5: Change into the directory that contains the unsorted '.bam' files.
    # Samtools will use files in this directory in the next step.
    # Upon execution of this script, the 'echo' command will print that working directory has been changed in this step.

    cd $PWD/../reads_aligned_bam
    echo ""
    echo "Your working directory has been updated to:" $PWD
    echo ""

    # Step 5.6: Sort the '.bam' files using the 'samtools sort' command.
    # By default, this command will create '.bam' output files because the input files are in the '.bam' format.
    # This command also creates temporary '.tmp.0000.bam' files that will get deleted automatically prior to completion of this command.
    
    for i in *aligned.bam; do samtools sort -o $PWD/../unfiltered_reads_aligned_sorted_bam_bai/"${i/_aligned/_aligned_sorted}" "$i"; done 

    # Step 5.7: Change into directory that contains the aligned and sorted but unfiltered '.bam' files.
    # Upon execution of this script, the 'echo' command will print that working directory has been changed in this step.

    cd $PWD/../unfiltered_reads_aligned_sorted_bam_bai
    echo ""
    echo "Your working directory has been updated to:" $PWD
    echo ""

    # Step 5.8: Index the aligned and sorted but unfiltered '.bam' files using the 'samtools index' command.
    # This command generates '.bai' output files that are required to exist in order fro IGV to visualize '.bam' files.
    # The '.bai' files will be placed into the same directory, so that they can be easily found later by IGV (which searches for '.bai' files associated with '.bam' files in the same directory).
    # The '-b' option is used to produce '.bai' output files.

    for i in *aligned_sorted.bam; do samtools index -b "$i" "${i/_aligned_sorted.bam/_aligned_sorted.bam.bai}"; done

    # Step 5.9: Create a directory into which the aligned, sorted and filtered (uniquely mapping) '.bam' files and their corresponding '.bai' files will be placed.

    mkdir -p $PWD/../filtered_reads_aligned_sorted_bam_bai
    
    # Step 5.10: Eliminate multi-mapping reads (filtering) using the 'samtools view' command.
    # This command will eliminate all multi-mapping reads and will only retain uniquely-mapping reads.
    # The '-q' option is set to '10', so that all reads with a MAPQ score < 10 will be eliminated. Only reads with a MAPQ score >= 10 will be retained.
    # The '-b' option is used to produce '.bam' output files.

    for i in *aligned_sorted.bam; do samtools view -bq 10 "$i" > $PWD/../filtered_reads_aligned_sorted_bam_bai/"${i/_aligned_sorted.bam/_aligned_sorted_filtered.bam}"; done

    # Step 5.11: Change into the directory that contains the 'aligned_sorted_filtered.bam' reads.
    # Samtools will use files in this directory to create their index '.bai' files in the next step.
    # Upon execution of this script, the 'echo' command will print that working directory has been changed in this step.

    cd $PWD/../filtered_reads_aligned_sorted_bam_bai
    echo ""
    echo "Your working directory has been updated to:" $PWD
    echo ""

    # Step 5.12: Print the number of uniquely mapping reads of each final '.bam' file.
    # Use the 'samtools view -c' command.
    # Use "echo" to also print the name of each file.

    echo ""
    echo "The following uniquely mapping read files have been generated:"
    echo ""

    for i in *aligned_sorted_filtered.bam; do echo ""; echo "Number of uniquely mapping reads for $i:"; samtools view -c $i; done
    
    # Step 5.13: Print the directory of the uniquely mapping reads.
 
    echo ""
    echo "Directory of uniquely-mapping reads:" $PWD
    echo ""

    # Step 5.14: Create index '.bai' files for the unique (filtered), aligned and sorted reads.
    # This is required in order for IGV to accept '.bam' files.
    # The '.bai' output files will be placed into the same folder in which the 'aligned_sorted_filtered.bam' files were deposited, so that IGV can find the '.bai' file for each '.bam' file easily.

    for i in *aligned_sorted_filtered.bam; do samtools index -b "$i" "${i/_aligned_sorted_filtered.bam/_aligned_sorted_filtered.bam.bai}"; done

##################################################################
### STEP 6: Delete all files redundant for downstream analysis ###
##################################################################

    # Step 6.1: Print a message telling the user that all intermediate/redundant folders will be deleted.
    # Print the list of folders that will be deleted.

    echo ""
    echo "To reduce memory usage of your local computer, the following folders will be deleted:"
    echo ""
    echo "1. 'fastQ_files'         | file size:" 
    du -hs $PWD/../fastQ_files 
    echo ""
    echo "2. 'fastQ_files_trimmed' | file size:" 
    du -hs $PWD/../fastQ_files_trimmed
    echo ""
    echo "3. 'reads_aligned_sam    | file size:" 
    du -hs $PWD/../reads_aligned_sam
    echo ""
    echo "4. 'reads_aligned_bam    | file size:" 
    du -hs $PWD/../reads_aligned_bam
    echo ""

    # Step 6.2: Delete the folders including their contents.
    # To disable this step, simply add '#' in front of every "rm -r" command.

    #rm -r $PWD/../fastQ_files
    #rm -r $PWD/../fastQ_files_trimmed
    #rm -r $PWD/../reads_aligned_sam
    #rm -r $PWD/../reads_aligned_bam

######################################################
### STEP 7: Normalize '.bam' files using deepTools ###   
######################################################    
    
    # Step 7.1: Create a master directory in which normalization will take place.

    mkdir -p $PWD/../normalization

    # Step 7.2: Create a directory into which the RPKM normalized '.bw' and '.bed' files will be placed.
    
    mkdir -p $PWD/../normalization/bw_bamCoverage-normalized-RPKM_filtered_aligned_sorted_reads
    mkdir -p $PWD/../normalization/bed_bamCoverage-normalized-RPKM_filtered_aligned_sorted_reads

    # Step 7.3: Use the 'bamCoverage' command to RPKM-normalize your '.bam' files into '.bw' files.
    # Use aligned (to the reference genome), sorted and filtered (i.e. uniquely mapping) reads as input.
    # This command requires only the ChIP '.bam' files to normalize according to genome-wide reads.
    # The ChIP sample '.bam' file should be referenced after '-b1'.
    # The output will be a '.bigwig' ('.bw') file that can be imported into IGV for visualization.
    # We normalize using 'reads per kilobase per million mapped reads (RPKM)'. Usage of this parameter DOES NOT require the effective (i.e. mappable) genome size as input.
    # The output will be a '.bigwig' ('.bw') file that can be imported into IGV for visualization. 
    # The bin size is set to '1' using '-bs 1'. This bin size can be used because the E. coli genome is small and does not require much computational power.

    for i in *aligned_sorted_filtered.bam; do bamCoverage --normalizeUsing RPKM -bs 1 -b "$i" -o $PWD/../normalization/bw_bamCoverage-normalized-RPKM_filtered_aligned_sorted_reads/"${i/_filtered.bam/_filtered_normalized-bamCoverage-RPKM.bw}"; done

    # Step 7.4: Use the same 'bamCoverage' command from the step above to RPKM-normalize your '.bam' files into '.bed' files.
    # The output will be a bedgraph ('.bed') file that needs to be converted into a '.bw' file for visualization in IGV. 
    # We normalize using 'reads per kilobase per million mapped reads (RPKM)'. Usage of this parameter DOES NOT require the effective (i.e. mappable) genome size as input.
    # The output will be a '.bigwig' ('.bw') file that can be imported into IGV for visualization. 
    # The bin size is set to '1' using '-bs 1'. This bin size can be used because the E. coli genome is small and does not require much computational power.

    for i in *aligned_sorted_filtered.bam; do bamCoverage --normalizeUsing RPKM --outFileFormat bedgraph -bs 1 -b "$i" -o $PWD/../normalization/bed_bamCoverage-normalized-RPKM_filtered_aligned_sorted_reads/"${i/_filtered.bam/_filtered_normalized-bamCoverage-RPKM.bed}"; done

    # Step 7.5: Change into the normalized file directory and output a list of all normalized files.

    cd $PWD/../normalization/bw_bamCoverage-normalized-RPKM_filtered_aligned_sorted_reads
    echo ""
    echo "Your working directory has been updated to:" $PWD
    echo ""

    echo ""
    echo "The following RPKM-normalized files have been generated:"
    echo ""

    for i in *.bw; do echo $i; done

########################################################################################
### STEP 8: Generate files showing the maximum value in 1 kb tiles across the genome ###   
########################################################################################

    # Step 8.1: Change into the directory that contains the reference genome file.

    cd $PWD/../../reference_genome
    echo ""
    echo "Your working directory has been updated to:" $PWD
    echo ""    
    
    # Step 8.2: Create an indexed version of the reference genome.
    # Use the 'samtools faidx' command to create the '.fai' file from the original '.fasta' reference genome file.

    samtools faidx *.fasta

    # Step 8.3: Create a directory into which the 1 kb max. value pipeline will take place.
    
    mkdir -p $PWD/../normalization/RPKM_max_values

    # Step 8.4: Create all other required accessory folders within the 'RPKM_max_values' directory.

    mkdir -p $PWD/../normalization/RPKM_max_values/support_files
    mkdir -p $PWD/../normalization/RPKM_max_values/bed_RPKM_max_values
    mkdir -p $PWD/../normalization/RPKM_max_values/bw_RPKM_max_values
    
    # Step 8.5: Create a truncated version of the '.fasta.fai' index reference genome file.
    # In the bioinformatics literature, this truncated file is referred to as the 'chrom.size' file.
    
    cut -f 1,2 *.fasta.fai > $PWD/../normalization/RPKM_max_values/support_files/bl21_chromosome_size.bed
    
    # Step 8.6: Change into the directory that contains the truncated reference genome '.bed' file. 
    
    cd $PWD/../normalization/RPKM_max_values/support_files
    echo ""
    echo "Your working directory has been updated to:" $PWD
    echo ""

    # Step 8.7: Split up the reference genome into segments of 1 kb.
    # Use the 'bedtools makewindows' command.

    bedtools makewindows -g bl21_chromosome_size.bed -w 1000 > bl21_1000bp_windows.bed

    # Step 8.8: Change into the directory that contains the RPKM-normalized '.bed' files. 
    
    cd $PWD/../../bed_bamCoverage-normalized-RPKM_filtered_aligned_sorted_reads
    echo ""
    echo "Your working directory has been updated to:" $PWD
    echo ""
    
    # Step 8.9: Using a RPKM-normalized '.bed' file (not .bw), take the max. value for each of the windows generated above.
    # Use the 'bedtools map' command.

    for i in *.bed; do bedtools map -a $PWD/../RPKM_max_values/support_files/bl21_1000bp_windows.bed -b $i -o max -c 4 > $PWD/../RPKM_max_values/bed_RPKM_max_values/"${i/_filtered_normalized-bamCoverage-RPKM.bed/_filtered_normalized-bamCoverage-RPKM_1000bp_max.bed}"; done

    # Step 8.10: Change into the directory that contains the 1000 bp window '.bed' files. 
    
    cd $PWD/../RPKM_max_values/bed_RPKM_max_values
    echo ""
    echo "Your working directory has been updated to:" $PWD
    echo ""

    # Step 8.11: Make 'bedGraphToBigWig' executable.

    chmod +x $PWD/../../../src/bedGraphToBigWig
    
    # Step 8.11: Convert the bedgraph ('.bed') files into '.bw' files.
    # This step is required to visualize the file.

    for i in *.bed; do $PWD/../../../src/./bedGraphToBigWig "$i" $PWD/../support_files/bl21_chromosome_size.bed $PWD/../bw_RPKM_max_values/"${i/_1000bp_max.bed/_1000bp_windows_max.bw}"; done

    # Step 8.12: Change into the normalized file directory and output a list of all 1000-bp-maximum-value-normalized files.

    cd $PWD/../bw_RPKM_max_values
    echo ""
    echo "Your working directory has been updated to:" $PWD
    echo ""

    echo ""
    echo "The following RPKM-normalized files (showing the maximum value per 1 kb window) have been generated:"
    echo ""

    for i in *.bw; do echo $i; done

    # Step 8.13: Print messages how the 1000-bp-maximum-value-normalized files can be visualized in IGV.

    echo ""
    echo "To visualize the 1000-bp-maximum-value-normalized files, do the following:"
    echo ""
    echo "Drag the 'AXXXX_XXXXXX_aligned_sorted_filtered_normalized-bamCoverage-RPKM_1000bp_windows_max.bw' files into IGV"
    echo ""
    echo "Make sure you use the correct reference genome (the same reference genome that you used for read mapping."
    echo ""
    echo "First, the scaling will not be shown correctly in IGV. To fix this, enable 'autoscale'." 
    echo "Then, click the '+' button once (top right corner) to zoom in one level (to the second largest zoom level)."
    echo "From there, set the coordinates in the top bar manually to 1-4,558,959 (for the standard BL21(DE3) reference genome)." 
    echo "This ensures that IGV will not perform any automatic binning (but in the fully zoomed-out view, it will show the incorrect binning)"

#########################################
### STEP 9: Concluding step of script ###   
#########################################

    # Step 9.1: Print a note saying that the '.bam' files should be binned in IGV (using 'count' option) after finishing this script (if visualizing the '.bam' files is wanted).
    
    echo ""
    echo ""
    echo "If you would like to visualize the '.bam' files, please manually convert the '.bam' files of the filtered reads into '.tdf' files using the IGV option 'count'." 
    echo ""
    echo "To find 'count', open IGV, navigate to the tab 'Tools' and select 'Run igvtools...'."
    echo ""
    echo "Set the 'Window Size' to 1."
    echo ""
    echo ""

    # Step 9.2: Stop the timer and print how long the execution of this script has taken.

    ELAPSED_TIME=$(($SECONDS - $START_TIME))

    echo "Runtime of this analysis: $(($ELAPSED_TIME/60)) min $(($ELAPSED_TIME%60)) sec"  

    sleep 0.25

    # Step 7.5: Remind the user to save the terminal output after running this script.

    echo ""
    echo "REMINDER: Please save the Terminal output in a file."
    echo "The file will contain the read filtering details and error messages."
    echo "A folder named 'run_log' has been generated for this purpose."
    echo ""

    # Step 7.6: Print end note.

    echo ""
    echo ""
    echo "END OF SCRIPT"
    echo ""
    echo ""

#####################
### END OF SCRIPT ###
#####################

###########################################################################################################################################################################################################################################################################################################
