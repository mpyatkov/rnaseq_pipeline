##################################################################################
#Andy Rampersaud, 06.01.15
#This README contains information about running the TopHat_Paired_End job
##################################################################################
#Goal: 		use TopHat2 to map paired end RNA-Seq reads
#Input:		*.fastq.gz for read1 data and *.fastq.gz for read2 data
#Output:	multiple files (main file: *primary_unique.bam) 
#NOTE: This pipeline allows for processing of samples in parallel. See the setup_TopHat_Paired_End.sh for details.  
#---------------------------------------------------------------------------------
#The order of script calls:
#1. TopHat_Paired_End.sh (sources the setup_TopHat_Paired_End.sh)
#2. TopHat_Paired_End.qsub is called by TopHat_Paired_End.sh
#3. Wait until all jobs have completed running
#4. TopHat_Paired_End_Summary.sh then summarizes the jobs for the data set
#---------------------------------------------------------------------------------
#The template folder "TopHat_Paired_End" contains the required scripts for running this TopHat_Paired_End job
#Location of TopHat_Paired_End:
#/restricted/projectnb/waxmanlab/routines/RNASeq_Scripts
#1) Copy this folder to your corresponding "Scripts" folder in your data set dir
#2) Update the setup_TopHat_Paired_End.sh as needed
#3) Navigate (cd) to the directory where you have all of your files from step (1) and change permissions by running the following on the terminal:
chmod 700 *.sh
chmod 700 *.qsub
#4) Run the TopHat_Paired_End.sh:
./TopHat_Paired_End.sh
#5) As mentioned above, wait until all jobs have completed running. Then run TopHat_Paired_End_Summary.sh:
./TopHat_Paired_End_Summary.sh
#This should create a text file summarizing the TopHat_Paired_End job
#Note:
#If you want to know each job's duration run the following command in the job folder:
grep 'elapsed' *.o*
#---------------------------------------------------------------------------------
#Job output:
#A "aligner" folder will be created within each sample specific folder
#Output file descriptions:
#1) accepted_hits.bam -> The mapped reads for your sample. If your reads are mappable to several places, this file will have up to 20 of those.
#2) Several mapping statistics:
#a) statistics_for_all_accepted_reads.txt contains the mapping statistics of all mappable reads, which include multi-mappable reads.
#b) statistics_for_primary_reads.txt contains the mapping statistics of all mappable reads, ***excluding*** multi-mappable reads.
#c) statistics_for_nonprimary_reads.txt contains the mapping statistics for multi-mappable reads.
#d) statistics_for_unmapped_reads.txt contains the mapping statistics for unmapped reads.
#e) numunique_hits.txt contains the mapping statistics of all mappable reads, ***excluding*** multi-mappable reads (identical to b)). This is not applicable for paired-end reads.
#---------------------------------------------------------------------------------
#primary.bam -> The mapped reads that only contain the best mappable place for each of your read that was able to be mapped.
#junctions.bed. A UCSC BED track of junctions reported by TopHat. Each junction consists of two connected BED blocks, where each block is as long as the maximal overhang of any read spanning the junction. The score is the number of alignments spanning the junction.
#insertions.bed and deletions.bed. UCSC BED tracks of insertions and deletions reported by TopHat. 
#---------------------------------------------------------------------------------
#primary_unique.bam -> The subset of the primary.bam where each read only has 1 reported alignment
#In other words, these are the uniquely mapped reads
#The uniquely mapped reads is the read set we want to use for both counting and visualization
#---------------------------------------------------------------------------------
#TopHat_Paired_End_Summary.sh:
#The output of this script is a TopHat2_Stats.txt file that summarizes the TopHat2_Stats for each sample
#It also generates a TopHat2_Stats_SpliceReads.txt file that summarizes the number of reads with skipped regions and splice junction reads
#See TopHat_Paired_End.qsub for additional details
#---------------------------------------------------------------------------------
##################################################################################
