##################################################################################
#Andy Rampersaud, 06.02.15
#This README contains information about running the Read_Strandness job
##################################################################################
#Goal: 		Infer RNA-seq experiment design from SAM/BAM file. (strand specificity)
#Input:		*_primary_unique.bam 
#Output:	*Read_Strandness.txt file
#NOTE: This pipeline allows for processing of samples in parallel. See the setup_Read_Strandness.sh for details.  
#---------------------------------------------------------------------------------
#The order of script calls:
#1. Read_Strandness.sh (sources the setup_Read_Strandness.sh)
#2. Read_Strandness.qsub is called by Read_Strandness.sh
#3. Wait until all jobs have completed running
#4. Read_Strandness_Summary.sh then summarizes the jobs for the data set
#---------------------------------------------------------------------------------
#The template folder "Read_Strandness" contains the required scripts for running this Read_Strandness job
#Location of Read_Strandness:
#/restricted/projectnb/waxmanlab/routines/RNASeq_Scripts
#---------------------------------------------------------------------------------
#1) Copy this folder to your corresponding "Scripts" folder in your data set dir
#2) Update the setup_Read_Strandness.sh as needed
#3) Navigate (cd) to the directory where you have all of your files from step (1) and change permissions by running the following on the terminal:
chmod 700 *.sh
chmod 700 *.qsub
#4) Run the Read_Strandness.sh:
./Read_Strandness.sh
#5) As mentioned above, wait until all jobs have completed running. Then run Read_Strandness_Summary.sh:
./Read_Strandness_Summary.sh
#This should create a text file summarizing the Read_Strandness job
#Note:
#If you want to know each job's duration run the following command in the job folder:
grep 'elapsed' *.o*
#---------------------------------------------------------------------------------
#Job output:
#A text file describing the strand specificity will be created
#You'll need this information before running the UCSC_BigWig
#---------------------------------------------------------------------------------
##################################################################################
