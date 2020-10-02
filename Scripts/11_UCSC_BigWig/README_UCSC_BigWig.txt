##################################################################################
#Andy Rampersaud, 08.18.15
#This README contains information about running the UCSC_BigWig job
##################################################################################
#Goal: 		convert RNA-seq data from BAM format into wiggle format for UCSC Browser visualization
#Input:		*_primary_unique.bam 
#Output:	*.bw file(s)
#NOTE: This pipeline allows for processing of samples in parallel. See the setup_UCSC_BigWig.sh for details.  
#---------------------------------------------------------------------------------
#The order of script calls:
#1. UCSC_BigWig.sh (sources the setup_UCSC_BigWig.sh)
#2. UCSC_BigWig.qsub is called by UCSC_BigWig.sh
#3. Wait until all jobs have completed running
#4. UCSC_BigWig_Summary.sh then summarizes the jobs for the data set
#---------------------------------------------------------------------------------
#The template folder "UCSC_BigWig" contains the required scripts for running this UCSC_BigWig job
#Location of UCSC_BigWig:
#/restricted/projectnb/waxmanlab/routines/RNASeq_Scripts
#---------------------------------------------------------------------------------
#1) Copy this folder to your corresponding "Scripts" folder in your data set dir
#2) Update the setup_UCSC_BigWig.sh as needed
#3) Navigate (cd) to the directory where you have all of your files from step (1) and change permissions by running the following on the terminal:
# chmod 700 *.sh
# chmod 700 *.qsub
#4) Run the UCSC_BigWig.sh:
# ./UCSC_BigWig.sh
#5) As mentioned above, wait until all jobs have completed running. Then run UCSC_BigWig_Summary.sh:
# ./UCSC_BigWig_Summary.sh
#This script will copy BAM and *.bw files to the the lab server (waxmanlabvm.bu.edu)
#Note:
#If you want to know each job's duration run the following command in the job folder:
# grep 'elapsed' *.o*
#---------------------------------------------------------------------------------
#Job output:
#An output folder named "UCSC_BigWig" will be saved into your "tophat2" folder
#Note that job log files may indicate the following:
#/bin/sh: wigToBigWig: command not found
#/bin/sh: wigToBigWig: command not found
#Please ignore the above messages and confirm that the "UCSC_BigWig" (saved into your "tophat2" folder) folder exists and that your *.bw file(s) are present
#---------------------------------------------------------------------------------
##################################################################################
