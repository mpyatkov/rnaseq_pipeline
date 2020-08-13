##################################################################################
#Andy Rampersaud, 05.19.15
#This README contains information about running the Read_Length job
##################################################################################
#Goal: 		confirm our expected read length
#Input:		*.fastq.gz for read1 data and *.fastq.gz for read2 data
#Output:	text file summary of read count with corresponding read length(s)
#NOTE: This pipeline allows for processing of samples in parallel. See the setup_Read_Length.sh for details.  
#---------------------------------------------------------------------------------
#The order of script calls:
#1. Read_Length.sh (sources the setup_Read_Length.sh)
#2. Read_Length.qsub is called by Read_Length.sh
#3. Wait until all jobs have completed running
#4. Read_Length_Summary.sh then summarizes the jobs for the data set
#---------------------------------------------------------------------------------
#The template folder "Read_Length" contains the required scripts for running this Read_Length job
#Location of Read_Length:
#/restricted/projectnb/waxmanlab/routines/RNASeq_Scripts
#1) Copy this folder to your corresponding "Scripts" folder in your data set dir
#2) Update the setup_Read_Length.sh as needed
#3) Navigate (cd) to the directory where you have all of your files from step (1) and change permissions by running the following on the terminal:
chmod 700 *.sh
chmod 700 *.qsub
#4) Run the Read_Length.sh:
./Read_Length.sh
#5) As mentioned above, wait until all jobs have completed running. Then run Read_Length_Summary.sh:
./Read_Length_Summary.sh
#This should create a text file summarizing the Read_Length job
#Note:
#If you want to know each job's duration run the following command in the job folder:
grep 'elapsed' *.o*
#---------------------------------------------------------------------------------
#Job output:
#Read_Length_Stats.txt file
#Example below:
#---------------------------------------------------------------------------------
#SAMPLE_ID 	FASTQ_File_Name 	Read_Length(bp) 	Read_Count
#G83_M1	Waxman-TP17_CGATGT_L007_R1_001.fastq	99	37247804
#G83_M1	Waxman-TP17_CGATGT_L007_R2_001.fastq	99	37247804
#G83_M2	Waxman-TP18_ACAGTG_L007_R1_001.fastq	99	39785793
#G83_M2	Waxman-TP18_ACAGTG_L007_R2_001.fastq	99	39785793
#G83_M3	Waxman-TP19_GCCAAT_L007_R1_001.fastq	99	29319981
#G83_M3	Waxman-TP19_GCCAAT_L007_R2_001.fastq	99	29319981
#G83_M4	Waxman-TP20_CAGATC_L007_R1_001.fastq	99	27293627
#G83_M4	Waxman-TP20_CAGATC_L007_R2_001.fastq	99	27293627
#---------------------------------------------------------------------------------
#Read_Length_Summary.sh:
#This script will create the above summary text file and store in a "Stats/Read_Length" folder
#---------------------------------------------------------------------------------
##################################################################################
