##################################################################################
#Andy Rampersaud, 06.05.15
#This README contains information about running Rename_Folders
#Rename_Folders does not submit jobs to the compute nodes
#Result: Sample specific folders will have recognizable names
##################################################################################
#The sequencing facility will typically provide raw data files in the following folder structure:
#---------------------------------------------------------------------------------
#Sample_Waxman-TP17	
#Sample_Waxman-TP18	
#Sample_Waxman-TP19	
#Sample_Waxman-TP20	
#---------------------------------------------------------------------------------
#Experimental information for each sample is usually in an email or lab ChIPSeq_samples_*.xlsx file
#After collecting this information, we simply want to rename these folders so that we have the following:
#---------------------------------------------------------------------------------
#G83_M1	
#G83_M2	
#G83_M3	
#G83_M4	
#---------------------------------------------------------------------------------
#The above is much easier to understand and keep track of progress moving through the analysis pipeline
#---------------------------------------------------------------------------------
#The order of script calls:
#1. Rename_Folders.sh (sources the setup_Rename_Folders.sh)
#---------------------------------------------------------------------------------
#The template folder "Rename_Folders" contains the required scripts for running this Rename_Folders job
#Location of Rename_Folders:
#/restricted/projectnb/waxmanlab/routines/RNASeq_Scripts
#---------------------------------------------------------------------------------
#1) Copy this folder to your corresponding "Scripts" folder in your data set dir
#2) Update the setup_Rename_Folders.sh as needed
#3) Navigate (cd) to the directory where you have all of your files from step (1) and change permissions by running the following on the terminal:
chmod 700 *.sh
#4) Run the Rename_Folders.sh:
./Rename_Folders.sh
#---------------------------------------------------------------------------------
#Job output:
#Check that folder have been properly renamed
#---------------------------------------------------------------------------------
##################################################################################
