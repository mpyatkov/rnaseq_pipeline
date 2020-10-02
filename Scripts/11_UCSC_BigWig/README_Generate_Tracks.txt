##################################################################################
#Andy Rampersaud, 08.14.15
#This README contains information about running Generate_Tracks
#Generate_Tracks does not submit jobs to the compute nodes
#The output text file should be saved to the local computer to allow upload to the UCSC Browser
##################################################################################
#Goal: 		Create track lines to provide the UCSC Browser to view data
#Input:		Text file tables 
#Output:	*_Tracks.txt (Directly uploaded to the UCSC Browser)
#NOTE: This pipeline allows for processing of samples in parallel. See the setup_Generate_Tracks.sh for details.  
#---------------------------------------------------------------------------------
#The UCSC tracks generated will be using information from two text files:
#(1) BigBed_Name_Color.txt: (Nothing to change/update)
#This text file is used the add the male and female mouse liver chromatin state maps
#More than likely for a typical UCSC Session you'll want the chromatin state maps with your screen shots
#(2) Sample_Labels_Color.txt: (Needs to be changed/updated)
#This text file is used to indicate information about your samples
#Example of Sample_Labels_Color.txt:
#---------------------------------------------------------------------------------
#Sample_DIR	Sample_ID	Description	Color
#Sample_Waxman-TP17	G83_M1	Male 8wk-pool 1	0,0,255
#Sample_Waxman-TP18	G83_M2	Male 8wk-pool 2	0,0,255
#Sample_Waxman-TP19	G83_M3	Female 8wk-pool 1	255,0,0
#Sample_Waxman-TP20	G83_M4	Female 8wk-pool 2	255,0,0
#---------------------------------------------------------------------------------
#The Color will determine the bigWig track color
#The BAM file of read pile-ups will be colored by strand
#---------------------------------------------------------------------------------
#The order of script calls:
#1. Generate_Tracks.sh (sources the setup_Generate_Tracks.sh)
#---------------------------------------------------------------------------------
#The template folder "Generate_Tracks" contains the required scripts for running this Generate_Tracks job
#Location of Generate_Tracks:
#/restricted/projectnb/waxmanlab/routines/RNASeq_Scripts
#---------------------------------------------------------------------------------
#1) Copy this folder to your corresponding "Scripts" folder in your data set dir
#2) Update the setup_Generate_Tracks.sh as needed
#3) Navigate (cd) to the directory where you have all of your files from step (1) and change permissions by running the following on the terminal:
# chmod 700 *.sh
#4) Run the Generate_Tracks.sh:
# ./Generate_Tracks.sh
#---------------------------------------------------------------------------------
##################################################################################
#Job output:
#The output of this job will be a "UCSC_Track_Lines" folder within the "10_Generate_Tracks" folder
#Within the "UCSC_Track_Lines" folder:
#(1) G*_Tracks_PileUp.txt:		Use these track lines to view read pile-ups (BAM files)
#(2) G*_Tracks_Wiggle.txt:		Use these tracks lines to view wiggle tracks (autoScale=off viewLimits=0.0:100.0)
#(3) G*_Tracks_Wiggle_autoON.txt:	Use these tracks lines to view wiggle tracks (autoScale=on)
#My reason for the above 3 sets of track lines:
#Read pile-ups are most informative to see relativly few reads (less informative for read dense regions)
#When viewing wiggle files it's useful to have a session with (autoScale=on) to allow one to record an appropriate y-axis value for all the tracks
#Using this recorded y-axis value, one could then modify the G*_Tracks_Wiggle.txt such that (viewLimits=0.0:100.0) is the recorded y-axis value (do a find and replace all function in a text editor)
#This facilitates screen shots of highly differentially expresed genes
#---------------------------------------------------------------------------------
#Instructions for using track lines:
#Once you have your text files of track lines, save these text files to your local computer
#Go to the UCSC Genome Browser:
#http://genome.ucsc.edu/
#Click on "Genome Browser" -> click on "manage custom tracks" button -> click on "add custom tracks" button -> click "Browse" button -> Find your text file of track lines -> click "Submit" button
#Hopefully all the tracks load successfully (don't see any errors)
#Once you have viewed your tracks, it's a good idea to save it as a session in your UCSC account
#---------------------------------------------------------------------------------
##################################################################################
