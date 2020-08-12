##################################################################################
#Andy Rampersaud, 02.02.16
#This README contains information about running the bamCorrelate job
##################################################################################
#Goal: 		use the bamCorrelate program to evaluate correlation between samples or replicates
#Input:		BAM file(s) of reads
#Output:	Output_bamCorrelate folder (see list below) 
#---------------------------------------------------------------------------------
#The order of script calls:
#1. bamCorrelate.sh (sources the setup_bamCorrelate.sh)
#2. bamCorrelate.qsub is called by bamCorrelate.sh
#3. Wait until all jobs have completed running
#4. Check out the results in the "Output" folder
#---------------------------------------------------------------------------------
#The template folder "bamCorrelate" contains the required scripts for running this bamCorrelate job
#Location of bamCorrelate:
#/restricted/projectnb/waxmanlab/routines/
#---------------------------------------------------------------------------------
#1) Copy this folder to your corresponding "Scripts" folder in your data set dir
#2) Update the setup_bamCorrelate.sh as needed
#3) Input_Samples.txt as needed
#4) Navigate (cd) to the directory where you have all of your files from step (1) and change permissions by running the following on the terminal:
chmod 700 *.sh
chmod 700 *.qsub
#5) Run the bamCorrelate.sh:
./bamCorrelate.sh
#7) As mentioned above, wait until all jobs have completed running. Then check out the results in the "Output" folder
#Note:
#If you want to know each job's duration run the following command in the job folder:
grep 'elapsed' *.o*
#---------------------------------------------------------------------------------
##################################################################################
#Job output:
#A "Output_bamCorrelate_*" folder will be created within the script folder
#The "Output_bamCorrelate_*" folder will contain:
#(1) heatmap.png:
#	a) Image file of the resulting heatmap of all pairwise comparisons between samples
#(1) scatterplot.png:
#	a) Image file of the resulting scatterplot of all pairwise comparisons between samples
#	b) The heatmap is more informative than the scatterplot 
#---------------------------------------------------------------------------------
##################################################################################
