##################################################################################
#Andy Rampersaud, 09.28.15
#This README contains information about running the Extract_Counts job
##################################################################################
#Goal: 		use htseq-count to count RNA-Seq reads in defined counting regions
#Input:		*_primary.bam, *.gtf 
#Output:	*_featureCounts.out and featureCounts.out.summary 
#NOTE: This pipeline allows for processing of samples in parallel. See the setup_Extract_Counts.sh for details.  
#---------------------------------------------------------------------------------
#The order of script calls:
#1. Extract_Counts.sh (sources the setup_Extract_Counts.sh)
#2. Extract_Counts.qsub is called by Extract_Counts.sh
#3. Wait until all jobs have completed running
#4. Extract_Counts_Summary.sh then summarizes the jobs for the data set
#---------------------------------------------------------------------------------
#The template folder "Extract_Counts" contains the required scripts for running this Extract_Counts job
#Location of Extract_Counts:
#/restricted/projectnb/waxmanlab/routines/RNASeq_Scripts
#---------------------------------------------------------------------------------
#Extract_Counts runs with 4 different GTF files:
#(1) GTF File: genes.gtf			
#(2) GTF File: RefSeq_GeneBody.gtf
#(3) GTF File: Intron_Only_Regions.gtf
#(4) GTF File: Exon_Only_Regions.gtf
#---------------------------------------------------------------------------------
#1) Copy this folder to your corresponding "Scripts" folder in your data set dir
#2) Update the setup_Extract_Counts.sh as needed
#3) Navigate (cd) to the directory where you have all of your files from step (1) and change permissions by running the following on the terminal:
chmod 700 *.sh
chmod 700 *.qsub
#4) Run the Extract_Counts.sh:
./Extract_Counts.sh
#5) As mentioned above, wait until all jobs have completed running. Then run Extract_Counts_Summary.sh:
./Extract_Counts_Summary.sh
#This should create a text file summarizing the Extract_Counts job
#Note:
#If you want to know each job's duration run the following command in the job folder:
grep 'elapsed' *.o*
#---------------------------------------------------------------------------------
#Job output:
#A count folder will be created within each sample specific folder
#Here is the pairing between counting method and count folder name:
#--------------------------
#(1) GTF File: genes.gtf, Count_DIR_Name: Illumina_GTF			
#(2) GTF File: RefSeq_GeneBody.gtf, Count_DIR_Name: RefSeq_Exon_GTF
#(3) GTF File: Intron_Only_Regions.gtf, Count_DIR_Name: RefSeq_Intron_GTF
#(4) GTF File: Exon_Only_Regions.gtf, Count_DIR_Name: RefSeq_Exon_Only_GTF
#--------------------------
#Within each of these count folder(s) there will be a corresponding *_featureCounts.out file (two columns: gene symbol and read count).
#---------------------------------------------------------------------------------
#Note:
#The featureCount program has the following option:
#---------------------------------------------------------------------------------
#    -M        	If specified, multi-mapping reads/fragments will be counted (ie.
#              	a multi-mapping read will be counted up to N times if it has N
#              	reported mapping locations). The program uses the `NH' tag to
#              	find multi-mapping reads.
#---------------------------------------------------------------------------------
#This (-M) option is *not* being used in the qsub script, therefore multi-mapping reads/fragments will *NOT* be counted
#We only want to count uniquely mapped reads
#---------------------------------------------------------------------------------
#Output description:
#http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf
#See above PDF for more details
#---------------------------------------------------------------------------------
#The qsub script processes the featureCount output file to a more simple output file
#---------------------------------------------------------------------------------
#The Extract_Counts_Summary.sh will create a separate summary table for the special counters mentioned above
#---------------------------------------------------------------------------------
##################################################################################
