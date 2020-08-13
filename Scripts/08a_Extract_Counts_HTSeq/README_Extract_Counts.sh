##################################################################################
#Andy Rampersaud, 08.19.15
#This README contains information about running the Extract_Counts job
##################################################################################
#Goal: 		use htseq-count to count RNA-Seq reads in defined counting regions
#Input:		*_primary.bam, *.gtf file
#Output:	*_HTSeq.out file
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
#Within each of these count folder(s) there will be a corresponding *_HTSeq.out file (two columns: gene symbol and read count).
#---------------------------------------------------------------------------------
#Note:
#htseq-count will only count uniquely mapped reads (what we want)
#Based on:
#http://www-huber.embl.de/users/anders/HTSeq/doc/history.html
#"htseq-count now skips reads that are non-uniquely mapped according to the ‘NH’ optional field"
#My command: "grep -w "NH:i:1"" is most likely what HT-Seq count is doing regardless of which BAM file is being used.
#---------------------------------------------------------------------------------
#Output description:
#http://www-huber.embl.de/users/anders/HTSeq/doc/count.html
#The script outputs a table with counts for each feature, followed by the special counters, which count reads that were not counted for any feature for various reasons. The names of the special counters all start with a double underscore, to facilitate filtering. (Note: The double unscore was absent up to version 0.5.4). The special counters are:

#    __no_feature: reads (or read pairs) which could not be assigned to any feature (set S as described above was empty).
#    __ambiguous: reads (or read pairs) which could have been assigned to more than one feature and hence were not counted for any of these (set S had mroe than one element).
#    __too_low_aQual: reads (or read pairs) which were skipped due to the -a option, see below
#    __not_aligned: reads (or read pairs) in the SAM file without alignment
#    __alignment_not_unique: reads (or read pairs) with more than one reported alignment. These reads are recognized from the NH optional SAM field tag. (If the aligner does not set this field, multiply aligned reads will be counted multiple times, unless they getv filtered out by due to the -a option.)

#Important: The default for strandedness is yes. If your RNA-Seq data has not been made with a strand-specific protocol, this causes half of the reads to be lost. Hence, make sure to set the option --stranded=no unless you have strand-specific data!
#---------------------------------------------------------------------------------
#The Extract_Counts_Summary.sh will create a separate summary table for the special counters mentioned above
#---------------------------------------------------------------------------------
##################################################################################
