##################################################################################
#Andy Rampersaud, 09.29.15
#This README contains information about running the DiffExp job
##################################################################################
#Goal: 		use R packages to identify differential genes based on RNA-Seq read counts
#Input:		Count file(s) (other files as well)
#Output:	(A) Output_DiffExp folder containing:
#			(1) DiffExp_v2_*.txt (1 file)
#			(2) DiffExp_v2_*_forSEGEXUpload.txt (for both DESeq and EdgeR) (2 files)
#			(3) Up_Genes.txt and Down_Genes.txt (for both DESeq and EdgeR) (4 files)
#			(4) Venn diagrams (*.png files) comparing DESeq and EdgeR DE genes (2 files)
#			(5) *_Venn_Tables.txt: table version of the Venn diagrams (1 file)
#NOTE: This pipeline allows for processing of samples in parallel. See below for details.  
#---------------------------------------------------------------------------------
#The order of script calls:
#1. DiffExp.sh (sources the setup_DiffExp.sh)
#2. DiffExp.qsub is called by DiffExp.sh
#3. Wait until all jobs have completed running
#4. DiffExp_Summary.sh then summarizes the job output
#---------------------------------------------------------------------------------
#The template folder "DiffExp" contains the required scripts for running this DiffExp job
#Location of DiffExp:
#/restricted/projectnb/waxmanlab/routines/RNASeq_Scripts/
#---------------------------------------------------------------------------------
#DiffExp_1 runs with 3 different GTF files:
#DiffExp_1a: Gene body (full exon region) analysis
#DiffExp_1b: Exon_Only_Regions analysis
#DiffExp_1c: Intron_Only_Regions analysis
#---------------------------------------------------------------------------------
#1) Copy this folder to your corresponding "Scripts" folder in your data set dir
#2) Update the setup_DiffExp.sh as needed
#3) Update the Condition_1.txt as needed
#4) Update the Condition_2.txt as needed
#5) Navigate (cd) to the directory where you have all of your files from step (1) and change permissions by running the following on the terminal:
chmod 700 *.sh
chmod 700 *.qsub
#6) Run the DiffExp.sh:
./DiffExp.sh
#7) As mentioned above, wait until all jobs have completed running. Then run DiffExp_Summary.sh:
./DiffExp_Summary.sh
#Note:
#If you want to know each job's duration run the following command in the job folder:
grep 'elapsed' *.o*
#---------------------------------------------------------------------------------
#Replicates per condition
#To handle multiple replicates per condition there are corresponding "Condition_1.txt" and "Condition_2.txt" text files
#The purpose of these text files is to indicate which samples belong to each condition
################################################
#The text file is formatted like the following:
#----------------------------------------------
#Sample_DIR	Sample_ID	Description
#Sample_Waxman-TP17	G83_M1	Male 8wk-pool 1
#Sample_Waxman-TP18	G83_M2	Male 8wk-pool 2
#Sample_Waxman-TP19	G83_M3	Female 8wk-pool 1
#Sample_Waxman-TP20	G83_M4	Female 8wk-pool 2	
#----------------------------------------------
#The 1st column: The Sample_DIR name
#The 2nd column: Waxman Lab Sample_ID 
#The 3rd column: Sample's description 
################################################
#As mentioned in the setup_DiffExp.sh:
#Assumptions for this job to run correctly:
#1. You have already run the TopHat_Paired_End job
#2. You have already run the Extract_Counts_* job(s)
#3. Your data is organized in the following way:
#You have a data set dir such as:
#/projectnb/wax-es/aramp10/G83_Samples
#Within this dir you have sample specific folders such as:
#G83_M1
#G83_M2
#G83_M3
#G83_M4
#Within each sample specific folder you have a *_R1_*.fastq.gz file and *_R2_*.fastq.gz such as:
#Waxman-TP17_CGATGT_L007_R1_001.fastq.gz
#Waxman-TP17_CGATGT_L007_R2_001.fastq.gz
#Within each sample specific folder you have a "tophat2" folder containing output files
#Within each "tophat2" folder you have a "HTSeq" and/or "featureCounts" folder containing output files
#---------------------------------------------------------------------------------
##################################################################################
#Job output:
#A "Output_DiffExp_1*" folder will be created within the script folder
#The "Output_DiffExp_1*" folder will contain:
#(1) DiffExp_v2_*.txt (1 file)
#(2) DiffExp_v2_*_forSEGEXUpload.txt (for both DESeq and EdgeR) (2 files)
#(3) Up_Genes.txt and Down_Genes.txt (for both DESeq and EdgeR) (4 files)
#(4) DiffExp_*_Venn_Diagrams.pdf (1 file)
#---------------------------------------------------------------------------------
#(1) DiffExp_v2_RefSeq_*.txt:
#1) id, which is your gene id
#2) accession numbers according to Refseq (If you see two accession for one gene, it means that that genes have multiple isoforms, multiple entries in ucsc mm9 annotation)
#3) raw counts for each each replicate
#4) normalized counts for each replicate (RPKM)
#5) normalized counts for each replicate (see the paper for how the normalization is done via "Geometric Normalization")
#6) baseMean, baseMeanA, baseMeanB columns from DESeq
#7) foldChange, log2FoldChange, pval, padj/FDR (for both DESeq and EdgeR)
#Please see the DESeq paper for full output columns explanations:
#Anders S and Huber W (2010). “Differential expression analysis for sequence count data.” Genome Biology, 11, pp. R106.
#---------------------------------------------------------------------------------
#(2) DiffExp_v2_RefSeq_*_forSEGEXUpload.txt:
#This file is a subset of the DiffExp_v2_RefSeq_*.txt
#It's specifically used for upload to the SEGEX database
#---------------------------------------------------------------------------------
#(3) Up_Genes.txt and Down_Genes.txt
#These text files are subsets of the DiffExp_v2_RefSeq_*_forSEGEXUpload.txt
#The Up_Genes.txt: significant differential genes with positive fold change (|FC| > 2 and padj < 0.05)
#The Down_Genes.txt: significant differential genes with negative fold change (|FC| > 2 and padj < 0.05)
#These files are meant to spot check that the DiffExp job worked
#These files are sorted by |FC| so the most "differential" genes will be in these text files
#---------------------------------------------------------------------------------
#(4) DiffExp_*_Venn_Diagrams.pdf (1 file)
#This PDF file shows the Venn diagram comparison between the Up_Genes.txt and Down_Genes.txt between R packages
#---------------------------------------------------------------------------------
##################################################################################
#---------------------------------------------------------------------------------
#After running the DiffExp_Summary.sh script:
#A "Summary_Differential_Expression" folder will be generated
#This folder contains 
#(1) SEGEX_Upload_Files folder
#	(a) Combined folder: Combined SEGEX Upload files for both DESeq and EdgeR (2 files)
#	(b) Individual folder: SEGEX Upload files for each counting method (Gene body (full exon region), Exon_Only_Regions, Intron_Only_Regions) (6 files)
#(2) DE_Gene_Counts.txt: 		summarizes the number of DE genes for each comparison (1 file)
#(3) DE_Genes_Venn_R_Package.pdf: 	Montage of Venn diagrams comparing DE genes between R packages (1 file)
#(4) DE_Gene_Venn_Tables.txt:		table version of the Venn diagrams (1 file)
#(5) DE_Genes_Venn_Count_Method.pdf:	Montage of Venn diagrams comparing DE genes between counting methods (1 file)
#---------------------------------------------------------------------------------
##################################################################################
