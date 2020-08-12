#!/bin/bash
##################################################################################
#Andy Rampersaud, 02.22.16
#Adapted from diff_exp scripts by Tisha Melia
##################################################################################
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
#Within each sample specific folder you have a "fastq" folder with:
#Files: *_R1_*.fastq.gz and *_R2_*.fastq.gz such as:
#Waxman-TP17_CGATGT_L007_R1_001.fastq.gz
#Waxman-TP17_CGATGT_L007_R2_001.fastq.gz
#Within each sample specific folder you have a "tophat2" folder containing output files
#Within each "tophat2" folder you have a "HTSeq" folder containing output files
##################################################################################
#Fill in the following information:
##################################################################################
#Information about the conditions being compared
#The differential expression will output the fold change in the following
#format: CONDITION_2/CONDITION_1
#--------------------------------------------------------------------------------
#Note about CONDITION_*_NAME:
#Do not use *spaces* for CONDITION_*_NAME (this will prevent the job from running correctly)
#Please use underscores in place of spaces
#--------------------------------------------------------------------------------
CONDITION_1_NAME=TEMPLATE
CONDITION_2_NAME=TEMPLATE
##################################################################################
#Comparison Number
#If you are doing multiple comparisons within the same dataset, it's helpful to number each comparison
#For example:
#Comparison: STAT5_High_VS_STAT5_Low, Output Folders: Output_DiffExp_1a, Output_DiffExp_1b, Output_DiffExp_1c
#Comparison: Hypox_30min_VS_Hypox_male, Output Folders: Output_DiffExp_2a, Output_DiffExp_2b, Output_DiffExp_2c
#Comparison: Hypox_90min_VS_Hypox_male, Output Folders: Output_DiffExp_3a, Output_DiffExp_3b, Output_DiffExp_3c
COMPAR_NUM=TEMPLATE
##################################################################################
##################################################################################
#Please: DO NOT EDIT CODE BELOW
#Additional variables will be populated from if statements
##################################################################################
##################################################################################
#Couting program used
#Note that the DiffExp.qsub scripts are custom for the counting program used
#No options for this variable (constant value)
COUNT_PROGRAM="featureCounts"
##################################################################################
