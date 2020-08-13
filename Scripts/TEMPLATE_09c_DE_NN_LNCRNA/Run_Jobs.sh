#!/bin/bash -l
##################################################################################
#Andy Rampersaud, 03.11.16
#06.29.2015: Tisha edited to accomodate edgeR results
#This script would be used to run DiffExp.qsub in parallel for different FASTQ files
#Way to run script:
#Usage: ./Run_Jobs.sh
#Example: 
#./Run_Jobs.sh

set -o errexit
set -o pipefail
set -o nounset

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

(
    set +eu
    module load anaconda2
    source activate RNAseq
)
##################################################################################
#Remove *.o files from previous jobs
rm -rf *.o* *.e* Output_* Summary_Differential_Expression

#Source job-specific variables:
source setup_DiffExp.sh

SCRIPT_DIR="$(pwd)"

echo "You are running differential expression for liver lncRNA"
echo "This analysis only makes sense if you have liver samples"
#---------------------------------------------------------------------------------
#Check that each variable prints a value to the terminal:
echo "-----------------------"
echo "Start of variable list:"
echo "-----------------------"
echo "DATASET_DIR:"
echo ${DATASET_DIR}
echo "DATASET_LABEL:"
echo ${DATASET_LABEL}
echo "GTF_FILES_DIR:"
echo ${GTF_FILES_DIR}
echo "CONDITION_1_NAME:"
echo ${CONDITION_1_NAME}
echo "CONDITION_2_NAME:"
echo ${CONDITION_2_NAME}
echo "COMPAR_NUM:"
echo ${COMPAR_NUM}
echo "COUNT_PROGRAM:"
echo ${COUNT_PROGRAM}
echo "SCRIPT_DIR:"
echo ${SCRIPT_DIR}
echo "TIME_LIMIT:"
echo ${TIME_LIMIT}
echo "-----------------------"
echo "End of variable list"
echo "-----------------------"
#---------------------------------------------------------------------------------
echo "-----------------------"
echo "Start of qsub commands:"
echo "-----------------------"
#Gene Annotation file to use
#---------------------------------------------------------------------------------
#For read mapping purposes we want the full gene body information (use the RefSeq_GeneBody.gtf)
#Presence of splice junctions should be based on the full gene structure for any isoform of a gene symbol
#For read counting purposes we can use either 
#1. RefSeq_GeneBody.gtf
#2. Intron_Only_Regions.gtf
#3. Exon_Only_Regions.gtf
#---------------------------------------------------------------------------------
#Start loop over GTF files:
GTF_List="ncRNA_exon_for_counting.gtf ncRNA_genebodies_for_counting.gtf exonic_only_gene_models_ncRNA_for_counting.gtf intronic_only_gene_models_ncRNA_for_counting.gtf"
for GTF_FILE in $GTF_List
do
    #echo ${ANNOTATION_FILE}
    #---------------------------------------------------------------------------------
    ##################################################################################
    ##################################################################################
    #Please: DO NOT EDIT CODE BELOW
    #Additional variables will be populated from if statements
    ##################################################################################
    ##################################################################################
    #Length of counting regions
    #---------------------------------------------------------------------------------
    #Need this dir that contains the various *_Lengths.txt
    #Feel free to use my dir but it's better practice to have a copy in your own DATASET_DIR
    Lengths_DIR="${GTF_FILES_DIR}/lengths"
    #---------------------------------------------------------------------------------
    #The RPKM calculation requires information about the length of counting regions
    #The GENE_LENGTHS_FILE will be determined by the ANNOTATION_FILE
    #------------------------------------------
    #The OUTPUT_PREFIX will be determined by the ANNOTATION_FILE also
    #The "v2" refers to the version 2 of DiffExp data in the SEGEX database
    #Current relevant platforms in the SEGEX database:
    #Platform,Species,"Probe Count","Sequences Loaded"
    #1_Mouse_RNA-Seq_v1_BaseCounts_DEseq_23K-Illumina-RefSeq-genes
    #2_Mouse_RNA-Seq_v2_FPKM_24K-RefSeq-genes,Mouse-mm9,24197,0
    #The "v2" refers to the 2nd platform
    #Need OUTPUT_PREFIX for the differentialAnalysis.R
    #------------------------------------------
    #The COUNT_DIR will be determined by the ANNOTATION_FILE also
    #Need to know which HTSeq output folder stores the count files
    #------------------------------------------
    #Useful to have a DiffExp_Index
    #The DiffExp_Index will be determined by the ANNOTATION_FILE also
    #DiffExp_1a: Gene body (full exon region) analysis
    #DiffExp_1b: Exon_Only_Regions analysis
    #DiffExp_1c: Intron_Only_Regions analysis
    #------------------------------------------
    if [ "${GTF_FILE}" == "ncRNA_exon_for_counting.gtf" ];
    then
        GENE_LENGTHS_FILE="ncRNA_exon_for_counting_lengths.txt"
        OUTPUT_PREFIX="DiffExp_v2_LncRNA_ExonCollapsed"
        COUNT_DIR="LncRNA_Exon_Collapsed_GTF"
        DiffExp_Index="DiffExp_"${COMPAR_NUM}"d"
        COL_SUFFIX="LncRNA_ExonCollapsed"
    fi
    #------------------------------------------
    if [ "${GTF_FILE}" == "ncRNA_genebodies_for_counting.gtf" ];
    then
        GENE_LENGTHS_FILE="ncRNA_genebodies_for_counting_lengths.txt"
        OUTPUT_PREFIX="DiffExp_v2_LncRNA_GeneBody"
        COUNT_DIR="LncRNA_GeneBody_GTF"
        DiffExp_Index="DiffExp_"${COMPAR_NUM}"e"
        COL_SUFFIX="LncRNA_GeneBody"
    fi
    #------------------------------------------
    if [ "${GTF_FILE}" == "exonic_only_gene_models_ncRNA_for_counting.gtf" ];
    then
        GENE_LENGTHS_FILE="exonic_only_gene_models_ncRNA_for_counting_lengths.txt"
        OUTPUT_PREFIX="DiffExp_v2_LncRNA_Exonic_Only"
        COUNT_DIR="LncRNA_Exonic_Only_GTF"
        DiffExp_Index="DiffExp_"${COMPAR_NUM}"f"
        COL_SUFFIX="LncRNA_Exonic_Only"
    fi
    #------------------------------------------
    if [ "${GTF_FILE}" == "intronic_only_gene_models_ncRNA_for_counting.gtf" ];
    then
        GENE_LENGTHS_FILE="intronic_only_gene_models_ncRNA_for_counting_lengths.txt"
        OUTPUT_PREFIX="DiffExp_v2_LncRNA_Intronic_Only"
        COUNT_DIR="LncRNA_Intronic_Only_GTF"
        DiffExp_Index="DiffExp_"${COMPAR_NUM}"g"
        COL_SUFFIX="LncRNA_Intronic_Only"
    fi
    #------------------------------------------
    #---------------------------------------------------------------------------------
    #--------------------------------------------------------------------------------
    #Create a job name that's a function of the folder name:
    # dir_name and job_name are required in the next steps
    dir_name=$(basename $(pwd))
    step_num=$(echo ${dir_name} | cut -d'_' -f 1)
    job_name="Step_${step_num}"

    #Now use arguments in the PBS script call:
    (set -x; qsub -N "${job_name}_${DiffExp_Index}" -P ${PROJECT} -l h_rt=${TIME_LIMIT} DiffExp.qsub ${SCRIPT_DIR} ${DATASET_DIR} ${DATASET_LABEL} ${GTF_FILES_DIR} ${LNCRNA_ANNOTATION_FILE}  ${CONDITION_1_NAME} ${CONDITION_2_NAME} ${Lengths_DIR} ${GENE_LENGTHS_FILE} ${COUNT_DIR} ${OUTPUT_PREFIX} ${DiffExp_Index} ${COL_SUFFIX} ${COUNT_PROGRAM})

    #---------------------------
    #End loop over GTF files:
done
echo "-----------------------"
echo "End of qsub commands"
echo "-----------------------"
##################################################################################
