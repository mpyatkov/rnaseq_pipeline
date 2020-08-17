#!/bin/bash -l
##################################################################################
#Andy Rampersaud, 03.11.16
#06.29.2015: Tisha edited to accomodate edgeR results
#This script would be used to run DiffExp.qsub in parallel for different FASTQ files
#Way to run script:
#Usage: ./Run_Jobs.sh
#Example: 
#./Run_Jobs.sh
##################################################################################

set -o errexit
set -o pipefail
set -o nounset

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

set +eu
module load anaconda2
source activate RNAseq
set -eu

#Remove *.o files from previous jobs
#Need the 2>/dev/null to supress the "No such file or directory"
rm -rf *.o* *.e* Output_* Summary_Differential_Expression

#Source job-specific variables:
source setup_DiffExp.sh

SCRIPT_DIR="$(pwd)"
#Gene Annotation file to use

#For read mapping purposes we want the full gene body information (use the
#RefSeq_GeneBody.gtf)
#Presence of splice junctions should be based on the full gene structure for any
#isoform of a gene symbol
#For read counting purposes we can use either 
#1. RefSeq_GeneBody.gtf
#2. Intron_Only_Regions.gtf
#3. Exon_Only_Regions.gtf

#Start loop over GTF files:
GTF_List="RefSeq_GeneBody.gtf Exon_Only_Regions.gtf Intron_Only_Regions.gtf"
for GTF_file in $GTF_List
do
    ANNOTATION_FILE=$(basename $GTF_file)

    #Length of counting regions
    
    #Need this dir that contains the various *_Lengths.txt
    #Feel free to use my dir but it's better practice to have a copy in your own Dataset_DIR
    Lengths_DIR="${GTF_FILES_DIR}/lengths"
    
    #The RPKM calculation requires information about the length of counting regions
    #The GENE_LENGTHS_FILE will be determined by the ANNOTATION_FILE
    
    #The OUTPUT_PREFIX will be determined by the ANNOTATION_FILE also
    #The "v2" refers to the version 2 of DiffExp data in the SEGEX database
    #Current relevant platforms in the SEGEX database:
    #Platform,Species,"Probe Count","Sequences Loaded"
    #1_Mouse_RNA-Seq_v1_BaseCounts_DEseq_23K-Illumina-RefSeq-genes
    #2_Mouse_RNA-Seq_v2_FPKM_24K-RefSeq-genes,Mouse-mm9,24197,0
    #The "v2" refers to the 2nd platform
    #Need OUTPUT_PREFIX for the differentialAnalysis.R
    
    #The COUNT_DIR will be determined by the ANNOTATION_FILE also
    #Need to know which HTSeq output folder stores the count files
    
    #Useful to have a DiffExp_Index
    #The DiffExp_Index will be determined by the ANNOTATION_FILE also
    #DiffExp_1a: Gene body (full exon region) analysis
    #DiffExp_1b: Exon_Only_Regions analysis
    #DiffExp_1c: Intron_Only_Regions analysis
    
    if [ "${ANNOTATION_FILE}" == "RefSeq_GeneBody.gtf" ];
    then
        GENE_LENGTHS_FILE="Exon_Regions_Lengths.txt"
        OUTPUT_PREFIX="DiffExp_v2_GeneBody"
        COUNT_DIR="RefSeq_Exon_GTF"
        DiffExp_Index="DiffExp_"${COMPAR_NUM}"a"
        COL_SUFFIX="GeneBody"
    fi
    
    if [ "${ANNOTATION_FILE}" == "Exon_Only_Regions.gtf" ];
    then
        GENE_LENGTHS_FILE="Exon_Only_Regions_Lengths.txt"
        OUTPUT_PREFIX="DiffExp_v2_Exonic_Only"
        COUNT_DIR="RefSeq_Exon_Only_GTF"
        DiffExp_Index="DiffExp_"${COMPAR_NUM}"b"
        COL_SUFFIX="Exonic_Only"
    fi
    
    if [ "${ANNOTATION_FILE}" == "Intron_Only_Regions.gtf" ];
    then
        #Redefine GTF file because we need all gene symbols to be present:
        ANNOTATION_FILE="RefSeq_GeneBody.gtf"
        GENE_LENGTHS_FILE="Intron_Only_Regions_Lengths.txt"
        OUTPUT_PREFIX="DiffExp_v2_Intronic_Only"
        COUNT_DIR="RefSeq_Intron_GTF"
        DiffExp_Index="DiffExp_"${COMPAR_NUM}"c"
        COL_SUFFIX="Intronic_Only"
    fi
    
    # dir_name and job_name are required in the next steps
    dir_name=$(basename $(pwd))
    step_num=$(echo ${dir_name} | cut -d'_' -f 1)
    job_name="Step_${step_num}"
    
    #Now use arguments in the PBS script call:
    (set -x; qsub -N "${job_name}_${DiffExp_Index}" -P "${PROJECT}" -l h_rt="${TIME_LIMIT}" DiffExp.qsub ${SCRIPT_DIR} ${DATASET_DIR} ${DATASET_LABEL} ${GTF_FILES_DIR} ${ANNOTATION_FILE}  ${CONDITION_1_NAME} ${CONDITION_2_NAME} ${Lengths_DIR} ${GENE_LENGTHS_FILE} ${COUNT_DIR} ${OUTPUT_PREFIX} ${DiffExp_Index} ${COL_SUFFIX} ${COUNT_PROGRAM})

    #End loop over GTF files:
done
echo "End of qsub commands"
echo "-----------------------"
##################################################################################
