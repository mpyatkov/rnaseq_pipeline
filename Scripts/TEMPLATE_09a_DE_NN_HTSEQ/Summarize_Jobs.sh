#!/bin/bash
##################################################################################
#Andy Rampersaud, 03.12.16
#This script would be used to summarize DiffExp_* jobs
#Way to run script:
#Usage: 
#./DiffExp_Summary.sh
##################################################################################

set -o errexit
set -o pipefail
set -o nounset

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

SCRIPT_DIR="$(pwd)"

echo 'Loading required modules...'

module load R/3.6.0

#Source job-specific variables:
source setup_DiffExp.sh

#Need to get M_Num_Cond1_List:
#Text file has a header line to ignore:
tail -n +2 Condition_1.txt > Condition_1.temp
echo "------------------------------------------"

#Use a while loop to run jobs
while IFS=$'\t' read -r -a myArray
do
    #---------------------------
    ##Check that text file is read in properly:
    #echo 'Sample_DIR:'
    Sample_DIR=${myArray[0]}
    #echo $Sample_DIR
    #echo 'Sample_ID:'
    Sample_ID=${myArray[1]}
    #echo $Sample_ID
    #echo 'Description:'
    Description=${myArray[2]}
    #echo $Description

    #Split by underscore and get the 2nd part (M number)
    M_Num=$(echo $Sample_ID | cut -d'_' -f 2)
    #echo "M_Num:"
    #echo ${M_Num}
    #Append variable string to itself:
    M_Num_Cond1_List+=${M_Num}

done < Condition_1.temp

echo "M_Num_Cond1_List:"
echo ${M_Num_Cond1_List}

#Remove the temp file:
rm Condition_1.temp
echo "------------------------------------------"

#Need to get M_Num_Cond2_List:
#Text file has a header line to ignore:
tail -n +2 Condition_2.txt > Condition_2.temp
echo "------------------------------------------"
#Use a while loop to run jobs
while IFS=$'\t' read -r -a myArray
do

    ##Check that text file is read in properly:
    #echo 'Sample_DIR:'
    Sample_DIR=${myArray[0]}
    #echo $Sample_DIR
    #echo 'Sample_ID:'
    Sample_ID=${myArray[1]}
    #echo $Sample_ID
    #echo 'Description:'
    Description=${myArray[2]}
    #echo $Description

    #Split by underscore and get the 2nd part (M number)
    M_Num=$(echo $Sample_ID | cut -d'_' -f 2)
    #echo "M_Num:"
    #echo ${M_Num}
    #Append variable string to itself:
    M_Num_Cond2_List+=${M_Num}

done < Condition_2.temp

echo "M_Num_Cond2_List:"
echo ${M_Num_Cond2_List}

#Remove the temp file:
rm Condition_2.temp
echo "------------------------------------------"
SEGEX_Platform="DiffExp_v2"
echo "SEGEX_Platform:"
echo ${SEGEX_Platform}
echo "------------------------------------------"
#Check that each variable prints a value to the terminal:
echo "-----------------------"
echo "Start of variable list:"
echo "-----------------------"
echo "DATASET_LABEL:"
echo ${DATASET_LABEL}
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
echo "------------------------------------------"
Comparison_Info=${CONDITION_2_NAME}'.'${DATASET_LABEL}'.'${M_Num_Cond2_List}'.'${CONDITION_1_NAME}'.'${DATASET_LABEL}'.'${M_Num_Cond1_List}'.'${SEGEX_Platform}
echo "Comparison_Info:"
echo ${Comparison_Info}
echo "------------------------------------------"

##################################################################################
INPUT_DIR=$(pwd)
OUTPUT_DIR=${INPUT_DIR}/Summary_Differential_Expression

if [ ! -d ${OUTPUT_DIR} ]; then
    mkdir -p ${OUTPUT_DIR}
else

    #Remove dir:
    rm -r ${OUTPUT_DIR}
    #Make new dir:
    mkdir -p ${OUTPUT_DIR}
fi

SEGEX_Upload_DIR=${OUTPUT_DIR}/SEGEX_Upload_Files

if [ ! -d ${SEGEX_Upload_DIR} ]
then 
    mkdir -p ${SEGEX_Upload_DIR}
else 
    rm -r ${SEGEX_Upload_DIR}/*
fi

Combined_DIR=${SEGEX_Upload_DIR}/Combined

if [ ! -d ${Combined_DIR} ]
then 
    mkdir -p ${Combined_DIR}
else 
    rm ${Combined_DIR}/*.txt
fi

Individual_DIR=${SEGEX_Upload_DIR}/Individual

if [ ! -d ${Individual_DIR} ]
then 
    mkdir -p ${Individual_DIR}
else 
    rm ${Individual_DIR}/*.txt
fi

DESEQ_OUTPUT_FILE=${Combined_DIR}/Combined_forSEGEXUpload_DESeq'.'${COUNT_PROGRAM}'.'${Comparison_Info}'.txt'
EDGER_OUTPUT_FILE=${Combined_DIR}/Combined_forSEGEXUpload_EdgeR'.'${COUNT_PROGRAM}'.'${Comparison_Info}'.txt'

if [ -f ${DESEQ_OUTPUT_FILE} ]
then 
    rm ${DESEQ_OUTPUT_FILE}
else
    touch ${DESEQ_OUTPUT_FILE}
fi

if [ -f ${EDGER_OUTPUT_FILE} ]
then 
    rm ${EDGER_OUTPUT_FILE}
else
    touch ${EDGER_OUTPUT_FILE}
fi
######################
PNG_DIR=${OUTPUT_DIR}/PNG_Files
######################
if [ ! -d ${PNG_DIR} ]
then 
    mkdir -p ${PNG_DIR}
else 
    rm -r ${PNG_DIR}/*
fi
######################
DE_Text_DIR=${OUTPUT_DIR}/DE_Text_Files
######################
if [ ! -d ${DE_Text_DIR} ]
then 
    mkdir -p ${DE_Text_DIR}
else 
    rm -r ${DE_Text_DIR}/*
fi
######################
#For each output folder:
#Copy the *_forSEGEXUpload*.txt file to the OUTPUT_DIR
#Then paste the *_forSEGEXUpload*.txt files to get the combined file
#---------------------------
Output_list=Output_DiffExp_*
for Output in $Output_list
do
    echo "-----------------------"
    echo ${Output}
    echo "-----------------------"
    #---------------------------
    echo 'Copying *_forSEGEXUpload*.txt'
    cd ${Output}
    cp *_forSEGEXUpload*.txt ${OUTPUT_DIR}
    #---------------------------
    cd ..
done
#---------------------------
cd ${OUTPUT_DIR}
#---------------------------------------------------------------------------------
#Need to paste files together in the following order:
#<GeneBody> <Intronic_Only> <Exonic_Only>
GeneBody_File=$(ls *'DESeq_'${COUNT_PROGRAM}'.txt' | grep 'GeneBody')
Intronic_Only_File=$(ls *'DESeq_'${COUNT_PROGRAM}'.txt' | grep 'Intronic_Only')
Exonic_Only_File=$(ls *'DESeq_'${COUNT_PROGRAM}'.txt' | grep 'Exonic_Only')
echo 'Running paste command'
paste ${GeneBody_File} ${Intronic_Only_File} ${Exonic_Only_File} > ${DESEQ_OUTPUT_FILE}
#---------------------------------------------------------------------------------
#Need to paste files together in the following order:
#<GeneBody> <Intronic_Only> <Exonic_Only>
GeneBody_File=$(ls *'EdgeR_'${COUNT_PROGRAM}'.txt' | grep 'GeneBody')
Intronic_Only_File=$(ls *'EdgeR_'${COUNT_PROGRAM}'.txt' | grep 'Intronic_Only')
Exonic_Only_File=$(ls *'EdgeR_'${COUNT_PROGRAM}'.txt' | grep 'Exonic_Only')
echo 'Running paste command'
paste ${GeneBody_File} ${Intronic_Only_File} ${Exonic_Only_File} > ${EDGER_OUTPUT_FILE}
#---------------------------------------------------------------------------------
echo "Move files to Individual_DIR"
mv *_forSEGEXUpload_'DESeq_'${COUNT_PROGRAM}'.txt' ${Individual_DIR}
mv *_forSEGEXUpload_'TPM_DESeq_'${COUNT_PROGRAM}'.txt' ${Individual_DIR}
mv *_forSEGEXUpload_'EdgeR_'${COUNT_PROGRAM}'.txt' ${Individual_DIR}
mv *_forSEGEXUpload_'TPM_EdgeR_'${COUNT_PROGRAM}'.txt' ${Individual_DIR}
#---------------------------------------------------------------------------------
OUTPUT_TABLE=${OUTPUT_DIR}/DE_Gene_Counts'.'${COUNT_PROGRAM}'.'${Comparison_Info}'.txt'
######################
if [ -f ${OUTPUT_TABLE} ]
then 
    rm ${OUTPUT_TABLE}
else
    touch ${OUTPUT_TABLE}
fi
######################
OUTPUT_TABLE_2=${OUTPUT_DIR}/DE_Gene_Venn_Tables'.'${COUNT_PROGRAM}'.'${Comparison_Info}'.txt'
######################
if [ -f ${OUTPUT_TABLE_2} ]
then 
    rm ${OUTPUT_TABLE_2}
else
    touch ${OUTPUT_TABLE_2}
fi
######################
#---------------------------------------------------------------------------------
#Print header to output file:
echo 'DiffExp_Index' $'\t''Output_File' $'\t''Up_Genes_Count' $'\t''Down_Genes_Count' >> ${OUTPUT_TABLE}
#---------------------------------------------------------------------------------
cd ${INPUT_DIR}
#---------------------------
#---------------------------------------------------------------------------------
Output_list=Output_DiffExp_*
for Output in $Output_list
do
    echo "-----------------------"
    echo ${Output}
    echo "-----------------------"
    DiffExp_Index=${Output}
    cd ${Output}
    DESeq_Output_File=$(ls *DiffExp_v2_*DESeq*.txt)
    echo 'Count for Up_Genes_DESeq_*.txt'
    Up_Genes_Count=$(wc -l Up_Genes_DESeq_*.txt | awk '{print $1-1}')
    echo 'Count for Down_Genes_DESeq_*.txt'
    Down_Genes_Count=$(wc -l Down_Genes_DESeq_*.txt | awk '{print $1-1}')
    #---------------------------
    echo 'Print to output file'
    echo ${DiffExp_Index} $'\t'${DESeq_Output_File} $'\t'${Up_Genes_Count} $'\t'${Down_Genes_Count}  >> ${OUTPUT_TABLE}
    #---------------------------
    EdgeR_Output_File=$(ls *DiffExp_v2_*EdgeR*.txt)
    echo 'Count for Up_Genes_EdgeR_*.txt'
    Up_Genes_Count=$(wc -l Up_Genes_EdgeR_*.txt | awk '{print $1-1}')
    echo 'Count for Down_Genes_EdgeR_*.txt'
    Down_Genes_Count=$(wc -l Down_Genes_EdgeR_*.txt | awk '{print $1-1}')
    #---------------------------
    echo 'Print to OUTPUT_TABLE'
    echo ${DiffExp_Index} $'\t'${EdgeR_Output_File} $'\t'${Up_Genes_Count} $'\t'${Down_Genes_Count}  >> ${OUTPUT_TABLE}
    #---------------------------
    echo 'Print to OUTPUT_TABLE_2'
    cat *_Venn_Tables*.txt >> ${OUTPUT_TABLE_2} 
    #---------------------------
    echo 'Copy PNG files to PNG_DIR'
    cp *.png ${PNG_DIR}
    #---------------------------
    echo 'Copy Down_Genes*.txt and Up_Genes*.txt files to DE_Text_DIR'
    cp Down_Genes*.txt ${DE_Text_DIR}
    cp Up_Genes*.txt ${DE_Text_DIR}
    #---------------------------
    cd ..
done
#---------------------------------------------------------------------------------
OUTPUT_FILE_PDF=${OUTPUT_DIR}/DE_Genes_Venn_R_Package'.'${COUNT_PROGRAM}'.'${Comparison_Info}'.pdf'
######################
if [ -f ${OUTPUT_FILE_PDF} ]
then 
    rm ${OUTPUT_FILE_PDF}
else
    touch ${OUTPUT_FILE_PDF}
fi

echo "-----------------------"
echo 'Create montage of PNG files'
cd ${PNG_DIR}
montage -geometry 500x500 *.png ${OUTPUT_FILE_PDF}
echo "-----------------------"
cd ..
#Remove ${PNG_DIR}:
rm -r ${PNG_DIR}

cd ${DE_Text_DIR}
OUTPUT_FILE_PDF_2=${OUTPUT_DIR}/DE_Genes_Venn_Count_Method'.'${COUNT_PROGRAM}'.'${Comparison_Info}'.pdf'

######################
if [ -f ${OUTPUT_FILE_PDF_2} ]
then 
    rm ${OUTPUT_FILE_PDF_2}
else
    touch ${OUTPUT_FILE_PDF_2}
fi

######################
echo 'Copy Venn_Count_Method.R script to DE_Text_DIR'
cp ${SCRIPT_DIR}/Scripts/Venn_Count_Method.R . 
echo '#--------------------------------------------------------------------------'
echo 'Run Venn_Count_Method.R commands'
echo '#--------------------------------------------------------------------------'

#Since renaming the files to include ${Comparison_Info} in the qsub
#I need a different way to get the file names:
#DESeq down files:
DESeq_Down_File1=$(ls Down_Genes_DESeq_Exonic_Only*)
DESeq_Down_File2=$(ls Down_Genes_DESeq_GeneBody*)
DESeq_Down_File3=$(ls Down_Genes_DESeq_Intronic_Only*)
echo 'Input files for Venn_Count_Method.R:'
echo 'DESeq_Down_File1:'
echo ${DESeq_Down_File1}
echo 'DESeq_Down_File2:'
echo ${DESeq_Down_File2}
echo 'DESeq_Down_File3:'
echo ${DESeq_Down_File3}
echo '#--------------------------------------------------------------------------'
Rscript Venn_Count_Method.R ${DESeq_Down_File1} ${DESeq_Down_File2} ${DESeq_Down_File3} 
echo '#--------------------------------------------------------------------------'

#DESeq up files:
DESeq_Up_File1=$(ls Up_Genes_DESeq_Exonic_Only*)
DESeq_Up_File2=$(ls Up_Genes_DESeq_GeneBody*)
DESeq_Up_File3=$(ls Up_Genes_DESeq_Intronic_Only*)
echo 'Input files for Venn_Count_Method.R:'
echo 'DESeq_Up_File1:'
echo ${DESeq_Up_File1}
echo 'DESeq_Up_File2:'
echo ${DESeq_Up_File2}
echo 'DESeq_Up_File3:'
echo ${DESeq_Up_File3}
echo '#--------------------------------------------------------------------------'
Rscript Venn_Count_Method.R ${DESeq_Up_File1} ${DESeq_Up_File2} ${DESeq_Up_File3} 
echo '#--------------------------------------------------------------------------'

#EdgeR down files:
EdgeR_Down_File1=$(ls Down_Genes_EdgeR_Exonic_Only*)
EdgeR_Down_File2=$(ls Down_Genes_EdgeR_GeneBody*)
EdgeR_Down_File3=$(ls Down_Genes_EdgeR_Intronic_Only*)
echo 'Input files for Venn_Count_Method.R:'
echo 'EdgeR_Down_File1:'
echo ${EdgeR_Down_File1}
echo 'EdgeR_Down_File2:'
echo ${EdgeR_Down_File2}
echo 'EdgeR_Down_File3:'
echo ${EdgeR_Down_File3}
echo '#--------------------------------------------------------------------------'
Rscript Venn_Count_Method.R ${EdgeR_Down_File1} ${EdgeR_Down_File2} ${EdgeR_Down_File3} 
echo '#--------------------------------------------------------------------------'

#EdgeR up files:
EdgeR_Up_File1=$(ls Up_Genes_EdgeR_Exonic_Only*)
EdgeR_Up_File2=$(ls Up_Genes_EdgeR_GeneBody*)
EdgeR_Up_File3=$(ls Up_Genes_EdgeR_Intronic_Only*)
echo 'Input files for Venn_Count_Method.R:'
echo 'EdgeR_Up_File1:'
echo ${EdgeR_Up_File1}
echo 'EdgeR_Up_File2:'
echo ${EdgeR_Up_File2}
echo 'EdgeR_Up_File3:'
echo ${EdgeR_Up_File3}
echo '#--------------------------------------------------------------------------'
Rscript Venn_Count_Method.R ${EdgeR_Up_File1} ${EdgeR_Up_File2} ${EdgeR_Up_File3} 
echo '#--------------------------------------------------------------------------'
echo 'Create montage of PNG files'
montage -geometry 500x500 *.png ${OUTPUT_FILE_PDF_2}
echo '#--------------------------------------------------------------------------'
cd ..

#Remove ${DE_Text_DIR}:
rm -r ${DE_Text_DIR}
echo '#--------------------------------------------------------------------------'
echo 'Check out '${OUTPUT_DIR}
#echo 'Check out '${OUTPUT_TABLE}
#echo 'Check out '${OUTPUT_TABLE_2}
#echo 'Check out '${OUTPUT_FILE_PDF}
#echo 'Check out '${OUTPUT_FILE_PDF_2}
echo '#--------------------------------------------------------------------------'
##################################################################################
