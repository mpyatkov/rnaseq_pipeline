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

rm -rf SEGEX_Upload_Files
# if [ -d Summary_Differential_Expression -o -d SEGEX_Upload_Files ]
# then
#     echo "No summary required. Output folders exist"
#     exit 0
# fi

echo
echo 'Loading required modules...'
echo

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

set +eu
module load miniconda
conda activate --stack ${CONDA_DIR}/rlang361
set -eu

SCRIPT_DIR=$(pwd)

# list of samples for each condition (ex. G186_M1M2M3_G184_M1M2M3)
samples_line(){
    local fname=$1
    a=$(tail -n+2 $fname | cut -f1 | cut -d "_" -f 1 | sort | uniq)

    res=""
    for i in $a
    do
	group="${i}_$(grep ${i} $fname | cut -f 1 | cut -d "_" -f 2 | paste -s -d '')_"
	res+=$group
	# res+="_"
    done

    echo $res
}

M_Num_Cond1_List=$(samples_line Condition_1.txt)
M_Num_Cond1_List=${M_Num_Cond1_List%\_}
M_Num_Cond2_List=$(samples_line Condition_2.txt)
M_Num_Cond2_List=${M_Num_Cond2_List%\_}

#Need to get M_Num_Cond1_List:
# M_Num_Cond1_List=$(tail -n+2 Condition_1.txt | cut -f 2| cut -d "_" -f2 | paste -s -d "")
# M_Num_Cond2_List=$(tail -n+2 Condition_2.txt | cut -f 2| cut -d "_" -f2 | paste -s -d "")

CONDITION_1_NAME=TEMPLATE
CONDITION_2_NAME=TEMPLATE
DE_INDEX=TEMPLATE
COMPAR_NUM=TEMPLATE
COUNT_PROGRAM=TEMPLATE

# TODO: refactor this, get all information from Pipeline_Setup.py
# CONDITION_1_NAME=$(tail -n+2 Condition_1.txt | cut -f 3 | sort | uniq | head -1)
# CONDITION_2_NAME=$(tail -n+2 Condition_2.txt | cut -f 3 | sort | uniq | head -1)
# DE_INDEX=$(pwd | xargs -n1 basename | grep -Po "\K([0-9a-zA-Z]*)(?=_)" | head -1)
# COUNT_PROGRAM="$(../00_Setup_Pipeline/01_Pipeline_Setup.py --counter_by_DE_INDEX ${DE_INDEX})"

SEGEX_Platform="DiffExp_v2"

#Check that each variable prints a value to the terminal:
Comparison_Info=${CONDITION_2_NAME}'.'${M_Num_Cond2_List}'.'${CONDITION_1_NAME}'.'${M_Num_Cond1_List}'.'${SEGEX_Platform}

echo ${Comparison_Info}

INPUT_DIR=$(pwd)
# PREFIXNUM=$(echo ${INPUT_DIR} | grep -Po "_\K([0-9][0-9]?)(?=_)")

OUTPUT_DIR=${INPUT_DIR}/Summary_Differential_Expression
rm -rf ${OUTPUT_DIR} && mkdir -p ${OUTPUT_DIR}

SEGEX_Upload_DIR=${OUTPUT_DIR}/SEGEX_Upload_Files
rm -rf ${SEGEX_Upload_DIR} && mkdir -p ${SEGEX_Upload_DIR}

Combined_DIR=${SEGEX_Upload_DIR}/Combined
rm -rf ${Combined_DIR} && mkdir -p ${Combined_DIR}

Individual_DIR=${SEGEX_Upload_DIR}/Individual
rm -rf ${Individual_DIR} && mkdir -p ${Individual_DIR}

DESEQ_OUTPUT_FILE=${Combined_DIR}/${COMPAR_NUM}_Combined_forSEGEXUpload_DESeq'.'${COUNT_PROGRAM}'.'${Comparison_Info}'.txt'
EDGER_OUTPUT_FILE=${Combined_DIR}/${COMPAR_NUM}_Combined_forSEGEXUpload_EdgeR'.'${COUNT_PROGRAM}'.'${Comparison_Info}'.txt'

rm -rf ${DESEQ_OUTPUT_FILE} && touch ${DESEQ_OUTPUT_FILE}
rm -rf ${EDGER_OUTPUT_FILE} && touch ${EDGER_OUTPUT_FILE}

# PNG_DIR=${OUTPUT_DIR}/PNG_Files
# rm -rf ${PNG_DIR} && mkdir -p ${PNG_DIR}

DE_Text_DIR=${OUTPUT_DIR}/DE_Text_Files
rm -rf ${DE_Text_DIR} && mkdir -p ${DE_Text_DIR}

# For each output folder:
# Copy the *_forSEGEXUpload*.txt file to the OUTPUT_DIR
# Then paste the *_forSEGEXUpload*.txt files to get the combined file

segex_files=$(find ./Output* -name "*_forSEGEXUpload*.txt")
for sf in ${segex_files}
do
    fname=$(basename $sf)
    cp $sf ${OUTPUT_DIR}/${COMPAR_NUM}_${fname}
done

pushd ${OUTPUT_DIR}

#---------------------------------------------------------------------------------
#Need to paste files together in the following order:
#<GeneBody> <Intronic_Only> <Exonic_Only>
set +eu
ExonCollapsed_File=$(find . -name "*DESeq_${COUNT_PROGRAM}.txt" | grep "Exon_\|ExonCollapsed\|Exon_Collapsed" | grep -vi "only")
GeneBody_File=$(find . -name  "*DESeq_${COUNT_PROGRAM}.txt" | grep "GeneBody")
Intronic_Only_File=$(find . -name  "*DESeq_${COUNT_PROGRAM}.txt" | grep 'Intron\|Intronic_Only')
Exonic_Only_File=$(find . -name  "*DESeq_${COUNT_PROGRAM}.txt" | grep 'Exon_Only\|Exonic_Only')
echo "${ExonCollapsed_File} ${GeneBody_File} ${Intronic_Only_File} ${Exonic_Only_File} > ${DESEQ_OUTPUT_FILE}"
echo 'Running paste command'
(set +x; paste ${ExonCollapsed_File} ${GeneBody_File} ${Intronic_Only_File} ${Exonic_Only_File} > ${DESEQ_OUTPUT_FILE})

#---------------------------------------------------------------------------------
#Need to paste files together in the following order:
#<GeneBody> <Intronic_Only> <Exonic_Only>
ExonCollapsed_File=$(find . -name  "*EdgeR_${COUNT_PROGRAM}.txt" | grep "Exon_f\|ExonCollapsed\|Exon_Collapsed")
GeneBody_File=$(find . -name  "*EdgeR_${COUNT_PROGRAM}.txt" | grep 'GeneBody')
Intronic_Only_File=$(find . -name  "*EdgeR_${COUNT_PROGRAM}.txt" | grep 'Intron\|Intronic_Only')
Exonic_Only_File=$(find . -name  "*EdgeR_${COUNT_PROGRAM}.txt" | grep 'Exonic_Only')

echo 'Running paste command'
(set +x; paste ${ExonCollapsed_File} ${GeneBody_File} ${Intronic_Only_File} ${Exonic_Only_File} > ${EDGER_OUTPUT_FILE})
#---------------------------------------------------------------------------------
set -eu

IFPKMDIR=${Individual_DIR}/Segex_FPKM
mkdir ${IFPKMDIR}
mv *_forSEGEXUpload_'DESeq_'${COUNT_PROGRAM}'.txt' ${IFPKMDIR}
mv *_forSEGEXUpload_'EdgeR_'${COUNT_PROGRAM}'.txt' ${IFPKMDIR}

ITPMDIR=${Individual_DIR}/Segex_TPM
mkdir ${ITPMDIR}
mv *_forSEGEXUpload_'TPM_DESeq_'${COUNT_PROGRAM}'.txt' ${ITPMDIR}
mv *_forSEGEXUpload_'TPM_EdgeR_'${COUNT_PROGRAM}'.txt' ${ITPMDIR}

mv ${SEGEX_Upload_DIR} ../

# move to main dir
popd

OUTPUT_TABLE=${OUTPUT_DIR}/DE_Gene_Counts'.'${COUNT_PROGRAM}'.'${Comparison_Info}'.txt'
rm -rf ${OUTPUT_TABLE} && touch ${OUTPUT_TABLE}

OUTPUT_TABLE_2=${OUTPUT_DIR}/DE_Gene_Venn_Tables'.'${COUNT_PROGRAM}'.'${Comparison_Info}'.txt'
rm -rf ${OUTPUT_TABLE_2} && touch ${OUTPUT_TABLE_2}

#Print header to output file:
echo 'DiffExp_Index' $'\t''Output_File' $'\t''Up_Genes_Count' $'\t''Down_Genes_Count' >> ${OUTPUT_TABLE}

cd ${INPUT_DIR}

Output_list=Output_DiffExp_*
for Output in $Output_list
do
    echo ${Output}
    DiffExp_Index=${Output}

    cd ${Output}
    DESeq_Output_File=$(ls *DiffExp_v2_*DESeq*.txt)

    echo 'Count for Up_Genes_DESeq_*.txt'
    Up_Genes_Count=$(wc -l Up_Genes_DESeq_*.txt | awk '{print $1-1}')

    echo 'Count for Down_Genes_DESeq_*.txt'
    Down_Genes_Count=$(wc -l Down_Genes_DESeq_*.txt | awk '{print $1-1}')

    echo 'Print to output file'
    echo ${DiffExp_Index} $'\t'${DESeq_Output_File} $'\t'${Up_Genes_Count} $'\t'${Down_Genes_Count}  >> ${OUTPUT_TABLE}
    EdgeR_Output_File=$(ls *DiffExp_v2_*EdgeR*.txt)

    echo 'Count for Up_Genes_EdgeR_*.txt'
    Up_Genes_Count=$(wc -l Up_Genes_EdgeR_*.txt | awk '{print $1-1}')

    echo 'Count for Down_Genes_EdgeR_*.txt'
    Down_Genes_Count=$(wc -l Down_Genes_EdgeR_*.txt | awk '{print $1-1}')

    echo 'Print to OUTPUT_TABLE'
    echo ${DiffExp_Index} $'\t'${EdgeR_Output_File} $'\t'${Up_Genes_Count} $'\t'${Down_Genes_Count}  >> ${OUTPUT_TABLE}

    set +eu
    echo 'Print to OUTPUT_TABLE_2'
    cat *_Venn_Tables*.txt >> ${OUTPUT_TABLE_2} 

    # echo 'Copy PNG files to PNG_DIR'
    # cp *.png ${PNG_DIR}
    set -eu

    echo 'Copy Down_Genes*.txt and Up_Genes*.txt files to DE_Text_DIR'
    cp Down_Genes*.txt ${DE_Text_DIR}
    cp Up_Genes*.txt ${DE_Text_DIR}
    cd ..
done

# remove not required anymore segex files
find . -name "*_forSEGEXUpload*.txt" | grep Output | xargs -n1 rm -rf 

# remove not required Venn diagrams
find . -name "*_Venn_*.txt" | grep Output | xargs -n1 rm -rf 

# TODO: remove outdated code
# this code is not required anymore, all venn diagrams functionaly
# moved to step 12
# >>>>>>>>>>

# OUTPUT_FILE_PDF=${OUTPUT_DIR}/DE_Genes_Venn_R_Package'.'${COUNT_PROGRAM}'.'${Comparison_Info}'.pdf'
# rm -rf ${OUTPUT_FILE_PDF} && touch ${OUTPUT_FILE_PDF}

# echo 'Create montage of PNG files'
# cd ${PNG_DIR}
# montage -geometry 500x500 *.png ${OUTPUT_FILE_PDF}
# cd ..
# rm -r ${PNG_DIR}

# cd ${DE_Text_DIR}
# OUTPUT_FILE_PDF_2=${OUTPUT_DIR}/DE_Genes_Venn_Count_Method'.'${COUNT_PROGRAM}'.'${Comparison_Info}'.pdf'
# # rm -rf ${OUTPUT_FILE_PDF_2} && touch ${OUTPUT_FILE_PDF_2}

# echo 'Copy Venn_Count_Method.R script to DE_Text_DIR'
# cp ${SCRIPT_DIR}/Scripts/Venn_Count_Method.R . 
# echo 'Run Venn_Count_Method.R commands'

#Since renaming the files to include ${Comparison_Info} in the qsub
#I need a different way to get the file names:
#DESeq down files:

# venn_count() {
#     # direction Up, Down
#     local d=$1
#     # package deseq, edger
#     local p=$2
    
#     # f1=$(ls ${d}_Genes_${p}_Exon* | grep "Exon_\|ExonCollapsed\|Exon_Collapsed")
#     # f2=$(ls ${d}_Genes_${p}_GeneBody*)
#     # f3=$(ls ${d}_Genes_${p}_Exonic_Only*)
#     # f4=$(ls ${d}_Genes_${p}_Intronic_Only*)
#     set +eu
#     f1=$(find . -name "${d}_Genes_${p}*" | grep "ExonCollapsed")
#     f2=$(find . -name "${d}_Genes_${p}*" | grep "GeneBody")
#     f3=$(find . -name "${d}_Genes_${p}*" | grep "ExonicOnly\|ExonOnly")
#     f4=$(find . -name "${d}_Genes_${p}*" | grep "IntronicOnly\|Intron")
#     set -eu
#     (set -x;Rscript Venn_Count_Method.R ${f1} ${f2} ${f3} ${f4})
#     echo "---------------->>>"
# }

# venn_count Down DESeq
# venn_count Up DESeq
# venn_count Down EdgeR
# venn_count Up EdgeR

# echo 'Create montage of PNG files'
# montage -geometry 500x500 *.png ${OUTPUT_FILE_PDF_2}
# cd ..
# <<<<<<<<<<

rm -r ${DE_Text_DIR}
echo '#--------------------------------------------------------------------------'
echo 'Check out '${OUTPUT_DIR}
echo '#--------------------------------------------------------------------------'
##################################################################################
