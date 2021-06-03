#!/bin/bash -l

set -o errexit
set -o pipefail
set -o nounset

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

#Remove *.o files from previous jobs
rm -rf ./logs Output_* Summary_Differential_Expression SEGEX_Upload_Files
mkdir -p ./logs

# load CONDITION_{1,2}_NAME
# load COMPAR_NUM <-- COMPARISON NUMBER?
# COUNT_PROGRAM <-- COUNTER
# source setup_DiffExp.sh

SCRIPT_DIR="$(pwd)"

# TODO: refactor this, get all information from Pipeline_Setup.py
# CONDITION_1_NAME=$(tail -n+2 Condition_1.txt | cut -f 3 | sort | uniq | head -1)
# CONDITION_2_NAME=$(tail -n+2 Condition_2.txt | cut -f 3 | sort | uniq | head -1)
# DE_INDEX=$(pwd | xargs -n1 basename | grep -Po "0\K([0-9][a-z])(?=_)" | head -1) # only 9a
# DE_INDEX=$(pwd | xargs -n1 basename | grep -Po "\K([0-9a-zA-Z]*)(?=_)" | head -1) # 09a current
# COMPAR_NUM=$(pwd | xargs -n1 basename | grep -Po "_\K([0-9][0-9]?)(?=_)" | head -1)


# get all parameters from config files
CONDITION_1_NAME=TEMPLATE
CONDITION_2_NAME=TEMPLATE
DE_INDEX=TEMPLATE
COMPAR_NUM=TEMPLATE

gtflist=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --gtf_by_DE_INDEX ${DE_INDEX}))

LENGTHS_DIR="${GTF_FILES_DIR}/lengths"

# dir_name and job_name are required in the next steps
dir_name=$(basename $(pwd))
step_num=$(echo ${dir_name} | cut -d'_' -f 1)
job_name="Step_${step_num}"

for ((i=0;i< ${#gtflist[@]} ;i+=2));
do
    gtf_file=${gtflist[i]}
    counter=${gtflist[i+1]}

    # we assume that output_dir looks like: GTFID_FTTYPE_GTF (ex: LncRNA_GeneBody_GTF)
    
    # export information about gtf file
    eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export_gtf_by_name_and_counter ${gtf_file} ${counter})"

    GENE_LENGTHS_FILE="${ANNOTATION_FILE%\.gtf}_lengths.txt"
    OUTPUT_PREFIX="DiffExp_v2_${OUTPUT_DIR%\_GTF}"
    COUNT_DIR="${OUTPUT_DIR}"
    DiffExp_Index="DiffExp_${COMPAR_NUM}${GTF_INDEX}"
    COL_SUFFIX=$(echo "${OUTPUT_DIR}" | sed 's/_[A-Z]*$//' | sed 's/^[a-zA-Z]*_//')
    
    if [[ ${CHANGE_ANNOTATION_FILE} != 0 ]]
    then
	echo "change annotation file: ${ANNOTATION_FILE} --> ${CHANGE_ANNOTATION_FILE}"
        ANNOTATION_FILE=${CHANGE_ANNOTATION_FILE}
    fi
        
    (set -x; qsub -N "${job_name}_${DiffExp_Index}" -P "${PROJECT}" -l h_rt="01:00:00" DiffExp.qsub ${SCRIPT_DIR} ${DATASET_DIR} ${DATASET_LABEL} ${GTF_FILES_DIR} ${ANNOTATION_FILE}  ${CONDITION_1_NAME} ${CONDITION_2_NAME} ${LENGTHS_DIR} ${GENE_LENGTHS_FILE} ${COUNT_DIR} ${OUTPUT_PREFIX} ${DiffExp_Index} ${COL_SUFFIX} ${COUNTER})
    
    
done
echo "End of qsub commands"
echo "-----------------------"

