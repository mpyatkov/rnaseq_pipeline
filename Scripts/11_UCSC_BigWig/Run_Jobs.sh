#!/bin/bash -l
##################################################################################
#Andy Rampersaud, 08.07.17
#This script would be used to run UCSC_BigWig.pbs in parallel for different FASTQ files
#Way to run script:
#Usage: ./Run_Jobs.sh
#Example: 
#./Run_Jobs.sh
##################################################################################

set -o errexit
set -o pipefail
set -o nounset

#Remove *.o files from previous jobs
rm -rf ./logs && mkdir -p ./logs

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"


set +eu
module load anaconda2
source activate RNAseq
set -eu

# calculate STRAND_RULE from STRANDEDNESS
if [ ${STRANDEDNESS} -eq 0 ]
then
STRAND_RULE="none"
elif [ ${STRANDEDNESS} -eq 2 ]
then
STRAND_RULE="1++,1--,2+-,2-+"
elif [ ${STRANDEDNESS} -eq 1 ]
then
STRAND_RULE="1+-,1-+,2++,2--"
fi

SCRIPT_DIR="$(pwd)"

#Check that each variable prints a value to the terminal:
echo "-----------------------"
echo "Start of variable list:"
echo "-----------------------"
echo "DATASET_DIR:"
echo ${DATASET_DIR}
echo "STRAND_RULE:"
echo ${STRAND_RULE}
echo "BU_User:"
echo ${BU_USER}
echo "VM_DIR_UCSC:"
echo ${VM_DIR_UCSC}
echo "SCRIPT_DIR:"
echo ${SCRIPT_DIR}
echo "TIME_LIMIT:"
echo ${TIME_LIMIT}
echo "-----------------------"
echo "End of variable list"
echo "-----------------------"

dir_name=$(basename $(pwd))
step_num=$(echo ${dir_name} | cut -d'_' -f 1)
job_name="Step_${step_num}"

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

# loop over all samples
for ((i=0;i< ${#samples[@]} ;i+=3));
do
    Sample_DIR=${samples[i]}
    Sample_ID=${samples[i+1]}
    Description=${samples[i+2]}
    (set -x; qsub -N "${job_name}_${Sample_ID}" -P "${PROJECT}" -l h_rt="${TIME_LIMIT}" UCSC_BigWig.qsub ${Sample_ID} ${DATASET_DIR} ${STRAND_RULE} ${BU_USER} ${VM_DIR_UCSC} ${SCRIPT_DIR})

done 
#Remove the temp file:

echo "End of qsub commands"
echo "-----------------------"
##################################################################################
