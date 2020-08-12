#!/bin/bash

# set -o errexit
set -o pipefail
set -o nounset
# set -x
##################################################################################
#Andy Rampersaud, 02.22.16
#Adapted from tophat/paired_end_mapping scripts by Tisha Melia
#This script would be used to run TopHat_Paired_End.pbs in parallel for different FASTQ files
#Way to run script:
#Usage: ./Run_Jobs.sh
#Example: 
#./Run_Jobs.sh
##################################################################################

#Remove *.o files from previous jobs
rm -rf *.e* *.o* *.po* *.pe*

# we activated this module for case when we need run this script independently
module load anaconda2
source activate RNAseq

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

# calculate distance_bt_read_pair from BIOANALYZER_LEN, ADAPTOR_LEN, READ_LEN
fragment_len=$(echo "scale=4;${BIOANALYZER_LEN}-(2*${ADAPTOR_LEN})" | bc)
distance_bt_read_pair=$(echo "scale=4;${fragment_len}-(2*${READ_LEN})" | bc)

# calculate STRANDEDNESS for TOPHAT
if [ ${STRANDEDNESS} -eq 0 ]
then
    tophat_strandedness="fr-unstranded"
elif [ ${STRANDEDNESS} -eq 1 ]
then
    tophat_strandedness="fr-firststrand"
elif [ ${STRANDEDNESS} -eq 2 ]
then
    tophat_strandedness="fr-secondstrand"
fi

# dir_name and job_name are required in the next steps
dir_name=$(basename $(pwd))
step_num=$(echo ${dir_name} | cut -d'_' -f 1)
job_name="Step_${step_num}"

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

# rename directories from sample_dir to sample_id
for ((i=0;i< ${#samples[@]} ;i+=3));
do
    # sample_dir=${samples[i]}
    sample_id=${samples[i+1]}
    # description=${samples[i+2]}
    
    (set -x; qsub -N "${job_name}_${sample_id}" -P "${PROJECT}" -l h_rt="${TIME_LIMIT}" TopHat_Paired_End.qsub "${sample_id}" "${tophat_strandedness}" "${distance_bt_read_pair}")
done

echo "End of TOPHAT commands"
echo "-----------------------"
##################################################################################
