#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

##################################################################################
#Andy Rampersaud, 02.22.16
#Adapted from tophat/paired_end_mapping scripts by Tisha Melia
#This script would be used to run TopHat_Paired_End.pbs in parallel for different FASTQ files
#Way to run script:
#Usage: ./Run_Jobs.sh
#Example: 
# ./Run_Jobs.sh 0 -- allows to reuse previously calculated results
# ./Run_Jobs.sh 1 -- recalculate all

# FULL_RECALC equal to number after ':-' if parameter was not provided
FULL_RECALC=${1:-0}

#Remove *.o files from previous jobs
rm -rf *.po* *.pe*
rm -rf ./logs
mkdir -p ./logs

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

# calculate distance_bt_read_pair from BIOANALYZER_LEN, ADAPTOR_LEN, READ_LEN
fragment_len=$(echo "scale=4;${BIOANALYZER_LEN}-(2*${ADAPTOR_LEN})" | bc)
distance_bt_read_pair=$(echo "scale=4;${fragment_len}-(2*${READ_LEN})" | bc)

# dir_name and job_name are required in the next steps
dir_name=$(basename $(pwd))
step_num=$(echo ${dir_name} | cut -d'_' -f 1)
job_name="Step_${step_num}_${DATASET_LABEL}"

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

# rename directories from sample_dir to sample_id
for ((i=0;i< ${#samples[@]} ;i+=2));
do
    sample_id=${samples[i]}
    # description=${samples[i+1]}

    ## skip this step if BAM file exists
    result_path="${DATASET_DIR}/${sample_id}/aligner/${sample_id}_primary_unique.bam"
    if [ -f ${result_path} -a ${FULL_RECALC} -eq 0 ]; then
	echo "Skip recalculation for ${sample_id} sample. Results already obtained and FULL_RECALC flag set to 0."
	continue
    fi

    
    if [ ${DEFAULT_ALIGNER} -eq 0 ]; then
	echo "TopHat runned..."
	(set -x; qsub -N "${job_name}_${sample_id}" -P "${PROJECT}" -l h_rt="${TIME_LIMIT}" TopHat_Paired_End.qsub "${sample_id}" "${distance_bt_read_pair}")
    else
	echo "STAR runned..."
	(set -x; qsub -N "${job_name}_${sample_id}" -P "${PROJECT}" -l h_rt="${TIME_LIMIT}" STAR_Paired_End.qsub "${sample_id}")
    fi
done

echo "End of TOPHAT commands"
echo "-----------------------"
##################################################################################
