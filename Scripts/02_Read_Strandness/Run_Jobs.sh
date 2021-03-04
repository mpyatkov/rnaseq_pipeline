#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

##################################################################################
#Andy Rampersaud, 02.22.16
#This script would be used to run Read_Strandness.pbs in parallel for different FASTQ files
#Way to run script:
#Usage: ./Run_Jobs.sh
#Example: 
# ./Run_Jobs.sh has one parameter - FULL_RECALC, if specify nothing by
# default will be used FULL_RECALC=1 (full recalculation)
# ./Run_Jobs.sh 0 -- allows to reuse previously calculated resutls
##################################################################################

# FULL_RECALC=1 by default if nothing provided
FULL_RECALC=${1:-1}

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

#Remove *.o files from previous jobs
rm -rf ./logs Job_Summary
mkdir -p ./logs

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

# dir_name and job_name are required in the next steps
dir_name=$(basename $(pwd))
step_num=$(echo ${dir_name} | cut -d'_' -f 1)
job_name="Step_${step_num}"

# rename directories from sample_dir to sample_id
for ((i=0;i< ${#samples[@]} ;i+=3));
do
    # sample_dir=${samples[i]}
    sample_id=${samples[i+1]}
    # description=${samples[i+2]}

    # skip recalculation if FULL_RECALC=0
    result_path="${DATASET_DIR}/${sample_id}/Read_Strandness/${sample_id}_export.sh"
    if [ -f ${result_path} -a ${FULL_RECALC} -eq 0 ]; then
	echo "Skip recalculation for ${sample_id} sample. Results already obtained and FULL_RECALC flag set to 0."
	continue
    fi
    
    # (set -x; qsub -N "${job_name}_${sample_id}" -P "${PROJECT}" -l h_rt="${TIME_LIMIT}" Read_Strandness.qsub "${sample_id}")
    (set -x; qsub -N "${job_name}_${sample_id}" -P "${PROJECT}" -l h_rt="00:15:00" Read_Strandness.qsub "${sample_id}")
    
done

echo "End of 05_READ_STRANDEDNESS commands"
echo "-----------------------"
##################################################################################
