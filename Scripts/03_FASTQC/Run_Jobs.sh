#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

##################################################################################
#Andy Rampersaud, 02.23.16
#This script would be used to run FASTQC.pbs in parallel for different FASTQ files
#Way to run script:
#Usage: ./Run_Jobs.sh
#Example: 
# ./Run_Jobs.sh has one parameter - FULL_RECALC, if specify nothing by
# default will be used FULL_RECALC=0 (do not recalc)
# ./Run_Jobs.sh 1 -- recalculates all
##################################################################################

# FULL_RECALC=0 by default if nothing provided
FULL_RECALC=${1:-0}

# Skip this step if recalculation flag set to 0
if [ ${FULL_RECALC} -eq 0 ]; then
    echo "Recalculation is not required. FULL_RECALC set to 0."
    exit 0
fi

# remove *.o files from previous jobs
rm -rf ./logs
mkdir -p logs

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

# dir_name and job_name are required in the next steps
dir_name=$(basename $(pwd))
step_num=$(echo ${dir_name} | cut -d'_' -f 1)
job_name="Step_${step_num}"

# output dir to store text files:
output_dir="$(pwd)/Job_Summary"
#output_dir="${DATASET_DIR}/Scripts/${dir_name}/Job_Summary"
rm -rf "${output_dir}" && mkdir -p "${output_dir}"

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

# rename directories from sample_dir to sample_id
for ((i=0;i< ${#samples[@]} ;i+=2));
do
    sample_id=${samples[i]}
    # description=${samples[i+1]}

    (set -x; qsub -N "${job_name}_${sample_id}" -P "${PROJECT}" -l h_rt="${TIME_LIMIT}" FASTQC.qsub "${sample_id}" "${output_dir}")
done

echo "End of FASTQC commands"
echo "-----------------------"
##################################################################################
