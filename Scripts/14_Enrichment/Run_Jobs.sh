#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

#Clean directory
rm -rf ./logs && mkdir -p ./logs

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

rm -rf Job_Summary && mkdir -p Job_Summary

# dir_name and job_name are required in the next steps
dir_name=$(basename $(pwd))
step_num=$(echo ${dir_name} | cut -d'_' -f 1)
job_name="Step_${step_num}_${DATASET_LABEL}"

(set -x; qsub -N "${job_name}" -P "${PROJECT}" -l h_rt="02:00:00" enrichment.qsub)
