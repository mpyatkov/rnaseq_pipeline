#!/bin/bash -l

set -o errexit
set -o pipefail
set -o nounset

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

#Remove *.o files from previous jobs
rm -rf ./logs && mkdir -p ./logs

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

# dir_name and job_name are required in the next steps
dir_name=$(basename $(pwd))
step_num=$(echo ${dir_name} | cut -d'_' -f 1)
job_name="Step_${step_num}"

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

# loop over all samples
for ((i=0;i< ${#samples[@]} ;i+=3));
do
    # sample_dir=${samples[i]}
    sample_id=${samples[i+1]}
    # description=${samples[i+2]}

    (set -x; qsub -N "${job_name}_${sample_id}" -P "${PROJECT}" -l h_rt="${TIME_LIMIT}" CollectMetrics.qsub "${sample_id}")

done

echo "End of qsub commands"

##################################################################################
