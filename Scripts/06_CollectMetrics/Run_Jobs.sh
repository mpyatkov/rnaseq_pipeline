#!/bin/bash -l

set -o errexit
set -o pipefail
set -o nounset

# ./Run_Jobs.sh 0 -- allows to reuse previously calculated results
# ./Run_Jobs.sh 1 -- recalculate all

# FULL_RECALC equal to number after ':-' if parameter was not provided
FULL_RECALC=${1:-0}

# Skip this step if recalculation flag set to 0
if [ ${FULL_RECALC} -eq 0 ]; then
    echo "Recalculation is not required. FULL_RECALC set to 0."
    exit 0
fi

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

#Remove *.o files from previous jobs
rm -rf ./logs && mkdir -p ./logs


# dir_name and job_name are required in the next steps
dir_name=$(basename $(pwd))
step_num=$(echo ${dir_name} | cut -d'_' -f 1)
job_name="Step_${step_num}_${DATASET_LABEL}"

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

# loop over all samples
for ((i=0;i< ${#samples[@]} ;i+=2));
do
    sample_id=${samples[i]}
    # description=${samples[i+1]}

    (set -x; qsub -N "${job_name}_${sample_id}" -P "${PROJECT}" -l h_rt="${TIME_LIMIT}" CollectMetrics.qsub "${sample_id}")

done

echo "End of qsub commands"

##################################################################################
