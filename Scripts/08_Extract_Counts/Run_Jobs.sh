#!/bin/bash -l

set -o errexit
set -o pipefail
set -o nounset

rm -rf ./logs && mkdir -p ./logs

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

# load anaconda module for case when we need independent run
# set +eu
# module load anaconda2
# source activate RNAseq
# set -eu

# dir_name and job_name are required in the next steps
dir_name=$(basename $(pwd))
step_num=$(echo ${dir_name} | cut -d'_' -f 1)
job_name="${DATASET_LABEL}_${step_num}"

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))
gtflist=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --gtf_annotation_and_counter))

# loop over all samples
for ((i=0;i< ${#samples[@]} ;i+=3))
do
    sample_dir=${samples[i]}
    sample_id=${samples[i+1]}
    description=${samples[i+2]}
    
    for ((j=0; j< ${#gtflist[@]} ;j+=2))
    do
        gtf_file=${gtflist[j]}
        counter=${gtflist[j+1]}
        (set -x; qsub -N "${job_name}_${sample_id}" -P "${PROJECT}" -l h_rt="04:00:00" Extract_Counts.qsub ${sample_id} ${gtf_file} ${counter})
    done
    
done

echo "End of 8 step commands"
echo "-----------------------"

