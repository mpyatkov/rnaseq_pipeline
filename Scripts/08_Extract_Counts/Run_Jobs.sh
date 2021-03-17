#!/bin/bash -l

set -o errexit
set -o pipefail
set -o nounset

# ./Run_Jobs.sh 0 -- allows to reuse previously calculated results
# ./Run_Jobs.sh 1 -- recalculate all

# FULL_RECALC equal to number after ':-' if parameter was not provided
FULL_RECALC=${1:-0}

rm -rf ./logs && mkdir -p ./logs

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

# dir_name and job_name are required in the next steps
dir_name=$(basename $(pwd))
step_num=$(echo ${dir_name} | cut -d'_' -f 1)
job_name="Step_${step_num}"

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))
gtflist=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --gtf_annotation_and_counter))

# loop over all samples
for ((i=0;i< ${#samples[@]} ;i+=2))
do
    # sample_dir=${samples[i]}
    sample_id=${samples[i]}
    description=${samples[i+1]}
    
    for ((j=0; j< ${#gtflist[@]} ;j+=2))
    do
        gtf_file=${gtflist[j]}
        counter=${gtflist[j+1]}

	## if file with counts already exist and FULL_RECALC=0, use the following directory
	gtf_dirname=$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export_gtf_by_name_and_counter ${gtf_file} ${counter} | grep OUTPUT_DIR | cut -d"=" -f2)
	result_file="${DATASET_DIR}/${sample_id}/${counter}/${gtf_dirname}/${sample_id}_${counter}.out"

	if [ -f ${result_file} -a ${FULL_RECALC} -eq 0 ]; then
	    
	    echo "Skip recalculation for ${sample_id} sample for dir ${gtf_dirname}. Results already obtained."
	    continue
	else
	    (set -x; qsub -N "${job_name}_${sample_id}" -P "${PROJECT}" -l h_rt="24:00:00" Extract_Counts.qsub ${sample_id} ${gtf_file} ${counter})
	fi
    done
    
done

echo "End of 8 step commands"
echo "-----------------------"

