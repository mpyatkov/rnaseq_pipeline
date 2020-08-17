#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

##################################################################################
#Andy Rampersaud, 07.19.17
#This script would be used to summarize FASTQC jobs
#Way to run script:
#Usage: 
#./FASTQC_Summary.sh
##################################################################################

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

#Extract the folder name:
dir_name=$(basename $(pwd))

output_dir=${dir_name}/Job_Summary

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

# rename directories from sample_dir to sample_id
for ((i=0;i< ${#samples[@]} ;i+=3));
do
    sample_dir=${samples[i]}
    sample_id=${samples[i+1]}
    # description=${samples[i+2]}
    # do not exit if can't copy files
    set +x
    cp "${output_dir}/${sample_id}_FASTQC/"*_fastqc.html "${VM_DIR_FASTQC}"
    if [ $? -ne 0 ]
    then
        printf "ERROR: can't copy ${output_dir}/${sample_id}_FASTQC/ to ${VM_DIR_FASTQC}\n\n"
    fi
    set -x
done

echo '#--------------------------------------------------------------------------'
echo 'Check the URL for the waxmanlabvm HTML files!'
echo 'Load the following link:'
echo 'Note: use full paths (instead of using tilde notation)'
echo 'http://waxmanlabvm.bu.edu/waxmanlab/FASTQC/'
echo '#--------------------------------------------------------------------------'
##################################################################################
