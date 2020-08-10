#!/bin/bash

# prevent some unexpected behaviour
set -o errexit
set -o pipefail
set -o nounset

##################################################################################
# Andy Rampersaud, 02.16.16
# This script would be used to rename folders
# Way to run script:
# Usage: ./Rename_Folders.sh
# Example:
# ./Rename_Folders.sh
##################################################################################

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

printf "Start renaming folders\n\n"

# rename directories from sample_dir to sample_id
for ((i=0;i< ${#samples[@]} ;i+=3));
do
    sample_dir=${samples[i]}
    sample_id=${samples[i+1]}
    # description=${samples[i+2]}

    printf "Rename %s to %s\n" "${sample_dir}" "${sample_id}"
    mv "${DATASET_DIR}/${sample_dir}" "${DATASET_DIR}/${sample_id}"
done

printf "End renaming folders\n\n"
