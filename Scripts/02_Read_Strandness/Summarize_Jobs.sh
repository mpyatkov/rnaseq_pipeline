#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

##################################################################################
#Andy Rampersaud, 02.23.16
#This script would be used to summarize Read_Strandness statistics 
##################################################################################

# ./Summarize_Jobs.sh 0 -- allows to reuse previously calculated results
# ./Summarize_Jobs.sh 1 -- recalculate all

# FULL_RECALC equal to number after ':-' if parameter was not provided
FULL_RECALC=${1:-0}

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

# output directory
output_dir="$(pwd)/Job_Summary"
rm -rf "${output_dir}" && mkdir -p "${output_dir}"

output_file="${output_dir}/Read_Strandness_Stats.txt"
rm -rf "${output_file}" && touch "${output_file}"

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

# loop over samples
for ((i=0;i< ${#samples[@]} ;i+=2));
do
    sample_id=${samples[i]}
    # description=${samples[i+2]}

    # change directory, save location
    pushd "${DATASET_DIR}/${sample_id}/Read_Strandness"
    
    #Print $Sample_ID
    echo "${sample_id}:" >> $output_file
    
    #Concatenate file to $OUTPUT_FILE
    (set -x; cat "${sample_id}_Read_Strandness.txt" >> "$output_file")
    
    #Add a return line
    echo >> $output_file

    # go back
    popd
done
echo "--------------------"

