#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

##################################################################################
#Andy Rampersaud, 02.23.16
#This script would be used to summarize Read_Strandness statistics 
#Way to run script:
#Usage: 
#./Read_Strandness_Summary.sh
##################################################################################

# FULL_RECALC=1 by default if nothing provided
FULL_RECALC=${1:-0}

# Skip this step if recalculation flag set to 0
if [ ${FULL_RECALC} -eq 0 ]; then
    echo "Recalculation is not required. FULL_RECALC set to 0."
    exit 0
fi

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

