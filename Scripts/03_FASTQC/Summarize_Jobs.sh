#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

##################################################################################
#Andy Rampersaud, 07.19.17
#This script would be used to summarize FASTQC jobs
##################################################################################

# ./Summarize_Jobs.sh 0 -- allows to reuse previously calculated results
# ./Summarize_Jobs.sh 1 -- recalculate all

# FULL_RECALC equal to number after ':-' if parameter was not provided
FULL_RECALC=${1:-0}

# Skip this step if recalculation flag set to 0
if [ ${FULL_RECALC} -eq 0 ]; then
    echo "Recalculation is not required. FULL_RECALC set to 0."
    exit 0
fi

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

#Extract the folder name:
dir_name=$(pwd)

output_dir=${dir_name}/Job_Summary

# output file for 03_Read_length
output_file="${output_dir}/Read_Length_Stats.txt"
# header for 03_Read_length
printf "SAMPLE_ID\tFASTQ_File_Name\tRead_Length(bp)\tRead_Count\n" >> "${output_file}"

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

for ((i=0;i< ${#samples[@]} ;i+=2));
do
    sample_id=${samples[i]}
    # description=${samples[i+1]}
    
    # combine all sample in one file (from step 03_Read_Length)
    cat "${output_dir}/${sample_id}_read_length.txt" >> $output_file

    # get_sample_info return (PRJ_NAME, READ1, READ2) or EXCEPTION
    sample_info=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --get_sample_info ${sample_id}))
    project_name=${sample_info[0]}

    # create remote fastqc directory if not created
    remote_fastqc_path="${VM_DIR_UCSC}/INDEXED_PROJECTS/${project_name}/FASTQC/"
    mkdir -p ${remote_fastqc_path}
        
    set +eu
    cp "${output_dir}/${sample_id}_FASTQC/"*_fastqc.html "${remote_fastqc_path}"
    if [ $? -ne 0 ]; then
        printf "ERROR: can't copy ${output_dir}/${sample_id}_FASTQC/ to ${remote_fastqc_path}. Probably file already exist.\n\n"
    fi
    set -eu
done

echo '#--------------------------------------------------------------------------'
echo 'Check the URL for the waxmanlabvm HTML files!'
echo 'Load the following link:'
echo 'http://waxmanlabvm.bu.edu/TRACKS/INDEXED_PROJECTS/'
echo '#--------------------------------------------------------------------------'
##################################################################################
