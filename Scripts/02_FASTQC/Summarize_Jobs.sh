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
dir_name=$(pwd)

output_dir=${dir_name}/Job_Summary

# output file for 03_Read_length
output_file="${output_dir}/Read_Length_Stats.txt"
# header for 03_Read_length
printf "SAMPLE_ID\tFASTQ_File_Name\tRead_Length(bp)\tRead_Count\n" >> "${output_file}"

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

# rename directories from sample_dir to sample_id

# checking existing VM_DIR_FASTQC
# do not exit if an error occurs
set +eu
if [[ -d ${VM_DIR_FASTQC} ]]
then
    echo "VM_DIR_FASTQC exists!"
else
    echo "Attempt to create ${VM_DIR_FASTQC} directory..."
    mkdir -p "${VM_DIR_FASTQC}"
    if [ $? -ne 0 ]; then
        printf "\nWARNING: Ask somebody with rights to give you write access to $VM_DIR_FASTQC directory\n"
    else
        printf "\nDirectory ${VM_DIR_FASTQC} created\n"
    fi
fi
set -eu

for ((i=0;i< ${#samples[@]} ;i+=3));
do
    # sample_dir=${samples[i]}
    sample_id=${samples[i+1]}
    # description=${samples[i+2]}

    # combine all sample in one file (from step 03_Read_Length)
    cat "${output_dir}/${sample_id}_read_length.txt" >> $output_file

    # copy all fastqc files to VM_DIR_FASTQC webserver directory 
    # do not exit if can't copy files
    set +eu
    cp "${output_dir}/${sample_id}_FASTQC/"*_fastqc.html "${VM_DIR_FASTQC}"
    if [ $? -ne 0 ]
    then
        printf "ERROR: can't copy ${output_dir}/${sample_id}_FASTQC/ to ${VM_DIR_FASTQC}\n\n"
    fi
    set -eu
done

echo '#--------------------------------------------------------------------------'
echo 'Check the URL for the waxmanlabvm HTML files!'
echo 'Load the following link:'
echo 'Note: use full paths (instead of using tilde notation)'
echo 'http://waxmanlabvm.bu.edu/waxmanlab/FASTQC/'
echo '#--------------------------------------------------------------------------'
##################################################################################
