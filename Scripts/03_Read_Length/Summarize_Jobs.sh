#!/bin/bash
##################################################################################
#Andy Rampersaud, 02.23.16
#This script would be used to summarize Read_Length jobs
#Way to run script:
#Usage: 
#./Read_Length_Summary.sh
##################################################################################

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

# retrieve the dir name for this step:
dir_name=$(basename $(pwd))

# output dir to store text files:
output_dir="$(pwd)/Job_Summary"
output_file="${output_dir}/Read_Length_Stats.txt"

rm -rf ${output_file} && touch ${output_file}

cd "${output_dir}"

# print headers to the output file
printf "SAMPLE_ID\tFASTQ_File_Name\tRead_Length(bp)\tRead_Count\n" >> "${output_file}"

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

# rename directories from sample_dir to sample_id
for ((i=0;i< ${#samples[@]} ;i+=3));
do
    # sample_dir=${samples[i]}
    sample_id=${samples[i+1]}
    # description=${samples[i+2]}
    echo ${sample_id}
    cat "${sample_id}_read_length.txt" >> $output_file
    
done

echo "Check out ${OUTPUT_FILE}"
echo "#-----------------------------------------"

