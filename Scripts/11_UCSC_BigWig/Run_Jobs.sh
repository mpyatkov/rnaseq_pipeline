#!/bin/bash -l
##################################################################################
#Andy Rampersaud, 08.07.17
#This script would be used to run UCSC_BigWig.pbs in parallel for different FASTQ files
#Way to run script:
#Usage: ./Run_Jobs.sh
#Example: 
#./Run_Jobs.sh
##################################################################################

set -o errexit
set -o pipefail
set -o nounset

#Remove *.o files from previous jobs
rm -rf ./logs && mkdir -p ./logs

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

if [[ ${BIGWIG_ENABLE} == 0 ]]; then
    echo "BIGWIG_ENABLE=0"
    echo "BIGWIG files are not required. Exit. "
    exit 0
fi

if [ -f "./COMBINED_PAIRS.txt" ]; then
    rm COMBINED_PAIRS.txt
fi

SCRIPT_DIR="$(pwd)"

dir_name=$(basename $(pwd))
step_num=$(echo ${dir_name} | cut -d'_' -f 1)
job_name="Step_${step_num}"

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

# loop over all samples
for ((i=0;i< ${#samples[@]} ;i+=2));
do
    # sample_dir=${samples[i]}
    sample_id=${samples[i]}
    description=${samples[i+1]}

    # check if sample in db
    sample_info=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --get_sample_info ${sample_id}))
    # get path to the sample
    path_to_sample="${VM_DIR_UCSC}/INDEXED_PROJECTS/${sample_info[0]}/"

    # if Forward or Backward bigwig files do not exist - start calculation
    if [ ! -f "${path_to_sample}/${sample_id}.Forward.bw" -o ! -f "${path_to_sample}/${sample_id}.Reverse.bw" ]; then
	echo "start calc for ${sample_id}"
        (set -x; qsub -N "${job_name}_${sample_id}" -P "${PROJECT}" -l h_rt="${TIME_LIMIT}" Individual_BigWig.qsub ${sample_id})
    fi
done

# waiting for competion of previous jobs and starting combining files creation

# 1. we suppose that all individual bigwig files already located on server
# 2. we suppose that all combined files are from one project and have same
# name notations; for example group from two files G186_M1 and G186_M2
# will be combined to file G186_M1M2

# if at least one of combined does not exist run script which checks
# and recalculates required combined files
recalculate_combined=0

# remove file from previous run which will contains group names
if [ -f "COMBINED_PAIRS.txt" ]; then
    rm -rf COMBINED_PAIRS.txt
fi

groups=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --groups))
for group in "${groups[@]}"; do

    combined_name_project=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --combined_name_by_group ${group}))
    combined_name=${combined_name_project[0]}
    project=${combined_name_project[1]}
    description=${combined_name_project[2]}
    color=${combined_name_project[3]}
    path_on_server="${VM_DIR_UCSC}/INDEXED_PROJECTS/${project}"
    forward_bw="${combined_name}.Forward.combined.bw"
    reverse_bw="${combined_name}.Reverse.combined.bw"
    echo "Checking availability of ${combined_name} combined files inside the ${VM_DIR_UCSC}/INDEXED_PROJECTS/${project}"
    if [ ! -f "${path_on_server}/${forward_bw}" -o ! -f "${path_on_server}/${reverse_bw}" ]; then
	echo "${combined_name} not found. Recalculation is required."
	recalculate_combined=1
    fi

    server_dir_name="INDEXED_PROJECTS/${project}"
    # print "filename \t group_name for forward reads"
    printf "%s\t%s\t%s\t%s\n" "${forward_bw}" "${description}" "${server_dir_name}" "${color}" >> COMBINED_PAIRS.txt

    # print "filename \t group_name for reverse reads"
    printf "%s\t%s\t%s\t%s\n" "${reverse_bw}" "${description}" "${server_dir_name}" "${color}">> COMBINED_PAIRS.txt
    
done

if [ ${recalculate_combined} -eq 1 ]; then
    qsub -hold_jid "${job_name}*" -N "${job_name}_Combined" -P "${PROJECT}" -l h_rt="${TIME_LIMIT}" Combined_BigWig.qsub
fi

./UCSC_Tracks.sh

echo "End of qsub commands"
echo "-----------------------"
##################################################################################
