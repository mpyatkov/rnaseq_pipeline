#!/bin/bash -l
##################################################################################
# Calculate TACO tracks
# 1. run_assembler.qsub: run stringtie on sample_id -> sample_id.gtf
# 2. go to directory with group/all sample_id gtf files and run TACO
#    meta-assembler -> GROUP_NAME.gtf

#Example: 
#./Run_Jobs.sh
##################################################################################

# FULL_RECALC equal to number after ':-' if parameter was not provided
FULL_RECALC=${1:-0}

set -o errexit
set -o pipefail
set -o nounset

#Remove *.o files from previous jobs
rm -rf ./logs && mkdir -p ./logs

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

# set +eu
# if [ -f "./COMBINED_PAIRS.txt" ]; then
#     rm COMBINED_PAIRS.txt
# set -eu

rm -rf Job_Summary && mkdir -p Job_Summary

MAIN_DIR="$(pwd)"

dir_name=$(basename $(pwd))
step_num=$(echo ${dir_name} | cut -d'_' -f 1)
job_name="Step_${step_num}"

function create_track_for_bigbed () {
    # THIS SCRIPT USES GLOBAL VARS
    local output_prefix=$1
    local project=$2
    local color=${3:-'0,0,0'}
    local description=${4:-"All samples together (see ${project}/${output_prefix}_description.txt)"}

    trackline="track type=bigBed name='${output_prefix}' description='${description} / ${output_prefix}' color=${color} visibility=full itemRgb=on bigDataUrl=http://waxmanlabvm.bu.edu/TRACKS/INDEXED_PROJECTS/${project}/${output_prefix}.bb"

    output_dir="${VM_DIR_UCSC}/PERSONAL/${BU_USER}/${DATASET_LABEL}/TACO_Track_Lines"
    mkdir -p ${output_dir}
    echo "${trackline}" > "${output_dir}/${description}_${output_prefix}.txt"
}

if [ "${TACO_MODE}" -eq 0 ]; then
    echo "TACO_MODE=0. Step is not required."
    exit 0
    
elif [ "${TACO_MODE}" -eq 1 ]; then

    echo "Recalculate combined file for all samples"
    combined_name_project=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --combined_name_by_group "ALL"))
    combined_name=${combined_name_project[0]}
    project=${combined_name_project[1]}
    description=${combined_name_project[2]} # just 'ALL'
    color=${combined_name_project[3]}

    path_on_server="${VM_DIR_UCSC}/INDEXED_PROJECTS/${project}"

    # For the case when we need to combine all samples in one TACO
    # meta-assembly I calculate hash for all sample_id names to
    # prevent long output name. The hash trimmed to the first 7
    # characters. After that I contatenate it with the project name
    # (ex. G186_1aeb2dc) to get the final name for file. Next time we
    # can check if a file with the same hash was created on the
    # server.
    
    # take hash from the combined_name and get only 7 first symbols
    hash=$(echo ${combined_name} | md5sum | cut -d" " -f1 | cut -c1-7)
    combined_name_hash="${project}_${hash}"
    
    echo "Combined name: ${combined_name_hash}"
    
    if [ ! -f "${path_on_server}/${combined_name_hash}.bb" -o ${FULL_RECALC} -eq 1 ]; then
	samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

	# loop over all samples
	for ((i=0;i< ${#samples[@]} ;i+=2));
	do
	    sample_id=${samples[i]}
	    (set -x; qsub -N "${job_name}_S_${sample_id}" -P "${PROJECT}" -l h_rt="01:00:00" run_assembler.qsub ${sample_id} ${combined_name_hash})
	done
	
	# calculate TACO meta-assembly, copy to server
	(set -x; qsub -hold_jid "${job_name}_S_*" -N "${job_name}_TACO_ALL" -P "${PROJECT}" -l h_rt="03:00:00" run_taco.qsub  "ALL" ${combined_name_hash} ${project})

	# To store information about what real name of 'hashed' output
	# file (which also contains info about samples components) will be 
	# created file with description (#_description.txt)
	
	mkdir -p ${combined_name_hash}
	printf "%s\t%s\t%s\n" ${BU_USER} "$(date)" "${combined_name}" > ${combined_name_hash}/${combined_name_hash}_description.txt

    else
	echo -e "\nTACO meta-assembly files already exist:"
	echo "bigBed file: ${path_on_server}/${combined_name_hash}.bb"
	echo "GTF file:    ${path_on_server}/${combined_name_hash}.gtf"
	echo "Description of files:"
	echo "${path_on_server}/${combined_name_hash}_description.txt"
	echo -e "Tracks will be presented in ${VM_DIR_UCSC}/PERSONAL/${BU_USER}/${DATASET_LABEL}/TACO_Track_Lines/\n"

	echo "Copy files from server to Job_Summary directory:"
	pushd Job_Summary
	mkdir ${combined_name_hash}
	pushd ${combined_name_hash}
	cp "${path_on_server}/${combined_name_hash}.bb" ./
	cp "${path_on_server}/${combined_name_hash}.gtf" ./
	cp "${path_on_server}/${combined_name_hash}_description.txt" ./
	popd
	popd
    fi
    
    # create track
    create_track_for_bigbed ${combined_name_hash} ${project} # other two parameters will be default 
    
else
    echo "check availability bigBed files for different groups"
    groups=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --groups))
    for group in "${groups[@]}"; do
	echo "Processing ${group} group"
	combined_name_project=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --combined_name_by_group ${group}))
	combined_name=${combined_name_project[0]}
	project=${combined_name_project[1]}
	condname=${combined_name_project[2]} # Condition_Name
	color=${combined_name_project[3]}
	
	path_on_server="${VM_DIR_UCSC}/INDEXED_PROJECTS/${project}"
	
	if [ ! -f "${path_on_server}/${combined_name}.bb" -o ${FULL_RECALC} -eq 1 ]; then
	    echo "Track recalculating..."
	    
	    samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples_by_group ${group} "color_enable"))
	    
	    # loop over all samples
	    for ((i=0;i< ${#samples[@]} ;i+=3));
	    do
		sample_id=${samples[i]}
		# condname=${samples[i+1]}
		# color=${samples[i+2]}
	
		(set -x; qsub -N "${job_name}_${group}_${sample_id}" -P "${PROJECT}" -l h_rt="01:00:00" run_assembler.qsub ${sample_id} ${combined_name})
	    done

	    # calculate TACO meta-assembly, copy to server
	    (set -x; qsub -hold_jid "${job_name}_${group}_*" -N "${job_name}_TACO_${group}" -P "${PROJECT}" -l h_rt="01:00:00" run_taco.qsub  ${group} ${combined_name} ${project})
	    
	else
	    echo -e "\nTACO meta-assembly files already exist:"
	    echo "bigBed file: ${path_on_server}/${combined_name}.bb"
	    echo "GTF file:    ${path_on_server}/${combined_name}.gtf"
	    echo -e "Tracks will be presented in ${VM_DIR_UCSC}/PERSONAL/${BU_USER}/${DATASET_LABEL}/TACO_Track_Lines/\n"

	    echo "Copy files from server to Job_Summary directory"
	    pushd Job_Summary
	    mkdir ${combined_name}
	    pushd ${combined_name}
	    cp "${path_on_server}/${combined_name}.bb" ./
	    cp "${path_on_server}/${combined_name}.gtf" ./
	    popd
	    popd
	fi

	# create track line and copy to server
	create_track_for_bigbed ${combined_name} ${project} ${color} ${condname}
	
	# echo "${path_on_server}/${project}/${combined_name}.bb"

    done

fi

echo "End of qsub commands"
echo "-----------------------"
##################################################################################
