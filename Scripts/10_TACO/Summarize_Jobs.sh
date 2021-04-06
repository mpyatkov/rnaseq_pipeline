#!/bin/bash

# copy bigBed and description file to server
# copy TACO track to server
# clean Job_Summary directory from the individual GTF files

set -o errexit
set -o pipefail
set -o nounset

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

# checking if TACO is enabled
if [ "${TACO_ENABLE}" -eq 0 ]; then
    echo "TACO_ENABLE=0. Step is not required."
    exit 0
fi

function create_description(){
    # creates description for group of samples with user,date and
    # samples included in group

    local group=$1
    local outname_path=$2

    printf "user: %s\ndate: %s\n" "${BU_USER}" "$(date)" > ${outname_path}

    samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples_by_group ${group}))

    # loop over all samples
    for ((i=0;i< ${#samples[@]} ;i+=2));
    do
	sample_id=${samples[i]}
	condname=${samples[i+1]}
	printf "%s %s\n" "${sample_id}" "${condname}" >> ${outname_path}
    done
}

function create_track_for_bigbed () {
    # THIS SCRIPT USES GLOBAL VARS
    local output_prefix=$1
    local project=$2
    local color=${3:-'0,0,0'}
    local description=${4:-"All samples together (see ${project}/${output_prefix}_description.txt)"}

    trackline="track type=bigBed name='${output_prefix}' description='${description} / ${output_prefix}' color=${color} visibility=full itemRgb=on bigDataUrl=http://waxmanlabvm.bu.edu/TRACKS/INDEXED_PROJECTS/${project}/${output_prefix}.bb"

    echo "${trackline}"
}

# get names of all groups
groups=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --groups))

# add ALL group to groups
# if groups length is 1, do not to add "ALL"
if [ "${#groups[@]}" -ne 1 ]; then
    groups+=("ALL")
fi

OUTDIR="Job_Summary"
TRACKS_OUT="${OUTDIR}/${DATASET_LABEL}_TACO_Tracks.txt"
rm -rf ${TRACKS_OUT}

for group in "${groups[@]}"; do

    combined_name_project=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --combined_name_by_group ${group}))
    combined_name=${combined_name_project[0]}
    project=${combined_name_project[1]}
    description=${combined_name_project[2]} # 'ALL'/ Condition_Name
    color=${combined_name_project[3]}
    path_on_server="${VM_DIR_UCSC}/INDEXED_PROJECTS/${project}"
    
    if [ ${group} = "ALL" ];then
	hash=$(echo ${combined_name} | md5sum | cut -d" " -f1 | cut -c1-4)
	combined_name="${project}_${hash}"
    fi
    
    # create description for each group
    description_path="${OUTDIR}/${description}_${combined_name}"
    (set -x; mkdir -p ${description_path})
    create_description "${group}" "${description_path}/${combined_name}_description.txt"
    
    # add track to track_output.txt
    if [ ${group} = "ALL" ]; then
	trackline=$(create_track_for_bigbed ${combined_name} ${project} ${color} )
    else
	trackline=$(create_track_for_bigbed ${combined_name} ${project} ${color} ${description})
    fi

    echo "${trackline}" >> ${TRACKS_OUT}

    
    # copy bigBed and description files if it is required
    if [ ! -f "${VM_DIR_UCSC}/INDEXED_PROJECTS/${project}/${combined_name}.bb" ]; then
	echo "copy files to server"
	(set -x; cp "${OUTDIR}/${description}_${combined_name}/${description}_${combined_name}.bb" "${VM_DIR_UCSC}/INDEXED_PROJECTS/${project}/${combined_name}.bb")
	(set -x; cp "${OUTDIR}/${description}_${combined_name}/${combined_name}_description.txt" "${VM_DIR_UCSC}/INDEXED_PROJECTS/${project}/")
    fi

done

# copy tracks to server
(set -x; mkdir -p "${VM_DIR_UCSC}/PERSONAL/${BU_USER}/${DATASET_LABEL}/UCSC_Track_Lines")
(set -x; cp "${OUTDIR}/${DATASET_LABEL}_TACO_Tracks.txt" "${VM_DIR_UCSC}/PERSONAL/${BU_USER}/${DATASET_LABEL}/UCSC_Track_Lines")

# remove individual GTF files which is not required anymore
(set -x; rm -rf ${OUTDIR}/*.gtf)
