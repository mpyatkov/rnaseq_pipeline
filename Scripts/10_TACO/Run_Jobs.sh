#!/bin/bash -l
##################################################################################
# Calculate TACO tracks
# 1. run_assembler.qsub: run stringtie on sample_id -> sample_id.gtf
# 2. go to directory with group/all sample_id gtf files and run TACO
#    meta-assembler -> GROUP_NAME.gtf
# 3. compare GROUP_NAME.gtf with ${TACO_REFCOMP_REFERENCE} using taco_refcomp
# 4. repeat steps 1-3 with on/off ${TACO_STRINGTIE_REFERENCE}

#Example: 
#./Run_Jobs.sh
##################################################################################

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

#Remove *.o files from previous jobs
rm -rf ./logs && mkdir -p ./logs

function calculate_taco_output() {

    local ref_flag=$1

    # checking that reference file is provided
    if [ ${TACO_STRINGTIE_REFERENCE} = "NA" -a ${ref_flag} = "RefOn" ]; then
	echo "ERROR: inconsistent options presented: TACO_STRINGTIE_REFERENCE=NA and REF_FLAG=RefOn"
	echo "       please provide reference file using TACO_STRINGTIE_REFERENCE option"
	echo "       in file Pipeline_Setup.conf"
	exit 1
    fi
    
    OUTDIR="Job_Summary_${ref_flag}"
    rm -rf ${OUTDIR} && mkdir -p ${OUTDIR}

    ## TACO is enable

    # 1. Run StringTie for each sample to get GTF file with isoforms.

    samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

    # loop over all samples
    for ((i=0;i< ${#samples[@]} ;i+=2));
    do
	sample_id=${samples[i]}
	(set -x; qsub -N "${job_name}_S_${ref_flag}_${sample_id}" -P "${PROJECT}" -l h_rt="01:00:00" run_assembler.qsub ${sample_id} ${OUTDIR} ${ref_flag})
    done

    # 2. wait until StringTie completes all tasks. After that run TACO
    # meta-assembly for each group separately and for mega-group (ALL)
    # which represents all samples together

    # 3. wait until TACO is completed and run TACO_REFCOMP to calculate
    # pairwise differences between TACO_REFCOMP_REFERENCE and groups
    
    # script creates processed TSV and GTF (optional) in ${OUTDIR}/${description}_${combined_name} directory
    # removes GTF files if TACO_GTF_KEEP=0 and TACO_REFERENCE_GTF=0
    
    # get names of all groups
    groups=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --groups))

    # add ALL group to groups
    # if groups length is 1, do not to add "ALL"
    if [ "${#groups[@]}" -ne 1 ]; then
	groups+=("ALL")
    fi

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
	

	(set -x; qsub -hold_jid "${job_name}_S_${ref_flag}_*" -N "${job_name}_TACO_${ref_flag}_${group}" \
		      -P "${PROJECT}" -l h_rt="01:00:00" run_taco.qsub "${group}" "${combined_name}" "${description}" "${OUTDIR}" "${ref_flag}")

    done
}

##### MAIN #####

MAIN_DIR="$(pwd)"
dir_name=$(basename $(pwd))
step_num=$(echo ${dir_name} | cut -d'_' -f 1)
job_name="Step_${step_num}"

# There are two ways how we can calculate STRINGTIE GTF files (with
# reference GTF and without). Both cases lead to different isoforms
# output. If we provide reference GTF in most cases we will get
# already known isoforms. If we will not provide reference GTF,
# StringTie will try to create de-novo assembly with the potential new
# isoforms. 

reference_options=( "RefOn" "RefOff")

for ref_flag in ${reference_options[@]}; do
    echo "STRINGTIE with ${ref_flag} flag"
    calculate_taco_output ${ref_flag}
done

echo "End of qsub commands"
echo "-----------------------"
##################################################################################
