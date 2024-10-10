#!/usr/bin/env bash

# prevent some unexpected behaviour
set -o errexit
# set -o pipefail
# set -o nounset

##################################################################################
#Andy Rampersaud, 03.14.16
#This 02_Run_Pipeline.sh script is used to run all steps in the pipeline
#The purpose of this script is to stream-line the pipeline so that 
#	1. A single command executes all steps sequentially
#	2. Minimize user intervention by automatically moving to the next step
#Assumptions for this script to work: 
#	1. This script is located in the "00_Setup_Pipeline" folder
#	2. The 00_Setup_Pipeline folder is in a "Scripts"
#	3. The Scripts folder contains all of the pipeline steps to run
#Way to run script:
#Usage: ./02_Run_Pipeline.sh <OPTION>
#	<OPTION> = set of of options specified below
#Example: 
#User wants to run the full pipeline:
#./02_Run_Pipeline.sh FULL
##################################################################################
#---------------------------------------------------------------------------------

# Process system argument:
if [[ $# -ne 1 ]]
then
    echo "Usage: ./02_Run_Pipeline.sh <option>"
    echo "<option> = START - light version of FULL without recalculation already existed samples"
    echo "<option> = FULL - full pipeline run with recalculation of all steps"
    echo "<option> = DEONLY recalculate DE and summary directories (09abcd,12,13,99)"
    echo "<option> = TACO creates TACO tracks with isoforms (01,02,04 and 10 without full recalculation if that possible)"
    echo "<option> = BIGWIG creates UCSC tracks (01,02,04 and 11 without full recalculation if that possible)"
    echo "<option> = start_step (example. 05) start from specific step"
    echo "See 02_Run_Pipeline.sh for details."	
    exit 0
fi

# can be option or resume step number
START_STEP=$1

# Job time START_TIME in seconds
START_TIME=$(date +"%s")

# export and print all variables from Pipeline_Setup.conf
eval "$(./01_Pipeline_Setup.py --export)"

# The Scripts folder is 1 level up from 00_Setup_Pipeline:
SETUP_DIR=$(pwd)
cd ..
STEPS_DIR=$(pwd)

# log file for the whole process (also output to stdout)
LOG_FILE=${SETUP_DIR}/log.txt
rm -rf ${LOG_FILE}
# for explanation visit the following link:
# https://unix.stackexchange.com/a/145654
exec &> >(tee -a ${LOG_FILE})

# arguments parsing
# GENERATE: BOOL generate 09abcd folders 0/1
# RESUME: BOOL only for case when we need restart from particular step [0/1]
# INLCUDE: STRING list of steps which will be executed separated by "\|" 
# FULL_RECALC: 1/0 - if set to 0 that means that some steps
# can reuse already calculated results (ex. BAM files, Counts, Strandedness)

FULL_RECALC=0

if [[ "${START_STEP}" == "FULL" ]]; then
    GENERATE=1
    RESUME=0
    INCLUDE="ALL"
    FULL_RECALC=1
elif [[ "${START_STEP}" == "ONLYBAM" ]]; then
    GENERATE=1
    RESUME=0
    INCLUDE="01\|02\|03\|04"
    FULL_RECALC=1
elif [[ "${START_STEP}" == "START" ]]; then
    GENERATE=1
    RESUME=0
    INCLUDE="ALL"
elif [[ "${START_STEP}" == "DEONLY" ]]; then
    GENERATE=1
    RESUME=0
    INCLUDE="09\|12\|13\|14\|99"
elif [[ "${START_STEP}" == "TACO" ]]; then
    GENERATE=0
    RESUME=0
    INCLUDE="01\|02\|04\|10"
elif [[ "${START_STEP}" == "BIGWIG" ]]; then
    GENERATE=0
    RESUME=0
    INCLUDE="01\|02\|04\|11"
else
    # resume from step $START_STEP
    GENERATE=0
    RESUME=1
    INCLUDE="NONE"
fi

# generate 09abcd directories
if [[ $GENERATE -eq 1 ]]; then
    pushd ${SETUP_DIR}
    ./01_Pipeline_Setup.py --generate
    popd
fi

# get all steps (directories with numbers in the beginning)
ALL_PIPELINE_STEPS=$(find . -maxdepth 1 -type d -name '[[:digit:]]*' | sort -n | xargs -n1 basename | grep -Po "\K(^[0-9a-zA-Z]*)(?=_)" | uniq | sed '/00/d')

if [[ $INCLUDE == "ALL" ]]; then
    echo "running full pipeline"
    PIPELINE_STEPS=${ALL_PIPELINE_STEPS}
else
    if [[ $RESUME == 1 ]]; then
	echo "resume from step: ${START_STEP}"
	# resume from START_STEP
	PIPELINE_STEPS=$(printf "%s\n" "${ALL_PIPELINE_STEPS}" | sed -n -e '/'${START_STEP}'/,$p')
    else
	echo "calculate only steps: ${INCLUDE}"
	# run only required directories (VENNONLY, DEONLY options)
	PIPELINE_STEPS=$(printf "%s\n" "${ALL_PIPELINE_STEPS}" | grep $INCLUDE)
    fi
fi

# START THE PIPELINE

wait_until_completed() {
    # return nothing just wait until job completes
    # BU_USER is global var, exported from config

    local JOB_NAME=$1

    JOB_COUNT=$(qstat -r -u "${BU_USER}" | grep 'Full jobname:' | grep "${JOB_NAME}" | wc -l)

    echo "${BU_USER} running ${JOB_COUNT} jobs"
    echo "Periodically checking until jobs complete (please wait)..."

    # Use a while loop to check ${JOB_COUNT}
    while [[ ${JOB_COUNT} -ne 0 ]]
    do
	# Wait 01 minute then check ${JOB_COUNT} again
	sleep 1m
        
	# use qstat to extract information from server queue
	JOB_COUNT=$(qstat -r -u ${BU_USER} | grep 'Full jobname:' | grep ${JOB_NAME} | wc -l)
    done
}

logs_checking() {
    # check logs in current directory, exit with an error if "IAMOK" is not found
    # at the end of the log file

    number_of_runs=$(find ./logs/ -name "Step_*.o*" | wc -l)
    number_of_ok=$(find ./logs/ -name "Step_*.o*" | xargs grep -i "IAMOK" | wc -l)
    if [[ "${number_of_runs}" -ne "${number_of_ok}" ]]
    then
	echo "FOUND ERROR. CHECK LOGS ON STEP --> ${STEP} <--"
	exit 1
    fi
}

execute_wait_summarize_job(){
    # execute, wait, summarize the particular job

    local STEP_NUMBER=$1
    local FULL_RECALC=$2
    printf "Start: %s\n\n -------------------------\n" "${STEP_NUMBER}"

    search_body="${STEPS_DIR} -maxdepth 1 -name '${STEP_NUMBER}*' -type d"
    sub_steps=$(eval find ${search_body} | xargs -n1 basename)
    
    for sstep in ${sub_steps}; do
	pushd "${STEPS_DIR}/$sstep"
	./Run_Jobs.sh ${FULL_RECALC}
	popd
    done
    
    job_name="Step_${STEP_NUMBER}_${DATASET_LABEL}"
    
    # check cluster queue until the step completes
    wait_until_completed ${job_name}
    
    # checking and summarizing
    for sstep in ${sub_steps}; do
	pushd "${STEPS_DIR}/$sstep"
	
	# checking logs in current directory
	if [ -d "./logs" ]; then 
	    logs_checking
	else
	    echo "WARNING: logs directory does not exist. Analysis of logs will not be produced"
	fi
	
	# summarize if all logs contain the "IAMOK" at the end
	# prevent immediate exit if ./Summarize_Jobs.sh does not exist
	set +eu
	if [ -f Summarize_Jobs.sh ]; then
	    echo "Summarise results..."
	    ./Summarize_Jobs.sh ${FULL_RECALC}
	fi
	set -eu
	printf "Done: %s\n\n -------------------------\n" "${job_name}"
	popd
    done

    printf "END: %s\n\n -------------------------" "${STEP_NUMBER}"
}


#For loop over steps in the pipeline
for STEP in ${PIPELINE_STEPS}
do
    # steps consisting from multiple substeps will be run as one task
    # 09a -> 09a_1, 09a_2, 09a_3, ...
    execute_wait_summarize_job $STEP ${FULL_RECALC}
done

cd "${STEPS_DIR}"

echo 'Pipeline is done, check out your results!'
echo "=========================================================="
echo "Finished on : $(date)"

#Use to calculate job time:
#End_Time in seconds
END_TIME=$(date +"%s")
DIFF=$((END_TIME-START_TIME))
echo "$(($DIFF / 3600)) hours, $(((DIFF / 60) % 60)) minutes and $((DIFF % 60)) seconds elapsed."

