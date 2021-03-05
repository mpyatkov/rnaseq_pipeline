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
#Usage: ./02_Run_Pipeline.sh <Start_Step>
#	<Start_Step> = Needs to be either "FULL" or a specific pipeline step 
#Example: 
#User wants to run the full pipeline:
#./02_Run_Pipeline.sh FULL
#User wants to run a subset of the pipeline starting at a specific step (04_TopHat_Paired_End):
#./02_Run_Pipeline.sh 04_TopHat_Paired_End
##################################################################################
#---------------------------------------------------------------------------------

# Process system argument:
if [[ $# -ne 1 ]]
then
    echo "Usage: ./02_Run_Pipeline.sh <option>"
    echo "<option> = FULL - full pipeline run"
    echo "<option> = start_step (example. 05) start from specific step"
    echo "<option> = DEONLY recalculate DE and summary directories (09abcd,12,13,14)"
    echo "<option> = VENNONLY recalculate VENN (12)"
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

# arguments parsing
# GENERATE: BOOL generate 09abcd folders 0/1
# RESUME: BOOL only for case when we need restart from particular step [0/1]
# INLCUDE: STRING list of steps which will be executed separated by "\|" 
# FULL_RECALC: 1/0 - if set to 0 that means that some steps
# can reuse already calculated results (ex. BAM files, Counts, Strandedness)
FULL_RECALC=1

if [[ "${START_STEP}" == "FULL" ]]; then
    GENERATE=1
    RESUME=0
    INCLUDE="ALL"
if [[ "${START_STEP}" == "REFRESH" ]]; then
    GENERATE=1
    RESUME=0
    INCLUDE="ALL"
    FULL_RECALC=0
elif [[ "${START_STEP}" == "DEONLY" ]]; then
    GENERATE=1
    RESUME=0
    INCLUDE="09\|12\|13\|14"
elif [[ "${START_STEP}" == "VENNONLY" ]]; then
    GENERATE=0
    RESUME=0
    INCLUDE="12"
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
    printf "Start: %s\n\n -------------------------" "${STEP_NUMBER}"

    search_body="${STEPS_DIR} -maxdepth 1 -name '${STEP_NUMBER}*' -type d"
    sub_steps=$(eval find ${search_body} | xargs -n1 basename)
    
    for sstep in ${sub_steps}; do
	

	pushd "${STEPS_DIR}/$sstep"
	./Run_Jobs.sh ${FULL_RECALC}
	popd
    done
    
    job_name="Step_${STEP_NUMBER}"
    
    # check cluster queue until the step completes
    wait_until_completed ${job_name}
    
    # checking and summarizing
    for sstep in ${sub_steps}; do
	pushd "${STEPS_DIR}/$sstep"
	
	# checking logs in current directory
	logs_checking
	
	# summarize if all logs contain the "IAMOK" at the end
	# prevent immediate exit if ./Summarize_Jobs.sh does not exist
	set +eu
	./Summarize_Jobs.sh ${FULL_RECALC}
	set -eu
	printf "Done: %s\n\n -------------------------" "${job_name}"
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

echo 'Printing job duration for all steps...'
OUTPUT_FILE="${SETUP_PIPELINE_DIR}"/Pipeline_Runtime.txt

rm -rf "${OUTPUT_FILE}" && touch "${OUTPUT_FILE}"

# Print header to output file:
printf "Run time for each submitted job:\n" >> "${OUTPUT_FILE}"
printf "Job run times that deviate from the average should be inspected for possible errors (check the job log files)\n\n" >> "${OUTPUT_FILE}"

cd "${STEPS_DIR}"

#Also want to print the time to run this script:
echo '--------------------' >> "${OUTPUT_FILE}"
echo '02_Run_Pipeline.sh run time:' >> "${OUTPUT_FILE}"
echo "Check out: ${OUTPUT_FILE}"

echo 'Pipeline is done, check out your results!'
echo "=========================================================="
echo "Finished on : $(date)"

#Use to calculate job time:
#End_Time in seconds
END_TIME=$(date +"%s")
DIFF=$((END_TIME-START_TIME))
echo "$(($DIFF / 3600)) hours, $(((DIFF / 60) % 60)) minutes and $((DIFF % 60)) seconds elapsed."
echo "$((DIFF / 3600)) hours, $(((DIFF / 60) % 60)) minutes and $((DIFF % 60)) seconds elapsed." >> "${OUTPUT_FILE}"



