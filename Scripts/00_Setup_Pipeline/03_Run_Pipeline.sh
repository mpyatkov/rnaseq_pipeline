#!/usr/bin/env bash

# prevent some unexpected behaviour
set -o errexit
# set -o pipefail
# set -o nounset

##################################################################################
#Andy Rampersaud, 03.14.16
#This 03_Run_Pipeline.sh script is used to run all steps in the pipeline
#The purpose of this script is to stream-line the pipeline so that 
#	1. A single command executes all steps sequentially
#	2. Minimize user intervention by automatically moving to the next step
#Assumptions for this script to work: 
#	1. This script is located in the "00_Setup_Pipeline" folder
#	2. The 00_Setup_Pipeline folder is in a "Scripts"
#	3. The Scripts folder contains all of the pipeline steps to run
#Way to run script:
#Usage: ./03_Run_Pipeline.sh <Start_Step>
#	<Start_Step> = Needs to be either "FULL" or a specific pipeline step 
#Example: 
#User wants to run the full pipeline:
#./03_Run_Pipeline.sh FULL
#User wants to run a subset of the pipeline starting at a specific step (04_TopHat_Paired_End):
#./03_Run_Pipeline.sh 04_TopHat_Paired_End
##################################################################################
#---------------------------------------------------------------------------------

# Activate anaconda environment
set +eu
module load anaconda2
source activate RNAseq
set -eu

# Process system argument:
if [[ $# -ne 1 ]]
then
    echo "Usage: ./03_Run_Pipeline.sh <Start_Step>"
    echo "<Start_Step> = Needs to be either \"FULL\" or a specific pipeline step"
    echo "See 03_Run_Pipeline.sh for details."	
    exit 0
fi

START_STEP=$1

# Job time START_TIME in seconds
START_TIME=$(date +"%s")

# export and print all variables from Pipeline_Setup.conf
eval "$(./01_Pipeline_Setup.py --export)"

# if full pipeline than generate DE directories
if [[ "${START_STEP}" == "FULL" ]]
then
    ./01_Pipeline_Setup.py --generate
fi
# TODO: check existence of 09* directories


# The Scripts folder is 1 level up from 00_Setup_Pipeline:
cd ..
STEPS_DIR=$(pwd)

# get all steps (directories with numbers in the beginning)
ALL_PIPELINE_STEPS=$(find . -maxdepth 1 -type d -name '[[:digit:]]*' | sort -n | xargs -n1 basename | sed -n '/^[0-9]/p' | sed '/00_Setup_Pipeline/d' | sed '/Generate_Tracks/d')

# is full pipeline execution? Yes by default
ISFULL=1

if [[ "${START_STEP}" == "FULL" ]]
then
    PIPELINE_STEPS=${ALL_PIPELINE_STEPS}
else
    # discard lines before START_STEPS
    PIPELINE_STEPS=$(printf "%s\n" "${ALL_PIPELINE_STEPS}" | sed -n -e '/'${START_STEP}'/,$p')
    ISFULL=0
fi

# START PIPELINE

#For loop over steps in the pipeline
for STEP in ${PIPELINE_STEPS}
do
    echo "Running: ${STEP} ..." 
    
    # cd "${STEP}"
    pushd "${STEP}"
    
    #Run_Jobs:
    # TODO: check out the exit code of each step to prevent running pipeline
    # if job return something wrong create error message and stop pipeline

    ./Run_Jobs.sh
    # Need a way to periodically check that all jobs have completed 
    # Count the number of jobs submitted by the user
    # Once the count goes to zero then summarize this job and move to the next step
    # Use the "${BU_USER}" variable from the 01_Pipeline_Setup.sh
    # Need to omit the 1st 2 lines of qstat command:
    # JOB_COUNT=$(qstat -u ${BU_USER} | awk 'FNR>2 {print $0}' | wc -l)
    # echo ${JOB_COUNT} 

    # Split by "_" grab the 1st part (use cut command)
    STEP_NUM=$(echo "${STEP}" | cut -d'_' -f 1)
    # Create the Job_Name:
    JOB_NAME="Step_${STEP_NUM}"
    
    # Use ${JOB_NAME} for qsub -N parameter
    # Better to check: qstat -r -u ${BU_USER}
    # The (-r) option gives the Full jobname
    # Extract lines with "Full jobname:"
    # Extract lines with step-specific name
    # Count the number of lines (number of running jobs)
    JOB_COUNT=$(qstat -r -u "${BU_USER}" | grep 'Full jobname:' | grep "${JOB_NAME}" | wc -l)

    echo "${BU_USER} running ${JOB_COUNT} jobs"
    
    echo "Periodically checking until jobs complete (please wait)..."

    # Use a while loop to check ${JOB_COUNT}
    while [[ ${JOB_COUNT} -ne 0 ]]
    do
	#Wait 01 minute then check ${JOB_COUNT} again
	sleep 1m
        
	JOB_COUNT=$(qstat -r -u ${BU_USER} | grep 'Full jobname:' | grep ${JOB_NAME} | wc -l)
    done
    
    # check if step failed
    number_of_runs=$(find ./logs/ -name "Step_*.o*" | wc -l)
    number_of_ok=$(find ./logs/ -name "Step_*.o*" | xargs grep -i "IAMOK" | wc -l)
    if [[ "${number_of_runs}" -ne "${number_of_ok}" ]]
    then
	echo "FOUND ERROR. CHECK LOGS ON STEP --> ${STEP} <--"
	exit 1
    fi

    # prevent immediate exit if ./Summarize_Jobs.sh does not exist
    set +eu
    ./Summarize_Jobs.sh
    set -eu

    printf "Done: %s\n\n -------------------------" "${STEP}"

    # cd ..
    popd
done

# For the 12_Generate_Tracks - I just want to run this step 
# Run the step 12_Generate_Tracks and automatically go back:
(
    STEP_DIR=$(find . -maxdepth 1 -type d -name '[[:digit:]]*' | xargs -n1 basename | sed -n '/Generate_Tracks/p')

    echo "Running: ${STEP_DIR}..." 
    cd "${STEP_DIR}"
    ./Generate_Tracks.sh
    printf "Done: %s\n\n -------------------------" "${STEP_DIR}"
)


echo 'Printing job duration for all steps...'
OUTPUT_FILE="${SETUP_PIPELINE_DIR}"/Pipeline_Runtime.txt

rm -rf "${OUTPUT_FILE}" && touch "${OUTPUT_FILE}"

# Print header to output file:
printf "Run time for each submitted job:\n" >> "${OUTPUT_FILE}"
printf "Job run times that deviate from the average should be inspected for possible errors (check the job log files)\n\n" >> "${OUTPUT_FILE}"

cd "${STEPS_DIR}"

# set +eu
# # scan all steps again to extract information about time
# for STEP in ${ALL_PIPELINE_STEPS}
# do
#     # Extract job runtimes and collect in output file
#     (
#         cd "${STEP}"
#         echo "${STEP}" >> "${OUTPUT_FILE}"
#         grep 'elapsed' *.o* >> "${OUTPUT_FILE}"
#     )
# done
# set -eu

#Also want to print the time to run this script:
echo '--------------------' >> "${OUTPUT_FILE}"
echo '03_Run_Pipeline.sh run time:' >> "${OUTPUT_FILE}"
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



