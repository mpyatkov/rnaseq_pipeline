#!/bin/bash

##################################################################################
#Andy Rampersaud, 02.22.16
#This script would be used to run CollectRnaSeqMetrics.pbs in parallel for different samples
#Way to run script:
#Usage: ./Run_Jobs.sh
#Example: 
#./Run_Jobs.sh 
##################################################################################

#Remove *.o files from previous jobs
rm -rf *.o* *.e*
##################################################################################
#---------------------------------------------------------------------------------

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

# Calculate strand_specificity from STRANDEDNESS
#CollectRnaSeqMetrics: strand_specificity
if [ ${STRANDEDNESS} -eq 0 ]
then
strand_specificity="NONE"

elif [ ${STRANDEDNESS} -eq 2 ]
then
strand_specificity="FIRST_READ_TRANSCRIPTION_STRAND"

elif [ ${STRANDEDNESS} -eq 1 ]
then
strand_specificity="SECOND_READ_TRANSCRIPTION_STRAND"
fi

# dir_name and job_name are required in the next steps
dir_name=$(basename $(pwd))
step_num=$(echo ${dir_name} | cut -d'_' -f 1)
job_name="Step_${step_num}"

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

# loop over all samples
for ((i=0;i< ${#samples[@]} ;i+=3));
do
    sample_dir=${samples[i]}
    sample_id=${samples[i+1]}
    description=${samples[i+2]}

    (set -x; qsub -N "${job_name}_${sample_id}" -P "${PROJECT}" -l h_rt="${TIME_LIMIT}" CollectRnaSeqMetrics.qsub "${sample_id}" "${strand_specificity}")
done

echo "End of STEP 6 commands"
echo "-----------------------"
##################################################################################
