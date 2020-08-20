#!/bin/bash -l

exit 0
##################################################################################
#Andy Rampersaud, 02.22.16
#This script would be used to run bamCorrelate.qsub 
#Way to run script:
#Usage: ./Run_Jobs.sh
#Example: 
#./Run_Jobs.sh
##################################################################################
set -o errexit
set -o pipefail
set -o nounset

#Remove *.o* files from previous jobs

rm -rf *.o* *.e* *.pe* *.po*

##################################################################################
# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

(
    set +eu
    module load anaconda2
    source activate RNAseq
)
#bamCorrelate
zMin=0.6
zMax=1.0

SCRIPT_DIR="$(pwd)"

#Check that each variable prints a value to the terminal:
echo "-----------------------"
echo "Start of variable list:"
echo "-----------------------"
echo "DATASET_DIR:"
echo ${DATASET_DIR}
echo "DATASET_LABEL:"
echo ${DATASET_LABEL}
echo "zMin:"
echo ${zMin}
echo "zMax:"
echo ${zMax}
echo "SCRIPT_DIR:"
echo ${SCRIPT_DIR}
echo "TIME_LIMIT:"
echo ${TIME_LIMIT}
echo "-----------------------"
echo "End of variable list"
echo "-----------------------"

#---------------------------------------------------------------------------------
OUTPUT_DIR=$SCRIPT_DIR/Job_Summary
rm -rf $OUTPUT_DIR && mkdir -p $OUTPUT_DIR

#---------------------------------------------------------------------------------
echo "-----------------------"
echo "Start of qsub commands:"
echo "-----------------------"
#--------------------------------------------------------------------------------
#Create a job name that's a function of the folder name:

#Extract the folder name:
DIR_name=`basename ${SCRIPT_DIR}`

#Split by "_" grab the 1st part (use cut command)
Step_Num=$(echo ${DIR_name} | cut -d'_' -f 1)

#Create the Job_Name:
#Use ${Job_Name} for qsub -N parameter
Job_Name=$(echo 'Step_'${Step_Num})

#Need to check if Input/BED_Regions files are present
#Count the number of BED files (same way as above for log files):
BED_File_Count=`ls -1 ${SCRIPT_DIR}/Input/BED_Regions/*.bed 2>/dev/null  | wc -l`

echo "-----------------------"
echo 'BED_File_Count:'
echo ${BED_File_Count}
echo "-----------------------"
#---------------------------------------------------------------------------------
#bamCorrelate computes the overall similarity between two or more BAM
#files based on read coverage (number of reads) within genomic regions.
#The correlation analysis is performed for the entire genome by running
#the program in 'bins' mode, or for certain user selected regions in 'BED-file'
#mode. Because the computation of the coverage is time consuming the program
#outputs an intermediary file that can then be used with the 'plotCorrelation' tool
#for visualizing the correlation.
#---------------------------------------------------------------------------------
#If <bin_mode="On">: correlation analysis is performed for the entire genome
#If <bin_mode="Off">: correlation analysis is performed for user selected regions in 'BED-file'
#---------------------------------------------------------------------------------

#If the user did not provide BED files (${BED_File_Count} = 0): just submit the full genome scan job:
if [ ${BED_File_Count} == 0 ]
then
    echo
    echo 'The user did not provide Input/BED_Regions/*.bed files so job(s) sumitted:'
    echo '1. full genome scan job'
    echo

    #Initialize bin_mode variable:
    bin_mode="On"
    echo "-----------------------"
    echo "bin_mode:"
    echo ${bin_mode}
    echo "-----------------------"
    #Now use arguments in the qsub script call:
    (set -x; qsub -N ${Job_Name}'_bamCor_Full' -P wax-dk -l h_rt=${TIME_LIMIT} bamCorrelate.qsub ${DATASET_DIR} ${DATASET_LABEL} ${bin_mode} ${SCRIPT_DIR} ${OUTPUT_DIR})

fi
#End of ${BED_File_Count} if statement

#If the user did provide BED files (${BED_File_Count} != 0): submit a job for each BED file but also submit the the full genome scan job:
if [ ${BED_File_Count} != 0 ]
then
    echo
    echo 'The user did provide Input/BED_Regions/*.bed files so job(s) sumitted:'
    echo '1. full genome scan job'
    echo '2. separate job for each BED file'
    echo

    #Initialize bin_mode variable:
    bin_mode="On"
    echo "-----------------------"
    echo "bin_mode:"
    echo ${bin_mode}
    echo "-----------------------"

    #Now use arguments in the qsub script call:
    (set -x; qsub -N ${Job_Name}'_bamCor_Full' -P wax-dk -l h_rt=${TIME_LIMIT} bamCorrelate.qsub ${DATASET_DIR} ${DATASET_LABEL} ${bin_mode} ${SCRIPT_DIR} ${OUTPUT_DIR} ${zMin} ${zMax})

    #Initialize bin_mode variable:
    bin_mode="Off"
    echo "-----------------------"
    echo "bin_mode:"
    echo ${bin_mode}
    echo "-----------------------"

    #Need a COUNTER variable to generate an output log per job
    COUNTER=0
    BED_INPUT_DIR=$SCRIPT_DIR/Input/BED_Regions
    BED_List=$BED_INPUT_DIR/*
    for BED_file in $BED_List
    do
        COUNTER=$((COUNTER+1))

        #Format COUNTER to have 2 digits
        JOB_COUNTER=$(printf "%0*d\n" 2 $COUNTER)
        echo "-----------------------"
        echo 'JOB_COUNTER:'
        echo "-----------------------"
        echo ${JOB_COUNTER}
        BED_file_name_BED=`basename $BED_file`
        BED_file_name=${BED_file_name_BED%\.bed}
        echo ${BED_file_name}
        #---------------------------

        #Note:
        #If you need more than 9 command line arguments, you can use the shift command. 
        #See details in the bamCorrelate.qsub

        #Now use arguments in the qsub script call:
	(set -x; qsub -N ${Job_Name}'_bamCor_'${JOB_COUNTER} -P wax-dk -l h_rt=${TIME_LIMIT} bamCorrelate.qsub ${DATASET_DIR} ${DATASET_LABEL} ${bin_mode} ${SCRIPT_DIR} ${OUTPUT_DIR} ${zMin} ${zMax} ${JOB_COUNTER} ${BED_file_name})
        #End BED_file for loop
    done
fi
#End of ${BED_File_Count} if statement
#---------------------------------------------------------------------------------
echo "-----------------------"
echo "End of qsub commands"
echo "-----------------------"
##################################################################################
