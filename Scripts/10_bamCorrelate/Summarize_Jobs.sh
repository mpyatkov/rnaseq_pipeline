#!/bin/bash
##################################################################################
# Andy Rampersaud, 02.16.16
#This script would be used to summarize bamCorrelate jobs
#Way to run script:
#Usage: 
#./bamCorrelate_Summary.sh
##################################################################################

set -o errexit
set -o pipefail
set -o nounset

#Minor issue with jobs using multiple cores
#Empty job output log files (*.pe* and *.po*) are created
#Remove them if they exist:
rm -rf *.po* *.pe*

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

SCRIPT_DIR="$(pwd)"

#---------------------------------------------------------------------------------
#Check that each variable prints a value to the terminal:
echo "-----------------------"
echo "Start of variable list:"
echo "-----------------------"
echo "DATASET_DIR:"
echo ${DATASET_DIR}
echo "DATASET_LABEL:"
echo ${DATASET_LABEL}
echo "SCRIPT_DIR:"
echo ${SCRIPT_DIR}
echo "TIME_LIMIT:"
echo ${TIME_LIMIT}
echo "-----------------------"
echo "End of variable list"
echo "-----------------------"

######################
OUTPUT_FILE=${SCRIPT_DIR}/Job_Summary/ENCODE_Blacklist_Overlap.txt
######################
if [ -f $OUTPUT_FILE ]
then 
    rm $OUTPUT_FILE
else
    touch $OUTPUT_FILE
fi
######################
#---------------------------------------------------------------------------------
cd ${SCRIPT_DIR}
cd Output
Job_Output_List=Output_bamCorrelate_*
for folder in ${Job_Output_List}
do
    #For the "Full_Genome" scan there is no corresponding Overlap_Summary.txt
    if [[ ${folder} == *"Full_Genome"* ]]; then
	#echo "Full_Genome folder"
	#Continue : means jump to the next element
	continue
    else
	echo ${folder}
	sed -n '1,5p' ${folder}/*_Output/Overlap_Summary.txt >> $OUTPUT_FILE
	echo >> $OUTPUT_FILE
    fi
    #End of if statement for "Full_Genome" check
done
cd ${SCRIPT_DIR}

echo 'Check out '${OUTPUT_FILE}
echo '#--------------------------------------------------------------------------'
##################################################################################
