#!/bin/bash
##################################################################################
#Andy Rampersaud, 07.19.17
#This script would be used to summarize UCSC_BigWig statistics 
#Way to run script:
#Usage: 
#./UCSC_BigWig_Summary.sh
##################################################################################
#---------------------------------------------------------------------------------

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

#---------------------------------------------------------------------------------
#Check that each variable prints a value to the terminal:
echo "-----------------------"
echo "Start of variable list:"
echo "-----------------------"
echo "DATASET_DIR:"
echo ${DATASET_DIR}
echo "Sample_Labels_DIR:"
echo ${SAMPLE_LABELS_DIR}
echo "VM_DIR_UCSC:"
echo ${VM_DIR_UCSC}
echo "SCRIPT_DIR:"
echo ${SCRIPT_DIR}
echo "TIME_LIMIT:"
echo ${TIME_LIMIT}
echo "-----------------------"
echo "End of variable list"
echo "-----------------------"
#---------------------------------------------------------------------------------
##################################################################################
INPUT_DIR=$(pwd)
cd $INPUT_DIR
################################################
#A text file (Sample_Labels.txt) is needed to run this script
SCRIPT_DIR=$(pwd)

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

# loop over all samples
for ((i=0;i< ${#samples[@]} ;i+=3));
do
    Sample_DIR=${samples[i]}
    Sample_DIR=${samples[i+1]}
    Description=${samples[i+2]}
    #---------------------------
    ##Check that text file is read in properly:
    #echo 'Sample_DIR:'
    Sample_DIR=${myArray[0]}
    #echo 'Sample_ID:'
    Sample_ID=${myArray[1]}
    #echo $Sample_ID
    #echo 'Description:'
    Description=${myArray[2]}
    #echo $Description
    #---------------------------
    echo
    echo $Sample_ID
    #Copy over BAM and *.bw files:
    echo 'Copying over BAM, bigWig, and bigBed files to waxmanlabvm...'
    cp ${DATASET_DIR}/${Sample_ID}/fastq/tophat2/$Sample_ID'_primary_unique.bam' ${VM_DIR_UCSC}
    cp ${DATASET_DIR}/${Sample_ID}/fastq/tophat2/$Sample_ID'_primary_unique.bam.bai' ${VM_DIR_UCSC}
    cp ${DATASET_DIR}/${Sample_ID}/fastq/tophat2/UCSC_BigWig/*'.bw' ${VM_DIR_UCSC}
done 
cd $INPUT_DIR

##################################################################################
#waxmanlabvm:
#If you want the number of files in a dir (check that all files transferred over)
#ls -1 | wc -l
##################################################################################

echo 'UCSC files for data visualization should now be copied to the VM.'
echo '#--------------------------------------------------------------------------'
##################################################################################
