#!/bin/bash
##################################################################################
#Andy Rampersaud, 07.19.17
#This script would be used to summarize UCSC_BigWig statistics 
#Way to run script:
#Usage: 
#./UCSC_BigWig_Summary.sh
##################################################################################
#---------------------------------------------------------------------------------

set -o errexit
set -o pipefail
set -o nounset

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

if [[ ${BIGWIG_ENABLE} == 0 ]]; then
    echo "BIGWIG_ENABLE=0"
    echo "BIGWIG files are not required. Exit. "
    exit 0
fi

#---------------------------------------------------------------------------------
#Check that each variable prints a value to the terminal:
echo "-----------------------"
echo "Start of variable list:"
echo "-----------------------"
echo "VM_DIR_UCSC:"
echo ${VM_DIR_UCSC}
echo "-----------------------"
echo "End of variable list"
echo "-----------------------"

##################################################################################
# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

set +eu
# checking existing VM_DIR_UCSC
if [[ -d ${VM_DIR_UCSC} ]]
then
    echo "VM_DIR_UCSC_FASTQC exists..."
else
    echo "Need to create VM_DIR_UCSC!"
    mkdir -p ${VM_DIR_UCSC}
    if [[ $? -ne 0 ]] ; then
        printf "\nWARNING: Ask somebody with rights to give you write access to $VM_DIR_UCSC directory\n\n"
    else
        echo "Directory ${VM_DIR_UCSC} created."
    fi
fi
set -eu

# loop over all samples
for ((i=0;i< ${#samples[@]} ;i+=3));
do
    Sample_DIR=${samples[i]}
    Sample_ID=${samples[i+1]}
    Description=${samples[i+2]}
    #Copy over BAM and *.bw files:
    echo 'Copying over BAM, bigWig, and bigBed files to waxmanlabvm...'
    set +eu
    cp ${DATASET_DIR}/${Sample_ID}/aligner/$Sample_ID'_primary_unique.bam' ${VM_DIR_UCSC}
    cp ${DATASET_DIR}/${Sample_ID}/aligner/$Sample_ID'_primary_unique.bam.bai' ${VM_DIR_UCSC}
    cp ${DATASET_DIR}/${Sample_ID}/UCSC_BigWig/*.bw ${VM_DIR_UCSC}
    set -eu
done 

# move all combined tracks to VM_DIR_UCSC
set +eu
mv *combined*.bw ${VM_DIR_UCSC}
set -eu

# Generate tracks for singular and combined bigwig files
./Generate_Tracks.sh

# copy track to ${VM_DIR_UCSC}
cp -a UCSC_Track_Lines ${VM_DIR_UCSC}

echo 'UCSC files for data visualization should now be copied to the VM.'
echo '#--------------------------------------------------------------------------'
##################################################################################
