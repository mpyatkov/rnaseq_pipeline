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
# INPUT_DIR=$(pwd)
# cd $INPUT_DIR

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

# loop over all samples
for ((i=0;i< ${#samples[@]} ;i+=3));
do
    Sample_DIR=${samples[i]}
    Sample_ID=${samples[i+1]}
    Description=${samples[i+2]}
    #Copy over BAM and *.bw files:
    echo 'Copying over BAM, bigWig, and bigBed files to waxmanlabvm...'
    set +eu
    cp ${DATASET_DIR}/${Sample_ID}/tophat2/$Sample_ID'_primary_unique.bam' ${VM_DIR_UCSC}
    cp ${DATASET_DIR}/${Sample_ID}/tophat2/$Sample_ID'_primary_unique.bam.bai' ${VM_DIR_UCSC}
    cp ${DATASET_DIR}/${Sample_ID}/tophat2/UCSC_BigWig/*'.bw' ${VM_DIR_UCSC}
    set -eu
done 
# cd $INPUT_DIR

./Generate_Tracks.sh

echo 'UCSC files for data visualization should now be copied to the VM.'
echo '#--------------------------------------------------------------------------'
##################################################################################
