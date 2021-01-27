#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

#Clean directory
rm -rf ./logs && mkdir -p ./logs

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

rm -rf Job_Summary && mkdir -p Job_Summary

# create table (SAMPLE_ID \t Condition_Name \t Condition_Name_SAMPLE_ID)
SAMPLE_FILE="SAMPLE_CONDNAME.tsv"
rm -rf ${SAMPLE_FILE}

groups=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --groups))
# for each group
for group in "${groups[@]}"; do
    # get samples array (triples of sample_dir sample_id description for each sample)
    samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples_by_group ${group}))
    for ((i=0;i< ${#samples[@]} ;i+=3));
    do
	SAMPLE_ID=${samples[i+1]}
	CONDNAME=${samples[i+2]}
	printf "%s\t%s\t%s\n" "${SAMPLE_ID}" "${CONDNAME}" "${CONDNAME}_${SAMPLE_ID}" >> ${SAMPLE_FILE}
    done
done


# extract all 09* dirs which contain ExonCollapsed subfolders
dedirs=$(find ../09* -iname "output*exoncoll*" | grep -Po '09[a-z]' | sort | uniq)

for d in $dedirs
do
    (set -x; qsub -N "Step_13_$d" -P "${PROJECT}" -l h_rt="04:00:00" correlation.qsub $d)
done

# (set -x; qsub -N "Step_13" -P "${PROJECT}" -l h_rt="04:00:00" correlation.qsub)

echo "End of 13_Correlation"
echo "-----------------------"
