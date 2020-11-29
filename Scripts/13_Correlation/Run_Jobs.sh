#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

#Clean directory
rm -rf ./logs && mkdir -p ./logs

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

rm -rf Job_Summary && mkdir -p Job_Summary

#dedirs=$(find ../ -maxdepth 1 -name "09*" -type d | xargs -n1  basename |cut -d "_" -f1 | sort | uniq)

# extract all 09* dirs which contain ExonCollapsed subfolders
dedirs=$(find ../09* -iname "output*exoncoll*" | grep -Po '09[a-z]' | sort | uniq)

for d in $dedirs
do
    (set -x; qsub -N "Step_13_$d" -P "${PROJECT}" -l h_rt="04:00:00" correlation.qsub $d)
done

# (set -x; qsub -N "Step_13" -P "${PROJECT}" -l h_rt="04:00:00" correlation.qsub)

echo "End of 13_Correlation"
echo "-----------------------"
