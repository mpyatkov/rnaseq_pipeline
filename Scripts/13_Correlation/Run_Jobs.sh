#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

#Clean directory
rm -rf ./logs && mkdir -p ./logs

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

(set -x; qsub -N "Step_13" -P "${PROJECT}" -l h_rt="04:00:00" correlation.qsub)

echo "End of 13_Correlation"
echo "-----------------------"
