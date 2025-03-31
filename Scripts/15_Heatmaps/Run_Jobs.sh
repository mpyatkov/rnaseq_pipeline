#!/bin/bash

module load R/4.4.3

output_folder=$(pwd)
input_folders=$(dirname "$output_folder")
results="${output_folder}/Job_Summary"

rm -rf "${results}" && mkdir -p "${results}"

# start R Scipts
Rscript ./Scripts/heatmap.R "${input_folders}" "${results}"
rm -rf Rplots.pdf
