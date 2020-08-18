#!/bin/bash
##################################################################################
#Andy Rampersaud, 02.23.16
#This script would be used to summarize CollectRnaSeqMetrics statistics 
#Way to run script:
#Usage: 
#./CollectRnaSeqMetrics_Summary.sh
##################################################################################

set -o errexit
set -o pipefail
set -o nounset

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

##################################################################################

output_dir="$(pwd)/Job_Summary"
rm -rf "${output_dir}" && mkdir -p "${output_dir}"

output_file="$output_dir/CollectRnaSeqMetrics_Stats.txt"
rm -rf "${output_file}" && touch "${output_file}"

echo SAMPLE_ID $'\t'DESCRIPTION $'\t'PCT_CODING_BASES $'\t'PCT_UTR_BASES $'\t'PCT_INTRONIC_BASES $'\t'PCT_INTERGENIC_BASES $'\t'PCT_MRNA_BASES $'\t'PCT_RIBOSOMAL_BASES >> $output_file

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

# loop over all samples
for ((i=0;i< ${#samples[@]} ;i+=3));
do
    sample_dir=${samples[i]}
    sample_id=${samples[i+1]}
    description=${samples[i+2]}

    echo "${sample_id}"

    #Need to cd to sample specific CollectRnaSeqMetrics folder
    pushd "${DATASET_DIR}/${sample_id}/tophat2/CollectRnaSeqMetrics"

    #Extract data from temp1.txt

    #Extract lines 7 and 8 from file
    sed -n 7,8p ${sample_id}'_metrics' > temp1.txt
    PCT_CODING_BASES=$(awk '{print $12}' temp1.txt | sed -n 2p)

    #echo $PCT_CODING_BASES
    PCT_UTR_BASES=$(awk '{print $13}' temp1.txt | sed -n 2p)

    #echo $PCT_UTR_BASES
    PCT_INTRONIC_BASES=$(awk '{print $14}' temp1.txt | sed -n 2p)

    #echo $PCT_INTRONIC_BASES
    PCT_INTERGENIC_BASES=$(awk '{print $15}' temp1.txt | sed -n 2p)

    #echo $PCT_INTERGENIC_BASES
    PCT_MRNA_BASES=$(awk '{print $16}' temp1.txt | sed -n 2p)

    #echo $PCT_MRNA_BASES
    PCT_RIBOSOMAL_BASES=$(awk '{print $11}' temp1.txt | sed -n 2p)

    #echo $PCT_RIBOSOMAL_BASES

    #Remove temp1.txt
    rm temp1.txt
    #Print to output file
    echo ${sample_id} $'\t'$description $'\t'$PCT_CODING_BASES $'\t'$PCT_UTR_BASES $'\t'$PCT_INTRONIC_BASES $'\t'$PCT_INTERGENIC_BASES $'\t'$PCT_MRNA_BASES $'\t'$PCT_RIBOSOMAL_BASES >> "${output_file}"
    popd
done

echo 'Check out '$output_file
echo '#-------------------------------------------------------------------------'
##################################################################################
