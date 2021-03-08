#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

##################################################################################
#Andy Rampersaud, 02.16.16
#This script would be used to summarize CollectInsertSizeMetrics statistics 
#Way to run script:
#Usage: 
#./CollectInsertSizeMetrics_Summary.sh
##################################################################################

# FULL_RECALC=1 by default if nothing provided
FULL_RECALC=${1:-0}

# Skip this step if recalculation flag set to 0
if [ ${FULL_RECALC} -eq 0 ]; then
    echo "Recalculation is not required. FULL_RECALC set to 0."
    exit 0
fi

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

# output directory
output_dir="$(pwd)/Job_Summary"
rm -rf "${output_dir}" && mkdir -p "${output_dir}"

# size metrics
size_metrics="${output_dir}/CollectInsertSizeMetrics_Stats.txt"
size_metrics_pdf="${output_dir}/CollectInsertSizeMetrics_Plots.pdf"
rm -rf "${size_metrics}" "${size_metrics_pdf}"
touch "${size_metrics}"

# header for size metrics 
echo Sample_ID $'\t'Description $'\t'MEDIAN_INSERT_SIZE $'\t'MEDIAN_ABSOLUTE_DEVIATION $'\t'MIN_INSERT_SIZE $'\t'MAX_INSERT_SIZE $'\t'MEAN_INSERT_SIZE $'\t'STANDARD_DEVIATION $'\t'READ_PAIRS $'\t'PAIR_ORIENTATION >> "${size_metrics}"

# --------------------

# rnaseq metrics
rna_metrics="${output_dir}/CollectRnaSeqMetrics_Stats.txt"
rm -rf "${rna_metrics}" && touch "${rna_metrics}"

# header for rnaseq metrics 
echo SAMPLE_ID $'\t'DESCRIPTION $'\t'CODING_BASES $'\t'UTR_BASES $'\t'PCT_INTRONIC_BASES $'\t'PCT_INTERGENIC_BASES $'\t'PCT_MRNA_BASES $'\t'RIBOSOMAL_BASES >> "${rna_metrics}"

collect_size_metrics() {
    local sample_id=$1
    local output_dir=$2
    local output_file=$3

    #Copy hist file to output_dir
    cp ${sample_id}'_hist.pdf' ${output_dir}

    #Extract data from temp1.txt
    #Extract lines 7 and 8 from file
    sed -n 7,8p $sample_id'_metrics' > temp1.txt
    MEDIAN_INSERT_SIZE=$(awk '{print $1}' temp1.txt | sed -n 2p)
    
    #echo $MEDIAN_INSERT_SIZE
    MEDIAN_ABSOLUTE_DEVIATION=$(awk '{print $2}' temp1.txt | sed -n 2p)
    
    #echo $MEDIAN_ABSOLUTE_DEVIATION
    MIN_INSERT_SIZE=$(awk '{print $3}' temp1.txt | sed -n 2p)
    
    #echo $MIN_INSERT_SIZE
    MAX_INSERT_SIZE=$(awk '{print $4}' temp1.txt | sed -n 2p)
    
    #echo $MAX_INSERT_SIZE
    MEAN_INSERT_SIZE=$(awk '{print $5}' temp1.txt | sed -n 2p)
    
    #echo $MEAN_INSERT_SIZE
    STANDARD_DEVIATION=$(awk '{print $6}' temp1.txt | sed -n 2p)
    
    #echo $STANDARD_DEVIATION
    READ_PAIRS=$(awk '{print $7}' temp1.txt | sed -n 2p)
    
    #echo $READ_PAIRS
    PAIR_ORIENTATION=$(awk '{print $8}' temp1.txt | sed -n 2p)
    
    #echo $PAIR_ORIENTATION
    #Remove temp1.txt
    rm temp1.txt
    
    #Print to output file
    echo $sample_id $'\t'$description $'\t'$MEDIAN_INSERT_SIZE $'\t'$MEDIAN_ABSOLUTE_DEVIATION $'\t'$MIN_INSERT_SIZE $'\t'$MAX_INSERT_SIZE $'\t'$MEAN_INSERT_SIZE $'\t'$STANDARD_DEVIATION $'\t'$READ_PAIRS $'\t'$PAIR_ORIENTATION >> "${output_file}"

}

collect_rnaseq_metrics() {

    local sample_id=$1
    local output_file=$2

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

}

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

# loop over all samples
for ((i=0;i< ${#samples[@]} ;i+=3));
do
    sample_dir=${samples[i]}
    sample_id=${samples[i+1]}
    description=${samples[i+2]}

    # Need to cd to sample specific CollectInsertSizeMetrics folder
    pushd "${DATASET_DIR}/${sample_id}/CollectInsertSizeMetrics"
    collect_size_metrics ${sample_id} ${output_dir} ${size_metrics}
    popd

    # Need to cd to sample specific CollectRnaSeqMetrics folder
    pushd "${DATASET_DIR}/${sample_id}/CollectRnaSeqMetrics"
    echo "---> $(pwd)"
    ls -l
    collect_rnaseq_metrics ${sample_id} ${rna_metrics}
    popd
done

pushd "${output_dir}"
pdftk *"_hist.pdf" cat output "${size_metrics_pdf}"

#Remove sample files:
rm *"_hist.pdf"
popd


