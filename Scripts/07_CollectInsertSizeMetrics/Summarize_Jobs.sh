#!/bin/bash

##################################################################################
#Andy Rampersaud, 02.16.16
#This script would be used to summarize CollectInsertSizeMetrics statistics 
#Way to run script:
#Usage: 
#./CollectInsertSizeMetrics_Summary.sh
##################################################################################

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

# output directory
output_dir="$(pwd)/Job_Summary"
rm -rf "${output_dir}" && mkdir -p "${output_dir}"

output_file=$output_dir/CollectInsertSizeMetrics_Stats.txt
output_file_pdf=$output_dir/CollectInsertSizeMetrics_Plots.pdf

rm -rf "${output_file}" "${output_file_pdf}"
touch "${output_file}" "${output_file_pdf}"

#Print header to output file:
echo Sample_ID $'\t'Description $'\t'MEDIAN_INSERT_SIZE $'\t'MEDIAN_ABSOLUTE_DEVIATION $'\t'MIN_INSERT_SIZE $'\t'MAX_INSERT_SIZE $'\t'MEAN_INSERT_SIZE $'\t'STANDARD_DEVIATION $'\t'READ_PAIRS $'\t'PAIR_ORIENTATION >> "${output_file}"

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

# loop over all samples
for ((i=0;i< ${#samples[@]} ;i+=3));
do
    sample_dir=${samples[i]}
    sample_id=${samples[i+1]}
    description=${samples[i+2]}

    #Need to cd to sample specific CollectInsertSizeMetrics folder
    pushd "${DATASET_DIR}/${sample_id}/fastq/tophat2/CollectInsertSizeMetrics"

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
    
    popd
done

pushd "${output_dir}"
pdftk *"_hist.pdf" cat output "${output_file_pdf}"

#Remove sample files:
rm *"_hist.pdf"
popd

echo 'Check out '$output_file
echo '#--------------------------------------------------------------------------'

