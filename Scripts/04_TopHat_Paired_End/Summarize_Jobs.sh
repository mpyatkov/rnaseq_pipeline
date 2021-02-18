#!/bin/bash
##################################################################################
# Andy Rampersaud, 02.23.16
#This script would be used to summarize TopHat2 mapping statistics 
#Way to run script:
#Usage: 
#./TopHat_Paired_End_Summary.sh
##################################################################################

set -o errexit
set -o pipefail
set -o nounset

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

if [ ${DEFAULT_ALIGNER} -eq 1 ]; then
    echo "Summarizing for STAR aligner in progress..."
    exit 0
fi

output_dir="$(pwd)/Job_Summary"
rm -rf "${output_dir}" && mkdir -p "${output_dir}"

output_file="${output_dir}/TopHat2_Stats.txt"
rm -rf "${output_file}" && touch "${output_file}"

output_file_bestmapped="${output_dir}/TopHat2_Stats_BestMapped.txt"
rm -rf "${output_file_bestmapped}" && touch "${output_file_bestmapped}"

output_file_uniquereads="${output_dir}/TopHat2_Stats_UniqueReads.txt"
rm -rf "${output_file_uniquereads}" && touch "${output_file_uniquereads}"

output_file_splicereads="${output_dir}/TopHat2_Stats_SpliceReads.txt"
rm -rf "${output_file_splicereads}" && touch "${output_file_splicereads}"

#-----------------------------------------------
# Minor issue with jobs using multiple cores
# Empty job output log files (*.pe* and *.po*) are created
# Remove them if they exist:
rm -rf *.pe* *.po*

# TODO: refactor this later
# Print header to output file:
echo 'Read_Mate' $'\t''SAMPLE_ID' $'\t''DESCRIPTION' $'\t''TOTAL_SEQUENCED_READS' $'\t''MULTIMAPPED_READS_ONLY' $'\t''MAPPED_READS_(BEST_LOCATION_MAPPED)' $'\t''BEST_LOCATION_MAPPED_RATIO' >> $output_file

# Print header to output file:
echo 'SAMPLE_ID' $'\t''DESCRIPTION' $'\t''Aligned_Pairs' $'\t''Concordant_pair_alignment_rate' $'\t''TOTAL_SEQUENCED_READS' $'\t''MAPPED_READS_(BEST_LOCATION_MAPPED)' $'\t''MAPPED_READS_(BEST_LOCATION_MAPPED)_RATIO' >> $output_file_bestmapped

# Print header to output file:
echo 'SAMPLE_ID' $'\t''DESCRIPTION' $'\t''Aligned_Pairs' $'\t''Concordant_pair_alignment_rate' $'\t''TOTAL_SEQUENCED_READS' $'\t''uniquely_mapped_reads' $'\t''uniquely_mapped_reads_RATIO' >> $output_file_uniquereads

# Print header to output file:
echo 'SAMPLE_ID' $'\t''DESCRIPTION' $'\t''Aligned_Pairs' $'\t''Concordant_pair_alignment_rate' $'\t''TOTAL_SEQUENCED_READS' $'\t''splice_junction_READS' $'\t''SPLICED_READ_RATIO' >> $output_file_splicereads


# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

# rename directories from sample_dir to sample_id
for ((i=0;i< ${#samples[@]} ;i+=3));
do
    sample_dir=${samples[i]}
    sample_id=${samples[i+1]}
    description=${samples[i+2]}

    pushd "${DATASET_DIR}/${sample_id}/aligner"

    #Summary for $output_file
    # TODO: refactor this latter

    #Get counts from *_align_summary.txt file
    TOTAL_SEQUENCED_READS_Left=$(grep 'Input' $sample_id'_align_summary.txt'  | awk '{print $3}' | awk 'NR==1 {print $1}')
    TOTAL_SEQUENCED_READS_Right=$(grep 'Input' $sample_id'_align_summary.txt'  | awk '{print $3}' | awk 'NR==2 {print $1}')

    #Get counts from *_align_summary.txt file
    MULTIMAPPED_READS_ONLY_Left=$(grep 'of these:' $sample_id'_align_summary.txt'  | awk '{print $3}' | awk 'NR==1 {print $1}')
    MULTIMAPPED_READS_ONLY_Right=$(grep 'of these:' $sample_id'_align_summary.txt'  | awk '{print $3}' | awk 'NR==2 {print $1}')

    #Get counts from *_align_summary.txt file
    MAPPED_READS_BEST_LOCATION_MAPPED_Left=$(grep 'Mapped'  $sample_id'_align_summary.txt'  | awk '{print $3}' | awk 'NR==1 {print $1}')
    MAPPED_READS_BEST_LOCATION_MAPPED_Right=$(grep 'Mapped'  $sample_id'_align_summary.txt'  | awk '{print $3}' | awk 'NR==2 {print $1}')

    #Calculate percentages
    BEST_LOCATION_MAPPED_RATIO_Left=$(echo "scale=4;$MAPPED_READS_BEST_LOCATION_MAPPED_Left/$TOTAL_SEQUENCED_READS_Left" | bc)
    BEST_LOCATION_MAPPED_RATIO_Right=$(echo "scale=4;$MAPPED_READS_BEST_LOCATION_MAPPED_Right/$TOTAL_SEQUENCED_READS_Right" | bc)
    #Print to output file
    echo 'Left reads:' $'\t'$sample_id $'\t'$description $'\t'$TOTAL_SEQUENCED_READS_Left $'\t'$MULTIMAPPED_READS_ONLY_Left $'\t'$MAPPED_READS_BEST_LOCATION_MAPPED_Left $'\t'$BEST_LOCATION_MAPPED_RATIO_Left >> $output_file
    echo 'Right reads:' $'\t'$sample_id $'\t'$description $'\t'$TOTAL_SEQUENCED_READS_Right $'\t'$MULTIMAPPED_READS_ONLY_Right $'\t'$MAPPED_READS_BEST_LOCATION_MAPPED_Right $'\t'$BEST_LOCATION_MAPPED_RATIO_Right >> $output_file
    
    #Summary for $OUTPUT_FILE_BestMapped
    
    #Need to get the number of Aligned pairs:
    Aligned_Pairs=$(grep 'Aligned pairs:'  $sample_id'_align_summary.txt'  | awk '{print $3}')

    #Get the "concordant pair alignment rate"
    Align_Rate=$(grep 'concordant pair alignment rate.'  $sample_id'_align_summary.txt'  | awk '{print $1}')
    
    #Need the TOTAL_SEQUENCED_READS = (left reads + right reads)
    TOTAL_SEQUENCED_READS=$(echo "scale=4;$TOTAL_SEQUENCED_READS_Left + $TOTAL_SEQUENCED_READS_Right" | bc)

    #Need the *_statistics_for_primary_reads.txt
    BEST_Mapped_Reads=$(grep 'total' $sample_id'_statistics_for_primary_reads.txt'  | awk '{print $1}')

    #Calculate percentages
    BEST_Mapped_Reads_RATIO=$(echo "scale=4;$BEST_Mapped_Reads/$TOTAL_SEQUENCED_READS" | bc)

    #Print to output file
    echo $sample_id $'\t'$description $'\t'$Aligned_Pairs $'\t'$Align_Rate $'\t'$TOTAL_SEQUENCED_READS $'\t'$BEST_Mapped_Reads $'\t'$BEST_Mapped_Reads_RATIO >> $output_file_bestmapped

    #Summary for $OUTPUT_FILE_UniqueReads
    #Need to summarize the uniquely mapped reads:
    uniquely_mapped_reads=$(grep 'total' $sample_id'_statistics_for_primary_unique_reads.txt'  | awk '{print $1}')

    #Calculate percentage:
    uniquely_mapped_reads_RATIO=$(echo "scale=4;$uniquely_mapped_reads/$TOTAL_SEQUENCED_READS" | bc)

    #Print to output file
    echo $sample_id $'\t'$description $'\t'$Aligned_Pairs $'\t'$Align_Rate $'\t'$TOTAL_SEQUENCED_READS  $'\t'$uniquely_mapped_reads $'\t'$uniquely_mapped_reads_RATIO >> $output_file_uniquereads

    #Summary for $OUTPUT_FILE_SpliceReads
    
    #Print 2nd line of the *_spliced_read_counts.txt:
    #Reminder: the following counts are from the single BAM file
    #SKIPPED_READS=$(awk 'NR == 2' $sample_id'_spliced_read_counts.txt' | awk '{print $2}')
    SPLICED_READS=$(awk 'NR == 2' $sample_id'_spliced_read_counts.txt' | awk '{print $3}')

    #----------------------------------------------
    #Calculate percentage:
    SPLICED_READS_RATIO=$(echo "scale=4;$SPLICED_READS/$TOTAL_SEQUENCED_READS" | bc)
    #Print to output file
    echo $sample_id $'\t'$description $'\t'$Aligned_Pairs $'\t'$Align_Rate $'\t'$TOTAL_SEQUENCED_READS  $'\t'$SPLICED_READS $'\t'$SPLICED_READS_RATIO >> $output_file_splicereads
    #----------------------------------------------
    popd
done

echo "Done Summarise TOPHAT"

