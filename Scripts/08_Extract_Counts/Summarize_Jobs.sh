#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

# FULL_RECALC=1 by default if nothing provided
FULL_RECALC=${1:-0}

# Skip this step if recalculation flag set to 0
if [ ${FULL_RECALC} -eq 0 ]; then
    echo "Recalculation is not required. FULL_RECALC set to 0."
    exit 0
fi

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

if [ ${DEFAULT_ALIGNER} -eq 1 ]; then
    echo "Summarizing for STAR aligner is not working at the moment. STAR output reports are different from TopHat output reports."
    exit 0
fi

# output directory
summary_dir="$(pwd)/Job_Summary"
rm -rf "${summary_dir}" && mkdir -p "${summary_dir}"

# gtflist
gtflist=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --gtf_annotation_and_counter))

# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

for ((i=0;i< ${#gtflist[@]} ;i+=2));
do
    gtf_file=${gtflist[i]}
    counter=${gtflist[i+1]}

    eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export_gtf_by_name_and_counter ${gtf_file} ${counter})"

    OUTPUT_FILE=${summary_dir}/${COUNTER}_Stats_${OUTPUT_DIR}.txt
    rm -rf ${OUTPUT_FILE} && touch ${OUTPUT_FILE}
    
    if [[ $COUNTER == "featureCounts" ]]
    then
	OUTPUT_FILE_SPECIAL="${summary_dir}/${COUNTER}_summary_${OUTPUT_DIR}.txt"
	rm -rf ${OUTPUT_FILE} && touch ${OUTPUT_FILE}
    fi

    ## REM_SPLICE used only because we need only intronic gtf files
    if [[ $REM_SPLICE_JUNC != 0 ]]
    then
        echo 'SAMPLE_ID' $'\t''DESCRIPTION' $'\t''Aligned_Pairs' $'\t''Read_Pair_IN_INTRONS' $'\t''Read_Pair_IN_INTRONS_RATIO' >> $OUTPUT_FILE
    else
        echo 'SAMPLE_ID' $'\t''DESCRIPTION' $'\t''Aligned_Pairs' $'\t''Read_Pair_IN_EXONS' $'\t''Read_Pair_IN_EXONS_RATIO' >> $OUTPUT_FILE
    fi
    
    for ((j=0;j<${#samples[@]} ;j+=3));
    do
        sample_dir=${samples[j]}
        sample_id=${samples[j+1]}
        description=${samples[j+2]}

        pushd "${DATASET_DIR}/${sample_id}/aligner"

        Aligned_Pairs=$(grep "Aligned pairs:"  "${sample_id}_align_summary.txt"  | awk '{print $3}')
        counter_DIR="${DATASET_DIR}/${sample_id}/${COUNTER}"
            
        pushd "${counter_DIR}/${OUTPUT_DIR}"
	

	# READS_IN_EXONS=$(awk '{n+=$2;} ; END {print n;}' "${sample_id}_${COUNTER}.out")
	if [[ $COUNTER == "htseq" ]]
	then
	    # need to remove last 5 lines before processing
	    READS_IN_EXONS=$(cat "${sample_id}_${COUNTER}.out" | head -n-5 | awk '{n+=$2;} ; END {print n;}')
        else
	    READS_IN_EXONS=$(awk '{n+=$2;} ; END {print n;}' "${sample_id}_${COUNTER}.out")
	fi

        # Calculate percentages
        READS_IN_EXONS_RATIO=$(echo "scale=4;$READS_IN_EXONS/$Aligned_Pairs" | bc)
        echo $sample_id $'\t'$description $'\t'$Aligned_Pairs $'\t'$READS_IN_EXONS $'\t'$READS_IN_EXONS_RATIO >> $OUTPUT_FILE
	popd
	popd
    done

    if [[ $COUNTER == "htseq" ]]
    then
	echo 'SAMPLE_ID' $'\t''DESCRIPTION' $'\t''Aligned_Pairs' $'\t''no_feature_count' $'\t''no_feature_RATIO' $'\t''ambiguous_count'$'\t''ambiguous_RATIO' $'\t''too_low_aQual_count'$'\t''too_low_aQual_RATIO' $'\t''not_aligned_count'$'\t''not_aligned_RATIO' $'\t''alignment_not_unique'$'\t''alignment_not_unique_RATIO' >> ${OUTPUT_FILE}
    fi

    for ((j=0; j<${#samples[@]} ;j+=3));
    do
        sample_dir=${samples[j]}
        sample_id=${samples[j+1]}
        description=${samples[j+2]}

        pushd "${DATASET_DIR}/${sample_id}/aligner"

        Aligned_Pairs=$(grep "Aligned pairs:"  "${sample_id}_align_summary.txt"  | awk '{print $3}')
        counter_DIR="${DATASET_DIR}/${sample_id}/${COUNTER}"
	pushd "${counter_DIR}/${OUTPUT_DIR}"

	if [[ $COUNTER == "htseq" ]]
	then

	    echo 'Getting no_feature count...'
            no_feature=$(grep '__no_feature' "${sample_id}_${COUNTER}.out"  | awk '{print $2}')
            
            echo 'Calculate percentage...'
            no_feature_RATIO=$(echo "scale=4;$no_feature/$Aligned_Pairs" | bc)
            
            echo 'Getting ambiguous count...'
            ambiguous=$(grep '__ambiguous' "${sample_id}_${COUNTER}.out"  | awk '{print $2}')
            
            echo 'Calculate percentage...'
            ambiguous_RATIO=$(echo "scale=4;$ambiguous/$Aligned_Pairs" | bc)
            
            echo 'Getting too_low_aQual count...'
            too_low_aQual=$(grep '__too_low_aQual' "${sample_id}_${COUNTER}.out"  | awk '{print $2}')
            echo 'Calculate percentage...'
            too_low_aQual_RATIO=$(echo "scale=4;$too_low_aQual/$Aligned_Pairs" | bc)
            
            echo 'Getting not_aligned count...'
            not_aligned=$(grep '__not_aligned' "${sample_id}_${COUNTER}.out"  | awk '{print $2}')
            echo 'Calculate percentage...'
            not_aligned_RATIO=$(echo "scale=4;$not_aligned/$Aligned_Pairs" | bc)
            
            echo 'Getting alignment_not_unique count...'
            alignment_not_unique=$(grep '__alignment_not_unique' "${sample_id}_${COUNTER}.out"  | awk '{print $2}')
            
            echo 'Calculate percentage...'
            alignment_not_unique_RATIO=$(echo "scale=4;$alignment_not_unique/$Aligned_Pairs" | bc)
            echo 'Printing to output file...'
            echo $sample_id $'\t'$description $'\t'$Aligned_Pairs $'\t'$no_feature $'\t'$no_feature_RATIO $'\t'$ambiguous $'\t'$ambiguous_RATIO $'\t'$too_low_aQual $'\t'$too_low_aQual_RATIO $'\t'$not_aligned $'\t'$not_aligned_RATIO $'\t'$alignment_not_unique $'\t'$alignment_not_unique_RATIO >> ${OUTPUT_FILE}

	else
	    
	    awk 'NR>1' "${sample_id}_${COUNTER}.out.summary" > temp1.txt
	                
            #Divide 2nd column by total number of reads:
            awk -v Aligned_Pairs=${Aligned_Pairs} '{ Ratio=($2)/(Aligned_Pairs) ; print $1"\t"$2"\t"Aligned_Pairs"\t"Ratio }' temp1.txt > temp2.txt
            #Need a header:
            echo 'Status' $'\t'$sample_id'_Read_Count' $'\t''Aligned_Pairs' $'\t''Ratio' >> Header.txt
            #Combine files:
            cat Header.txt temp2.txt > temp3.txt

            # if [ "${ANNOTATION_FILE}" == "RefSeq_GeneBody.gtf" ];
            # then
            #     sed 's/.out.*/:/' temp3.txt > temp4.txt
            #     mv temp4.txt temp3.txt
            # fi

            echo 'Concatenating summary file...'
            cat temp3.txt >> $OUTPUT_FILE_SPECIAL

            #Remove temp files:
            rm temp*.txt
            rm Header.txt

	fi
	popd
	popd
    done
done
