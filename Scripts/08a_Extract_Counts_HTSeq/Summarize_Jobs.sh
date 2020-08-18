#!/bin/bash

##################################################################################
#Andy Rampersaud, 03.11.16
#This script would be used to summarize Extract_Counts statistics 
#Way to run script:
#Usage: 
#./Extract_Counts_Summary.sh
##################################################################################

set -o errexit
set -o pipefail
set -o nounset

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

rm -rf *.po* *.pe*

# output directory
output_dir="$(pwd)/Job_Summary"
rm -rf "${output_dir}" && mkdir -p "${output_dir}"

#Start loop over GTF files:
gtf_list=${GTF_FILES_DIR}/*.gtf
for gtf_file in ${gtf_list}
do
    ANNOTATION_FILE=$(basename ${gtf_file})

    #Skip all of these steps if we are dealing with lncRNA gtf files
    if [ "${ANNOTATION_FILE}" == "ncRNA_exon_for_counting.gtf" -o "${ANNOTATION_FILE}" == "ncRNA_genebodies_for_counting.gtf" -o "${ANNOTATION_FILE}" == "intronic_only_gene_models_ncRNA_for_counting.gtf" -o "${ANNOTATION_FILE}" == "exonic_only_gene_models_ncRNA_for_counting.gtf" ];
    then
        echo "lncRNA GTF file is ignored by HTSeq, see FeatureCount result output"
    else
        # Need if statements to determine output_file depending on the
        # ANNOTATION_FILE used:
        
        if [ "${ANNOTATION_FILE}" == "genes.gtf" ];
        then
            output_file=${output_dir}/HTSeq_Stats_Illumina_GTF.txt
        fi
        #------------------------------------------
        if [ "${ANNOTATION_FILE}" == "RefSeq_GeneBody.gtf" ];
        then
            output_file=${output_dir}/HTSeq_Stats_RefSeq_Exon_GTF.txt
        fi
        #------------------------------------------
        if [ "${ANNOTATION_FILE}" == "Intron_Only_Regions.gtf" ];
        then
            output_file=${output_dir}/HTSeq_Stats_RefSeq_Intron_GTF.txt
        fi
        #------------------------------------------
        if [ "${ANNOTATION_FILE}" == "Exon_Only_Regions.gtf" ];
        then
            output_file=${output_dir}/HTSeq_Stats_RefSeq_Exon_Only_GTF.txt
        fi

        rm -rf ${output_file} && touch ${output_file}

        #Print header to output file:
        #------------------------------------------
        if [ "${ANNOTATION_FILE}" == "Intron_Only_Regions.gtf" ];
        then
            echo 'SAMPLE_ID' $'\t''DESCRIPTION' $'\t''Aligned_Pairs' $'\t''Read_Pair_IN_INTRONS' $'\t''Read_Pair_IN_INTRONS_RATIO' >> $output_file
        fi
        #------------------------------------------
        if [ "${ANNOTATION_FILE}" == "genes.gtf" ] || [ "${ANNOTATION_FILE}" == "RefSeq_GeneBody.gtf" ] || [ "${ANNOTATION_FILE}" == "Exon_Only_Regions.gtf" ];
        then
            echo 'SAMPLE_ID' $'\t''DESCRIPTION' $'\t''Aligned_Pairs' $'\t''Read_Pair_IN_EXONS' $'\t''Read_Pair_IN_EXONS_RATIO' >> $output_file
        fi


        # samples contains array of (sample_dir, sample_id, description) for each sample
        samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

        # loop over all samples
        for ((i=0;i< ${#samples[@]} ;i+=3));
        do
            sample_dir=${samples[i]}
            sample_id=${samples[i+1]}
            description=${samples[i+2]}
            
            echo $sample_id
            #Need to cd to sample specific tophat2 folder
            cd ${DATASET_DIR}/$sample_id/tophat2
            #------------------------------------------
            #Get counts from *_align_summary.txt file
            # echo 'Getting the total number of mapped reads...'
            # #(statistics_for_all_accepted_reads.txt)
            # TOTAL_MAPPED_READS=$(grep 'total' statistics_for_all_accepted_reads.txt  | awk '{print $1}')
            #------------------------------------------
            #http://www-huber.embl.de/users/anders/HTSeq/doc/count.html
            # For paired-end data, does htseq-count count reads or read pairs?
            # Read pairs. The script is designed to count “units of evidence” for gene expression. If both mates map to the same gene, this still only shows that one cDNA fragment originated from that gene. Hence, it should be counted only once.
            #------------------------------------------
            #----------------------------------------------
            echo 'Getting the Aligned_Pairs...'
            #Need to get the number of Aligned pairs:
            Aligned_Pairs=$(grep 'Aligned pairs:'  $sample_id'_align_summary.txt'  | awk '{print $3}')
            #------------------------------------------
            echo 'Getting READS_IN_EXONS...'
            #Need to cd to sample specific HTSeq folder
            #Need if statements to process the *_HTSeq.out file depending on the ANNOTATION_FILE used:
            HTSeq_DIR=${DATASET_DIR}/$sample_id/tophat2/HTSeq
            #------------------------------------------
            if [ "${ANNOTATION_FILE}" == "genes.gtf" ];
            then
                cd ${HTSeq_DIR}/Illumina_GTF
            fi
            #------------------------------------------
            if [ "${ANNOTATION_FILE}" == "RefSeq_GeneBody.gtf" ];
            then
                cd ${HTSeq_DIR}/RefSeq_Exon_GTF
            fi
            #------------------------------------------
            if [ "${ANNOTATION_FILE}" == "Intron_Only_Regions.gtf" ];
            then
                cd ${HTSeq_DIR}/RefSeq_Intron_GTF
            fi
            #------------------------------------------
            if [ "${ANNOTATION_FILE}" == "Exon_Only_Regions.gtf" ];
            then
                cd ${HTSeq_DIR}/RefSeq_Exon_Only_GTF
            fi
            #------------------------------------------
            #Need to remove the last 5 lines
            head -n -5 $sample_id'_HTSeq.out' > $sample_id'_HTSeq.temp1'
            #Get the sum of the 2nd coulumn
            READS_IN_EXONS=$(awk '{n+=$2;} ; END {print n;}' $sample_id'_HTSeq.temp1')
            #Remove temp* files
            rm $sample_id'_HTSeq.temp1'
            echo 'Getting READS_IN_EXONS_RATIO...'
            #Calculate percentages
            READS_IN_EXONS_RATIO=$(echo "scale=4;$READS_IN_EXONS/$Aligned_Pairs" | bc)
            echo 'Printing to output file...'
            echo $sample_id $'\t'$description $'\t'$Aligned_Pairs $'\t'$READS_IN_EXONS $'\t'$READS_IN_EXONS_RATIO >> $output_file
            #End loop over Sample_Labels:
        done 
        ##################################################################################
        #Need to summarize HTSeq special counters
        ##################################################################################
                
        SCRIPT_DIR=$(pwd)
        
        echo 'SAMPLE_ID' $'\t''DESCRIPTION' $'\t''Aligned_Pairs' $'\t''no_feature_count' $'\t''no_feature_RATIO' $'\t''ambiguous_count'$'\t''ambiguous_RATIO' $'\t''too_low_aQual_count'$'\t''too_low_aQual_RATIO' $'\t''not_aligned_count'$'\t''not_aligned_RATIO' $'\t''alignment_not_unique'$'\t''alignment_not_unique_RATIO' >> $output_file
        #------------------------------------------
        # loop over all samples
        for ((i=0;i< ${#samples[@]} ;i+=3));
        do
            sample_dir=${samples[i]}
            sample_id=${samples[i+1]}
            description=${samples[i+2]}
        
            echo $sample_id
            #Need to cd to sample specific tophat2 folder
            cd ${DATASET_DIR}/$sample_id/tophat2
            
            
            #Get counts from *_align_summary.txt file
            # echo 'Getting the total number of mapped reads...'
            # #(statistics_for_all_accepted_reads.txt)
            # TOTAL_MAPPED_READS=$(grep 'total' \
            # statistics_for_all_accepted_reads.txt | awk '{print $1}')
            #------------------------------------------
            # http://www-huber.embl.de/users/anders/HTSeq/doc/count.html
            # For paired-end data, does htseq-count count reads or read pairs?
            # Read pairs. The script is designed to count “units of evidence”
            # for gene expression. If both mates map to the same gene, this
            # still only shows that one cDNA fragment originated from that
            # gene. Hence, it should be counted only once.
            #----------------------------------------------
            echo 'Getting the Aligned_Pairs...'

            #Need to get the number of Aligned pairs:
            Aligned_Pairs=$(grep "Aligned pairs:"  "${sample_id}_align_summary.txt"  | awk '{print $3}')
            
            # Need to cd to sample specific HTSeq folder
            # Need if statements to process the *_HTSeq.out file depending on the
            # ANNOTATION_FILE used:

            HTSeq_DIR=${DATASET_DIR}/$sample_id/tophat2/HTSeq
            #------------------------------------------
            
            if [ "${ANNOTATION_FILE}" == "genes.gtf" ];
            then
                annotation_dir=${HTSeq_DIR}/Illumina_GTF
            fi

            if [ "${ANNOTATION_FILE}" == "RefSeq_GeneBody.gtf" ];
            then
                annotation_dir=${HTSeq_DIR}/RefSeq_Exon_GTF
            fi
            
            if [ "${ANNOTATION_FILE}" == "Intron_Only_Regions.gtf" ];
            then
                annotation_dir=${HTSeq_DIR}/RefSeq_Intron_GTF
            fi
            
            if [ "${ANNOTATION_FILE}" == "Exon_Only_Regions.gtf" ];
            then
                annotation_dir=${HTSeq_DIR}/RefSeq_Exon_Only_GTF
            fi

            # cd to annotation_dir
            pushd ${annotation_dir}
            
            echo 'Getting no_feature count...'
            no_feature=$(grep '__no_feature' $sample_id'_HTSeq.out'  | awk '{print $2}')
            
            echo 'Calculate percentage...'
            no_feature_RATIO=$(echo "scale=4;$no_feature/$Aligned_Pairs" | bc)
            
            echo 'Getting ambiguous count...'
            ambiguous=$(grep '__ambiguous' $sample_id'_HTSeq.out'  | awk '{print $2}')
            
            echo 'Calculate percentage...'
            ambiguous_RATIO=$(echo "scale=4;$ambiguous/$Aligned_Pairs" | bc)
            
            echo 'Getting too_low_aQual count...'
            too_low_aQual=$(grep '__too_low_aQual' $sample_id'_HTSeq.out'  | awk '{print $2}')
            echo 'Calculate percentage...'
            too_low_aQual_RATIO=$(echo "scale=4;$too_low_aQual/$Aligned_Pairs" | bc)
            
            echo 'Getting not_aligned count...'
            not_aligned=$(grep '__not_aligned' $sample_id'_HTSeq.out'  | awk '{print $2}')
            echo 'Calculate percentage...'
            not_aligned_RATIO=$(echo "scale=4;$not_aligned/$Aligned_Pairs" | bc)
            
            echo 'Getting alignment_not_unique count...'
            alignment_not_unique=$(grep '__alignment_not_unique' $sample_id'_HTSeq.out'  | awk '{print $2}')
            
            echo 'Calculate percentage...'
            alignment_not_unique_RATIO=$(echo "scale=4;$alignment_not_unique/$Aligned_Pairs" | bc)
            echo 'Printing to output file...'
            echo $sample_id $'\t'$description $'\t'$Aligned_Pairs $'\t'$no_feature $'\t'$no_feature_RATIO $'\t'$ambiguous $'\t'$ambiguous_RATIO $'\t'$too_low_aQual $'\t'$too_low_aQual_RATIO $'\t'$not_aligned $'\t'$not_aligned_RATIO $'\t'$alignment_not_unique $'\t'$alignment_not_unique_RATIO >> $output_file
            #End loop over Sample_Labels:

            #go back
            popd
        done 
        echo 'Check out '$output_file
        echo "--------------------"
    fi #End of if statement to ignore lncRNA gtf files
    #End loop over GTF files:
done

echo "Check out: ${output_dir}"
echo "--------------------"

