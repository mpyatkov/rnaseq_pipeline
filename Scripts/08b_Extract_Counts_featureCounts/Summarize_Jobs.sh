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
OUTPUT_DIR="$(pwd)/Job_Summary"
rm -rf "${OUTPUT_DIR}" && mkdir -p "${OUTPUT_DIR}"

#Start loop over GTF files:

gtf_list=${GTF_FILES_DIR}/*.gtf
for gtf_file in ${gtf_list}
do
    ANNOTATION_FILE=$(basename ${gtf_file})

    echo '#------------------------------------------'
    echo 'ANNOTATION_FILE:'
    echo ${ANNOTATION_FILE}
    echo '#------------------------------------------'

    #Need if statements to determine OUTPUT_FILE depending on the ANNOTATION_FILE used:
    #------------------------------------------
    if [ "${ANNOTATION_FILE}" == "genes.gtf" ];
    then
        OUTPUT_FILE=$OUTPUT_DIR/featureCounts_Stats_Illumina_GTF.txt
    fi
    
    if [ "${ANNOTATION_FILE}" == "RefSeq_GeneBody.gtf" ];
    then
        OUTPUT_FILE=$OUTPUT_DIR/featureCounts_Stats_RefSeq_Exon_GTF.txt
    fi
    
    if [ "${ANNOTATION_FILE}" == "Intron_Only_Regions.gtf" ];
    then
        OUTPUT_FILE=$OUTPUT_DIR/featureCounts_Stats_RefSeq_Intron_GTF.txt
    fi
    
    if [ "${ANNOTATION_FILE}" == "Exon_Only_Regions.gtf" ];
    then
        OUTPUT_FILE=$OUTPUT_DIR/featureCounts_Stats_RefSeq_Exon_Only_GTF.txt
    fi

    #Start of lncRNA gtf
    if [ "${ANNOTATION_FILE}" == "ncRNA_exon_for_counting.gtf" ];
    then
        OUTPUT_FILE=$OUTPUT_DIR/featureCounts_Stats_LncRNA_Exon_Collapsed_GTF.txt
    fi
    
    if [ "${ANNOTATION_FILE}" == "ncRNA_genebodies_for_counting.gtf" ];
    then
        OUTPUT_FILE=$OUTPUT_DIR/featureCounts_Stats_LncRNA_GeneBody_GTF.txt
    fi
    
    if [ "${ANNOTATION_FILE}" == "intronic_only_gene_models_ncRNA_for_counting.gtf" ];
    then
        OUTPUT_FILE=$OUTPUT_DIR/featureCounts_Stats_LncRNA_Exonic_Only_GTF.txt
    fi
    
    if [ "${ANNOTATION_FILE}" == "exonic_only_gene_models_ncRNA_for_counting.gtf" ];
    then
        OUTPUT_FILE=$OUTPUT_DIR/featureCounts_Stats_RefSeq_Intronic_Only_GTF.txt
    fi


    rm -rf ${OUTPUT_FILE} && touch ${OUTPUT_FILE}


    #Print header to output file:
    #------------------------------------------
    if [ "${ANNOTATION_FILE}" == "Intron_Only_Regions.gtf" -o "${ANNOTATION_FILE}" == "intronic_only_gene_models_ncRNA_for_counting.gtf" ];
    then
        echo 'SAMPLE_ID' $'\t''DESCRIPTION' $'\t''Aligned_Pairs' $'\t''Read_Pair_IN_INTRONS' $'\t''Read_Pair_IN_INTRONS_RATIO' >> $OUTPUT_FILE
    fi
    #------------------------------------------
    if [ "${ANNOTATION_FILE}" == "genes.gtf" -o "${ANNOTATION_FILE}" == "RefSeq_GeneBody.gtf" -o "${ANNOTATION_FILE}" == "Exon_Only_Regions.gtf" -o "${ANNOTATION_FILE}" == "ncRNA_exon_for_counting.gtf" -o "${ANNOTATION_FILE}" == "ncRNA_genebodies_for_counting.gtf" -o  "${ANNOTATION_FILE}" == "exonic_only_gene_models_ncRNA_for_counting.gtf" ];
    then
        echo 'SAMPLE_ID' $'\t''DESCRIPTION' $'\t''Aligned_Pairs' $'\t''Read_Pair_IN_EXONS' $'\t''Read_Pair_IN_EXONS_RATIO' >> $OUTPUT_FILE
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

        cd ${DATASET_DIR}/$sample_id/fastq/tophat2
        #------------------------------------------
        #Get counts from *_align_summary.txt file
        # echo 'Getting the total number of mapped reads...'
        # #(statistics_for_all_accepted_reads.txt)
        # TOTAL_MAPPED_READS=$(grep 'total' statistics_for_all_accepted_reads.txt  | awk '{print $1}')
        #------------------------------------------
        #For paired-end data, we used the (-p) option:
        #-p        	If specified, fragments (or templates) will be counted instead
        #              	of reads. This option is only applicable for paired-end reads.
        #              	The two reads from the same fragment must be adjacent to each
        #              	other in the provided SAM/BAM file.
        
        echo 'Getting the Aligned_Pairs...'

        #Need to get the number of Aligned pairs:
        Aligned_Pairs=$(grep 'Aligned pairs:'  $sample_id'_align_summary.txt'  | awk '{print $3}')
        
        echo 'Getting READS_IN_EXONS...'
        #Need to cd to sample specific featureCounts folder
        #Need if statements to process the *_featureCounts.out file depending on the ANNOTATION_FILE used:
        featureCounts_DIR=${DATASET_DIR}/$sample_id/fastq/tophat2/featureCounts
        
        if [ "${ANNOTATION_FILE}" == "genes.gtf" ];
        then
            cd ${featureCounts_DIR}/Illumina_GTF
        fi
        
        if [ "${ANNOTATION_FILE}" == "RefSeq_GeneBody.gtf" ];
        then
            cd ${featureCounts_DIR}/RefSeq_Exon_GTF
        fi

        if [ "${ANNOTATION_FILE}" == "Intron_Only_Regions.gtf" ];
        then
            cd ${featureCounts_DIR}/RefSeq_Intron_GTF
        fi
        
        if [ "${ANNOTATION_FILE}" == "Exon_Only_Regions.gtf" ];
        then
            cd ${featureCounts_DIR}/RefSeq_Exon_Only_GTF
        fi

        #Start of lncRNA gtf
        if [ "${ANNOTATION_FILE}" == "ncRNA_exon_for_counting.gtf" ];
        then
            cd ${featureCounts_DIR}/LncRNA_Exon_Collapsed_GTF
        fi

        if [ "${ANNOTATION_FILE}" == "ncRNA_genebodies_for_counting.gtf" ];
        then
            cd ${featureCounts_DIR}/LncRNA_GeneBody_GTF
        fi

        if [ "${ANNOTATION_FILE}" == "intronic_only_gene_models_ncRNA_for_counting.gtf" ];
        then
            cd ${featureCounts_DIR}/LncRNA_Intronic_Only_GTF
        fi

        if [ "${ANNOTATION_FILE}" == "exonic_only_gene_models_ncRNA_for_counting.gtf" ];
        then
            cd ${featureCounts_DIR}/LncRNA_Exonic_Only_GTF
        fi

        #Get the sum of the 2nd coulumn
        READS_IN_EXONS=$(awk '{n+=$2;} ; END {print n;}' $sample_id'_featureCounts.out')
        echo 'Getting READS_IN_EXONS_RATIO...'

        #Calculate percentages
        READS_IN_EXONS_RATIO=$(echo "scale=4;$READS_IN_EXONS/$Aligned_Pairs" | bc)
        echo 'Printing to output file...'
        echo $sample_id $'\t'$description $'\t'$Aligned_Pairs $'\t'$READS_IN_EXONS $'\t'$READS_IN_EXONS_RATIO >> $OUTPUT_FILE
        #End loop over Sample_Labels:
        popd
    done 

    ##################################################################################
    #Need to summarize featureCounts special counters
    ##################################################################################

    ##################################################################################
    #Need if statements to determine OUTPUT_FILE_SPECIAL depending on the ANNOTATION_FILE used:

    if [ "${ANNOTATION_FILE}" == "genes.gtf" ];
    then
        OUTPUT_FILE_SPECIAL=$OUTPUT_DIR/featureCounts_summary_Illumina_GTF.txt
    fi

    if [ "${ANNOTATION_FILE}" == "RefSeq_GeneBody.gtf" ];
    then
        OUTPUT_FILE_SPECIAL=$OUTPUT_DIR/featureCounts_summary_RefSeq_Exon_GTF.txt
    fi

    if [ "${ANNOTATION_FILE}" == "Intron_Only_Regions.gtf" ];
    then
        OUTPUT_FILE_SPECIAL=$OUTPUT_DIR/featureCounts_summary_RefSeq_Intron_GTF.txt
    fi

    if [ "${ANNOTATION_FILE}" == "Exon_Only_Regions.gtf" ];
    then
        OUTPUT_FILE_SPECIAL=$OUTPUT_DIR/featureCounts_summary_RefSeq_Exon_Only_GTF.txt
    fi

    #Start of lncRNA gtf
    if [ "${ANNOTATION_FILE}" == "ncRNA_exon_for_counting.gtf" ];
    then
        OUTPUT_FILE_SPECIAL=$OUTPUT_DIR/featureCounts_summary_LncRNA_Exon_Collapsed_GTF.txt
    fi

    if [ "${ANNOTATION_FILE}" == "ncRNA_genebodies_for_counting.gtf" ];
    then
        OUTPUT_FILE_SPECIAL=$OUTPUT_DIR/featureCounts_summary_LncRNA_GeneBody_GTF.txt
    fi

    if [ "${ANNOTATION_FILE}" == "intronic_only_gene_models_ncRNA_for_counting.gtf" ];
    then
        OUTPUT_FILE_SPECIAL=$OUTPUT_DIR/featureCounts_summary_LncRNA_Intronic_Only_GTF.txt
    fi

    if [ "${ANNOTATION_FILE}" == "exonic_only_gene_models_ncRNA_for_counting.gtf" ];
    then
        OUTPUT_FILE_SPECIAL=$OUTPUT_DIR/featureCounts_summary_LncRNA_Exonic_Only_GTF.txt
    fi

    rm -rf ${OUTPUT_FILE} && touch ${OUTPUT_FILE}

    
    # samples contains array of (sample_dir, sample_id, description) for each sample
    samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

    # loop over all samples
    for ((i=0;i< ${#samples[@]} ;i+=3));
    do
        sample_dir=${samples[i]}
        sample_id=${samples[i+1]}
        description=${samples[i+2]}

        cd ${DATASET_DIR}/$sample_id/fastq/tophat2
        #------------------------------------------
        #Get counts from *_align_summary.txt file
        # echo 'Getting the total number of mapped reads...'
        # #(statistics_for_all_accepted_reads.txt)
        # TOTAL_MAPPED_READS=$(grep 'total' statistics_for_all_accepted_reads.txt  | awk '{print $1}')
        #------------------------------------------
        #For paired-end data, we used the (-p) option:
        #-p        	If specified, fragments (or templates) will be counted instead
        #              	of reads. This option is only applicable for paired-end reads.
        #              	The two reads from the same fragment must be adjacent to each
        #              	other in the provided SAM/BAM file.

        echo 'Getting the Aligned_Pairs...'
        #Need to get the number of Aligned pairs:
        Aligned_Pairs=$(grep 'Aligned pairs:'  $sample_id'_align_summary.txt'  | awk '{print $3}')

        #Need to cd to sample specific featureCounts folder
        #Need if statements to process the *_featureCounts.out file depending on the ANNOTATION_FILE used:
        featureCounts_DIR=${DATASET_DIR}/$sample_id/fastq/tophat2/featureCounts

        if [ "${ANNOTATION_FILE}" == "genes.gtf" ];
        then
            cd ${featureCounts_DIR}/Illumina_GTF
        fi

        if [ "${ANNOTATION_FILE}" == "RefSeq_GeneBody.gtf" ];
        then
            cd ${featureCounts_DIR}/RefSeq_Exon_GTF
        fi

        if [ "${ANNOTATION_FILE}" == "Intron_Only_Regions.gtf" ];
        then
            cd ${featureCounts_DIR}/RefSeq_Intron_GTF
        fi

        if [ "${ANNOTATION_FILE}" == "Exon_Only_Regions.gtf" ];
        then
            cd ${featureCounts_DIR}/RefSeq_Exon_Only_GTF
        fi

        #Start of lncRNA gtf
        if [ "${ANNOTATION_FILE}" == "ncRNA_exon_for_counting.gtf" ];
        then
            annotation_dir=${featureCounts_DIR}/LncRNA_Exon_Collapsed_GTF
        fi

        if [ "${ANNOTATION_FILE}" == "ncRNA_genebodies_for_counting.gtf" ];
        then
            annotation_dir=${featureCounts_DIR}/LncRNA_GeneBody_GTF
        fi

        if [ "${ANNOTATION_FILE}" == "intronic_only_gene_models_ncRNA_for_counting.gtf" ];
        then
            annotation_dir=${featureCounts_DIR}/LncRNA_Intronic_Only_GTF
        fi

        if [ "${ANNOTATION_FILE}" == "exonic_only_gene_models_ncRNA_for_counting.gtf" ];
        then
            annotation_dir=${featureCounts_DIR}/LncRNA_Exonic_Only_GTF
        fi

        pushd "${annotation_dir}"
        
        #------------------------------------------
        #Omit header line:
        #------------------------------------------
        if [ "${ANNOTATION_FILE}" == "RefSeq_GeneBody.gtf" ];
        then
            cp $sample_id'_featureCounts.out.summary' temp1.txt
        else
            awk 'NR>1' $sample_id'_featureCounts.out.summary' > temp1.txt
        fi
        #------------------------------------------
        #Divide 2nd column by total number of reads:
        awk -v Aligned_Pairs=${Aligned_Pairs} '{ Ratio=($2)/(Aligned_Pairs) ; print $1"\t"$2"\t"Aligned_Pairs"\t"Ratio }' temp1.txt > temp2.txt
        #Need a header:
        echo 'Status' $'\t'$sample_id'_Read_Count' $'\t''Aligned_Pairs' $'\t''Ratio' >> Header.txt
        #Combine files:
        cat Header.txt temp2.txt > temp3.txt
        #------------------------------------------
        if [ "${ANNOTATION_FILE}" == "RefSeq_GeneBody.gtf" ];
        then
            sed 's/.out.*/:/' temp3.txt > temp4.txt
            mv temp4.txt temp3.txt
        fi
        #------------------------------------------
        echo 'Concatenating summary file...'
        cat temp3.txt >> $OUTPUT_FILE_SPECIAL
        #Remove temp files:
        rm temp*.txt
        rm Header.txt
        ##------------------------------------------
        #End loop over Sample_Labels:
        popd
    done 

    echo 'Check out '$OUTPUT_FILE
    echo 'Check out '$OUTPUT_FILE_SPECIAL
    #End loop over GTF files:
done
echo 'Check out: '${OUTPUT_DIR}
echo '#--------------------------------------------------------------------------'
##################################################################################
