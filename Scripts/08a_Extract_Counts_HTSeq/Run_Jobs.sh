#!/bin/bash -l

set -o errexit
set -o pipefail
set -o nounset

##################################################################################
#Andy Rampersaud, 02.22.16
#This script would be used to run Extract_Counts.pbs in parallel for different FASTQ files
#Way to run script:
#Usage: ./Run_Jobs.sh
#Example: 
#./Run_Jobs.sh 
##################################################################################

rm -rf ./logs && mkdir -p ./logs

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

# load anaconda module for case when we need independent run
set +eu
module load anaconda2
source activate RNAseq
set -eu

#We always count by gene_id when using HTSeq
feature_id=gene_id

# calculate MODE
if [[ ${MODE} -eq 0 ]]
then
    mode="union"
elif [[ ${MODE} -eq 1 ]]
then
    mode="intersection-strict"
elif [[ ${MODE} -eq 2 ]]
then
    mode="intersection-nonempty"
fi

# calculate strandedness_htseq from STRANDEDNESS
if [ ${STRANDEDNESS} -eq 0 ]
then
    strandedness_htseq="no"
elif [ ${STRANDEDNESS} -eq 2 ]
then
    strandedness_htseq="yes"
elif [ ${STRANDEDNESS} -eq 1 ]
then
    strandedness_htseq="reverse"
fi


# dir_name and job_name are required in the next steps
dir_name=$(basename $(pwd))
step_num=$(echo ${dir_name} | cut -d'_' -f 1)
job_name="Step_${step_num}"


# samples contains array of (sample_dir, sample_id, description) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

# loop over all samples
for ((i=0;i< ${#samples[@]} ;i+=3));
do
    sample_dir=${samples[i]}
    sample_id=${samples[i+1]}
    description=${samples[i+2]}

    #Start loop over GTF files:
    gtf_list=${GTF_FILES_DIR}/*.gtf

    for gtf_file in $gtf_list
    do
        ANNOTATION_FILE=$(basename "${gtf_file}")

        # Please: DO NOT EDIT CODE BELOW
        # Additional variables will be populated from if statements

        # Extracting counts:
        # How reads are counted
        #--------------------------
        # if EXONS_ONLY is set to be true, then, only reads that overlap exons
        # will be counted (recommended)
        # if EXONS_ONLY is set to be false,then, all reads that overlap a gene
        # will be counted
        # Need if statements to determine FEATURE_TYPE depending on the
        #------------------------------------------

        #Use the -o (logical "or")
        if [ "${ANNOTATION_FILE}" == "genes.gtf" -o "${ANNOTATION_FILE}" == "RefSeq_GeneBody.gtf" ];
        then
            EXONS_ONLY="TRUE"
            #If statement to determine FEATURE_TYPE:
            if [ $EXONS_ONLY == "TRUE" ] 
            then
	        FEATURE_TYPE="exon"
            else
	        FEATURE_TYPE="CDS"
            fi
        fi

        #------------------------------------------
        if [ "${ANNOTATION_FILE}" == "Intron_Only_Regions.gtf" ];
        then
            #The 3rd column in the Intron_Only_Regions.gtf is "Intronic_Only"
            FEATURE_TYPE="Intronic_Only"
        fi

        #------------------------------------------
        if [ "${ANNOTATION_FILE}" == "Exon_Only_Regions.gtf" ];
        then
            #The 3rd column in the Exon_Only_Regions.gtf is "Exonic_Only"
            FEATURE_TYPE="Exonic_Only"
        fi

        #Start of lncRNA gtf, ignore all of these gtf
        if [ "${ANNOTATION_FILE}" == "ncRNA_exon_for_counting.gtf" -o "${ANNOTATION_FILE}" == "ncRNA_genebodies_for_counting.gtf" -o "${ANNOTATION_FILE}" == "intronic_only_gene_models_ncRNA_for_counting.gtf" -o "${ANNOTATION_FILE}" == "exonic_only_gene_models_ncRNA_for_counting.gtf" ];
        then
            echo "lncRNA GTF file is ignored by HTSeq, see FeatureCount result output"
        else

            #Now use arguments in the PBS script call:
            (set -x; qsub -N "${job_name}_${sample_id}" -P "${PROJECT}" -l h_rt="${TIME_LIMIT}" Extract_Counts.qsub ${sample_id} ${strandedness_htseq} ${feature_id} ${ANNOTATION_FILE} ${FEATURE_TYPE} ${mode})
        fi #End of if statement to ignore lncRNA gtfs
        #End loop over GTF files:
    done
    
done

echo "End of 8a step commands"
echo "-----------------------"

##################################################################################
