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

rm -rf *.o* *.e*

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

# load anaconda in case if run this script independently
(
    set +eu
    module load anaconda2
    source activate RNAseq
)

# calculate strandedness_featurecount
if [ ${STRANDEDNESS} -eq 0 ]
then
    strandedness_featurecount="0"
elif [ ${STRANDEDNESS} -eq 2 ]
then
    strandedness_featurecount="1"
elif [ ${STRANDEDNESS} -eq 1 ]
then
    strandedness_featurecount="2"
fi

feature_id=gene_id

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

    #---------------------------
    #Gene Annotation file to use
    #---------------------------------------------------------------------------------
    #For read mapping purposes we want the full gene body information (use the
    #RefSeq_GeneBody.gtf)
    #Presence of splice junctions should be based on the full gene structure for
    #any isoform of a gene symbol
    #For read counting purposes we can use either 
    #1. RefSeq_GeneBody.gtf
    #2. Intron_Only_Regions.gtf
    #3. Exon_Only_Regions.gtf
    #---------------------------------------------------------------------------------
    #Start loop over GTF files:
    gtf_list=$GTF_FILES_DIR/*.gtf
    for gtf_file in ${gtf_list}
    do
        ANNOTATION_FILE=$(basename ${gtf_file})

        #Additional variables will be populated from if statements
    
        #How reads are counted
        #--------------------------
        # if EXONS_ONLY is set to be true, then, only reads that overlap exons
        # will be counted (recommended)
        # if EXONS_ONLY is set to be false,then, all reads that overlap a gene
        # will be counted
        # Need if statements to determine FEATURE_TYPE depending on the
        # ANNOTATION_FILE used:

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
        
        if [ "${ANNOTATION_FILE}" == "Intron_Only_Regions.gtf" ];
        then
            #The 3rd column in the Intron_Only_Regions.gtf is "Intronic_Only"
            FEATURE_TYPE="Intronic_Only"
        fi

        if [ "${ANNOTATION_FILE}" == "Exon_Only_Regions.gtf" ];
        then
            #The 3rd column in the Exon_Only_Regions.gtf is "Exonic_Only"
            FEATURE_TYPE="Exonic_Only"
        fi

        #Start of lncRNA gtf
        if [ "${ANNOTATION_FILE}" == "ncRNA_exon_for_counting.gtf" ];
        then
            #Counting lncRNA by their exon (exon_collapsed)
            FEATURE_TYPE="exon"
        fi

        if [ "${ANNOTATION_FILE}" == "ncRNA_genebodies_for_counting.gtf" ];
        then
            #Counting lncRNA by their gene body
            FEATURE_TYPE="gene_body"
        fi
        
        if [ "${ANNOTATION_FILE}" == "intronic_only_gene_models_ncRNA_for_counting.gtf" ];
        then
            #Counting lncRNA by their intronic only regions
            FEATURE_TYPE="Intronic_Only"
        fi

        if [ "${ANNOTATION_FILE}" == "exonic_only_gene_models_ncRNA_for_counting.gtf" ];
        then
            #Counting lncRNA by their exonic only regions
            FEATURE_TYPE="Exonic_Only"
        fi
        
        (set -x; qsub -N "${job_name}_${sample_id}" -P "${PROJECT}" -l h_rt="${TIME_LIMIT}" Extract_Counts.qsub ${sample_id} ${strandedness_featurecount} ${feature_id} ${ANNOTATION_FILE} ${FEATURE_TYPE})
        
        #End loop over GTF files:
    done
    #End loop over Sample_Labels:
done 
echo "End of 08b qsub commands"
echo "-----------------------"
##################################################################################
