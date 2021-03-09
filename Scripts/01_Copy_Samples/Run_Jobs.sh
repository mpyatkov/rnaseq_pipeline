#!/bin/bash

# create directories for samples using the FASTQ index file
# if sample dir exists check the availability of FASTQ and BAM files

set -o errexit
set -o pipefail
set -o nounset

#Remove *.o files from previous jobs
rm -rf ./logs && mkdir -p ./logs

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

## array of all samples; (sample_dir, sample_id; condname) for each sample
samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))

# dir_name and job_name are required in the next steps
dir_name=$(basename $(pwd))
step_num=$(echo ${dir_name} | cut -d'_' -f 1)
job_name="Step_${step_num}"

SCRIPT_DIR=$(pwd)

## loop over the all samples
for ((i=0;i< ${#samples[@]} ;i+=2));
do
    
    SAMPLE_ID=${samples[i]}
    SAMPLE_DIR="${DATASET_DIR}/${SAMPLE_ID}"

    # get_sample_info return (PRJ_NAME, READ1, READ2) or EXCEPTION
    sample_info=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --get_sample_info ${SAMPLE_ID}))

    ## 'create' sample dir, if already exist do nothing
    mkdir -p ${SAMPLE_DIR}

    pushd ${SAMPLE_DIR}

    ### CHECK IF FASTQ FILES IS CREATED
    
    ## if number of fastq files not equal 2, we will try to recreate them
    number_of_fastq=$(find . -maxdepth 1 -name "*.f*q.gz" | wc -l)

    ## if sample in db and it is required to refresh fastq files
    if [ ${number_of_fastq} -ne 2 ]; then
	echo "Recreating FASTQ links"
	## remove previous FASTQ files
	rm -rf ${SAMPLE_DIR}/*.f*q.gz

	## get reads from index and make links to fastq files
	read1=${sample_info[1]}
	read2=${sample_info[2]}
	
	if [ ! -f ${read1} -o ! -f ${read2} ]; then
	    echo "ERROR: cannot find one of the files with reads:\n${read1}\n${read2}\n"
	    exit 1
	fi
	ln -s $read1 ./ && ln -s $read2 ./
    fi

    ### CHECK IF WE CAN COPY PRECALCULATED BAM FILE FROM THE SERVER

    ## if bam is not presented try to get it from waxmalabvm
    if [ ! -f "${SAMPLE_DIR}/aligner/${SAMPLE_ID}_primary_unique.bam" ]; then

	## get the potential CRAM_PATH
	project_name=${sample_info[0]}
	CRAM_PATH="${VM_DIR_UCSC}/INDEXED_PROJECTS/${project_name}/${SAMPLE_ID}_primary_unique.cram"
	CRAI_PATH="${VM_DIR_UCSC}/INDEXED_PROJECTS/${project_name}/${SAMPLE_ID}_primary_unique.crai"

	## check if CRAM exists
	if [ -f ${CRAM_PATH} ]; then
	    echo "Creating link to CRAM file"
	    mkdir -p aligner
	    pushd aligner
	    rm -rf ./*.cra*
	    ln -s ${CRAM_PATH} ./
	    ln -s ${CRAI_PATH} ./
	    
	    (set -x; qsub -N "${job_name}_${SAMPLE_ID}" -P "${PROJECT}" -l h_rt="00:45:00" -wd ${SCRIPT_DIR} ${SCRIPT_DIR}/unpack_cram.qsub ${SAMPLE_ID})
	    popd
	fi
    fi
        
    popd

done

echo "End of 01_Copy_Samples commands"
echo "-----------------------"
