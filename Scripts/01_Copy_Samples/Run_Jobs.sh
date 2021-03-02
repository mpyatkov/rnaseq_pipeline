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
for ((i=0;i< ${#samples[@]} ;i+=3));
do
    
    SAMPLE_ID=${samples[i+1]}
    SAMPLE_DIR="${DATASET_DIR}/${SAMPLE_ID}"

    # to prevent program interruption when extracting info about sample
    # we will use (set +eu) option
    set +eu
    # TODO ERR: catch the errors and move them to the log file
    # sample_info=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --get_sample_info ${SAMPLE_ID} 2>./logs/python_log.err))
    sample_info=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --get_sample_info ${SAMPLE_ID}))

    # if exit status is 1 this means that script raise exception and
    # SAMPLE_ID was not found in the local or global indexes
    # last exit code --> '$?', 0 -- sample is in the index, 1 -- sample is not in the index
    exit_code=$?
    if [ ${exit_code} -ne 0 ]; then
	sample_not_in_db=true
	# TODO ERR: to show error message
	# echo ${sample_info[@]}
    else
	sample_not_in_db=false
    fi
    set -eu

    ## 'create' sample dir, if already exist do nothing
    mkdir -p ${SAMPLE_DIR}

    pushd ${SAMPLE_DIR}

    ### CHECK IF FASTQ FILES IS CREATED
    
    ## if number of fastq files not equal 2, we will try to recreate them
    number_of_fastq=$(find . -maxdepth 1 -name "*.f*q.gz" | wc -l)

    ## if number of fastq files are not equal 2 and sample_not_in_db than exit
    if [ ${number_of_fastq} -ne 2 -a ${sample_not_in_db} = true ];then
	echo "ERROR: Cannot find fastq files inside the ${SAMPLE_DIR}"
	exit 1
    fi

    ## if sample in db and it is reqired to refresh fastq files
    if [ ${number_of_fastq} -ne 2 ]; then
	echo "Recreating FASTQ links"
	## remove previous FASTQ files
	rm -rf ${SAMPLE_DIR}/*.f*q.gz

	## get reads from index and make links to fastq files
	read1=${sample_info[1]}
	read2=${sample_info[2]}
	ln -s $read1 ./ && ln -s $read2 ./
    fi

    ### CHECK IF WE CAN COPY PRECALCULATED BAM FILE FROM THE SERVER

    ## if bam is not presented try to get it from waxmalabvm
    if [ ! -f "${SAMPLE_DIR}/aligner/${SAMPLE_ID}_primary_unique.bam" ]; then
	## get the potential CRAM_PATH
	if [ ${sample_not_in_db} = true ]; then
	    CRAM_PATH="${VM_DIR_UCSC}/NON-INDEXED/${SAMPLE_ID}_primary_unique.cram"
	    CRAI_PATH="${VM_DIR_UCSC}/NON-INDEXED/${SAMPLE_ID}_primary_unique.crai"
	else
	    project_name=${sample_info[0]}
	    CRAM_PATH="${VM_DIR_UCSC}/INDEXED_PROJECTS/${project_name}/${SAMPLE_ID}_primary_unique.cram"
	    CRAI_PATH="${VM_DIR_UCSC}/INDEXED_PROJECTS/${project_name}/${SAMPLE_ID}_primary_unique.crai"
	fi

	## check if CRAM exists
	if [ -f ${CRAM_PATH} ]; then
	    echo "Creating link to CRAM file"
	    mkdir -p aligner
	    pushd aligner
	    rm -rf ./*.cra*
	    ln -s ${CRAM_PATH} ./
	    ln -s ${CRAI_PATH} ./
	    
	    (set -x; qsub -N "${job_name}_${SAMPLE_ID}" -P "${PROJECT}" -l h_rt="00:30:00" -wd ${SCRIPT_DIR} ${SCRIPT_DIR}/unpack_cram.qsub ${SAMPLE_ID})
	    popd
	fi
    fi
        
    popd

done

echo "End of 01_Copy_Samples commands"
echo "-----------------------"
