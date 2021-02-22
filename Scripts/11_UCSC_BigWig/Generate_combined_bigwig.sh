#!/bin/bash -l

# Specify which shell to use
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -e ./logs/
#$ -o ./logs/
#$ -pe omp 2

# This script should allow you to generate combined bigwig files for
# groups of previously obtained Forward and Reverse bigwig files

# At this step that all required Forward and Reverse bigwig files are located
# in the /net/waxman-server/mnt/data/waxmanlabvm_home/TRACKS/COMMON
# directory.

# set strict options
set -o errexit
set -o pipefail
set -o nounset

# Extract samples names and combine them to appropriate lines, ex: G186_M1M2M3
# This procedure works only if your sample folders called as "PROJECTNAME_SAMPLEID",
# because here we extract all SAMPLEIDs and combine them as one line
samples_line(){
    local fname=$1
    a=$(tail -n+2 $fname | cut -f2 | cut -d "_" -f 1 | sort | uniq)

    res=""
    for i in $a
    do
	group="${i}_$(grep ${i} $fname | cut -f 2 | cut -d "_" -f 2 | paste -s -d '')_"
	res+=$group
	# res+="_"
    done

    echo ${res%_}
}

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

# activate conda environment
set +eu
module load miniconda
conda activate --stack ${CONDA_DIR}/wig
set -eu

# dir name for storing bigwig files 
# UCSC_SAMPLE=UCSC_BigWig

# remove file from previous run which will contains group names
rm -rf COMBINED_PAIRS.txt

groups=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --groups))

# for each group
for group in "${groups[@]}"; do

    # get samples array (triples of sample_dir sample_id description for each sample)
    samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples_by_group ${group}))

    # create new array contains path to Forward and Reverse bigwig files
    forward=()
    reverse=()
    # need to extract sample names and combine them together (ex: M1M2M3)
    suffix=()

    group_name=""
    for ((i=0;i< ${#samples[@]} ;i+=3));
    do
	SAMPLE_ID=${samples[i+1]}

	# because we use --samples_by_group with a group name
	# 01_Pipeline_Setup.py give us condition_name column instead
	# of description column. Condition_name should be the same for all
	# members of one group
	CONDNAME=${samples[i+2]}
	group_name=$CONDNAME
	
	# TODO: at the moment this procedure supposes that we make
	# calculations for all groups, thus for all samples should be
	# created bigwig files. It is required to check if these files
	# are exist (forward and reverse). 
	
	# forward+=("${DATASET_DIR}/${SAMPLE_ID}/${UCSC_SAMPLE}/${SAMPLE_ID}.Forward.bw")
	# reverse+=("${DATASET_DIR}/${SAMPLE_ID}/${UCSC_SAMPLE}/${SAMPLE_ID}.Reverse.bw")

	forward+=("${VM_DIR_UCSC}/COMMON/${SAMPLE_ID}.Forward.bw")
	reverse+=("${VM_DIR_UCSC}/COMMON/${SAMPLE_ID}.Reverse.bw")
	
	suffix+=("${SAMPLE_ID}")
    done
    
    # make tmp file for names (it is easy to work with file than array)
    printf "%s\n" "${suffix[@]}" > tmpfile
    suffix_combined=$(samples_line tmpfile)
    rm tmpfile

    # skip the rest if combined file already exist in ${VM_DIR_UCSC}/COMMON/
    FORWARD_BW="${suffix_combined}.Forward.combined.bw"
    REVERSE_BW="${suffix_combined}.Reverse.combined.bw"

    # if both files created skip calculation
    if [ ! -f "${VM_DIR_UCSC}/COMMON/${FORWARD_BW}" -o ! -f "${VM_DIR_UCSC}/COMMON/${REVERSE_BW}" ]; then
	# run wiggletools to calculate average for bigwig files and
	# convert them back to bigwig again
	echo "Forward reads ${suffix_combined}"
	# Forward reads
	wiggletools mean ${forward[@]} > tmp_forward.wig
	wigToBigWig tmp_forward.wig "./Chrom_Sizes/mm9.chrom.sizes" ${FORWARD_BW}
	rm tmp_forward.wig

	echo "Reverse reads ${suffix_combined}"
	wiggletools mean ${reverse[@]} > tmp_reverse.wig
	wigToBigWig tmp_reverse.wig "./Chrom_Sizes/mm9.chrom.sizes" ${REVERSE_BW}
	rm tmp_reverse.wig
    fi

    # add filename to COMBINED_PAIRS.txt (this file is required for track lines(next step))
    
    # print "filename \t group_name for forward reads"
    printf "%s\t%s\n" "${FORWARD_BW}" "${group_name}" >> COMBINED_PAIRS.txt

    # print "filename \t group_name for reverse reads"
    printf "%s\t%s\n" "${REVERSE_BW}" "${group_name}" >> COMBINED_PAIRS.txt

done

# move all combined files to the server if it is required
echo "move files to the server"
set +eu
mv *.bw ${VM_DIR_UCSC}/COMMON/
set -eu

# generate track lines
./Generate_Tracks.sh

echo "-----"
echo "IAMOK"
