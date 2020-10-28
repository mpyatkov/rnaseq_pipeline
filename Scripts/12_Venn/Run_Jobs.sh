#!/bin/bash -l

# creating 3 types of Venn diagrams

# venn_for_one_gtf_set: 
# pairwise comparison of DE results (ex. 09a_1 vs 09a_2)

# individual_comparison:
# create 2 venn diagrams for each individual folder (ex.09a_1)
# ex. venn diagram will be created for each feature (ExonCollapsed/IntronOnly...) set
# 1 type of venn - deseq/edger for each feature separately
# 2 type of venn - all feature together

set -o errexit
set -o pipefail
set -o nounset

## create venn diagrams for comparisons

module load R/3.6.0

#Remove *.o files from previous jobs
# rm -rf ./logs && mkdir -p ./logs

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

SCRIPT_DIR="$(pwd)"
OUTDIR=Job_Summary
rm -rf ${OUTDIR} && mkdir -p ${OUTDIR}

# 09a, 09b, 09c, 09d just different sets of GTF files and counters
# 09a - 3 refseq gtf + htseq counter
# 09b - 3 refseq gtf + featureCount counter
# 09c - 4 lnc15k gtf + featureCount counter
# 09d - 2 refseq_lnc76k gtf + featureCount counter

# get number of venn comparisons
venn_number=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --venn_number))

# for pairwise comparison (ex. 09a_1 vs 09a_2,...)
venn_for_one_gtf_set() {
    local DE_INDEX=$1
    local OUTDIR_PREFIX=$2
    local DETITLE=$3
    # not individual diagrams
    venn_individ=0
    output_dir="${OUTDIR}/DE_Comparisons/${OUTDIR_PREFIX}"

    # loop through the all venn comparisons
    mkdir -p ${output_dir}

    pushd ${output_dir}
    for((i=0;i<venn_number;i+=1));
    do
        # array of comparisons which involved in current venn
        comparisons=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --venn_comparisons_by_ix ${i}))

        # outname for pdf file
        outname=${OUTDIR_PREFIX}_$(IFS=_ ; echo "${comparisons[*]}")
        
        # copy ExonCollapsed file from each comparison directory 
        for comparison_id in "${comparisons[@]}"
        do
	    # echo "${comparison_id}"
	    search_body="${SCRIPT_DIR}/../${DE_INDEX}_DE_${comparison_id}_* -name 'DiffEx*ExonCollapsed*txt'"
	    defiles=$(eval find ${search_body})
	    for fname in $defiles
	    do
	        # make a link instead of copy
	        # ln -s $fname "./${comparison_id}_$(basename $fname)"
	        cp $fname "./${comparison_id}_$(basename $fname)"
	    done
        done

        # ${DATASET_LABEL} - exported from config file
        cp ${SCRIPT_DIR}/Venn.R ./
        Rscript Venn.R "${DATASET_LABEL}_${DETITLE}" "${outname}" ${venn_individ}
        rm *_DiffExp_v2* Venn.R
    done
    popd
    # mv ${OUTDIR_PREFIX} ./Job_Summary
}

# find ../ -name "09*" -maxdepth 1 -type d | xargs -n1 basename
individual_comparison() {
    local current_dir=$1
    outname=${current_dir}
    venn_individ=1

    mkdir -p ${current_dir}
    pushd ${current_dir}

    # COPY ALL DiffEx_v2 files to the current_dir
    # extract comparison number (09d_1_... --> 1)
    prefix=$(echo ${current_dir} | grep -Po "_\K([0-9][0-9]?)(?=_)" | head -n1 |tr -d '\n')
    echo $current_dir
    # search body for extracting all files started with DiffEx_v2... (ExonCollapsed, IntronOnly,...)
    search_body="${SCRIPT_DIR}/../${current_dir} -name 'DiffExp_v2*txt'"

    # DE files only from ${current_dir}
    defiles=$(eval find ${search_body})
    
    # for each fname make 2 plots
    for fname in $defiles
    do
	cp ${SCRIPT_DIR}/Venn.R ./
	cp $fname "./${prefix}_$(basename $fname)"
    done
    
    # make individual plots
    Rscript Venn.R "${DATASET_LABEL}_${current_dir}" "${outname}" ${venn_individ}
    rm *_DiffExp_v2* Venn.R

    popd
}

# wrapper for individual_comparison function
individual_comparisons_wrapper() {
    individual_dirs=$OUTDIR/Individual

    mkdir -p ${individual_dirs}

    pushd ${individual_dirs}

    search_body="${SCRIPT_DIR}/../ -maxdepth 1 -name '09*' -type d"
    all_dirs=$(eval find ${search_body} | xargs -n1 basename)

    for dir in ${all_dirs}
    do
	individual_comparison $dir
    done
    popd
}

# create individual venn diagrams for each 09 directory
individual_comparisons_wrapper

# create venn diagrams for each set of gtf+counter (09a,b,c,d)
venn_for_one_gtf_set 09a 12a_Venn_RefSeq24_HTSeq_ExonCollapsed RefSeq24k
venn_for_one_gtf_set 09b 12b_Venn_RefSeq24_featureCounts_ExonCollapsed RefSeq24k
venn_for_one_gtf_set 09c 12c_Venn_LncRNA15k_featureCounts_ExonCollapsed LncRNA15k
venn_for_one_gtf_set 09d 12d_Venn_RefSeqLncRNA76k_featureCounts_ExonCollapsed RefSeqLncRNA76k
