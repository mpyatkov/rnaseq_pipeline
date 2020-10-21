#!/bin/bash -l

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

venn_for_one_gtf_set() {
    local DE_INDEX=$1
    local OUTDIR_PREFIX=$2
    local DETITLE=$3
    
    # loop through the all venn comparisons
    mkdir ${OUTDIR_PREFIX}

    pushd ${OUTDIR_PREFIX}
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
        cp ../Venn.R ./
        Rscript Venn.R "${DATASET_LABEL}_${DETITLE}" "${outname}.pdf"
        rm *_DiffExp_v2* Venn.R
    done
    popd
    mv ${OUTDIR_PREFIX} ./Job_Summary
}

# create venn diagrams for each set of gtf+counter (09a,b,c,d)
venn_for_one_gtf_set 09a 12a_Venn_RefSeq24_HTSeq_ExonCollapsed RefSeq24k
venn_for_one_gtf_set 09b 12b_Venn_RefSeq24_featureCounts_ExonCollapsed RefSeq24k
venn_for_one_gtf_set 09c 12c_Venn_LncRNA15k_featureCounts_ExonCollapsed LncRNA15k
venn_for_one_gtf_set 09d 12d_Venn_RefSeqLncRNA76k_featureCounts_ExonCollapsed RefSeqLncRNA76k


