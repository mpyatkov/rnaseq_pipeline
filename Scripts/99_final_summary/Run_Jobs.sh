#!/usr/bin/bash -l

# This script summarizes data from the directories located above:
# - Copy "combined" pdf files from step 13
# - Copy all Segex files 
# - Summarize info about up/down genes by different features

set -o errexit
set -o pipefail
set -o nounset

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

set +eu
module load miniconda
conda activate --stack ${CONDA_DIR}/rlang361
conda activate --stack ${CONDA_DIR}/multiqc
set -eu

rm -rf *.txt
rm -rf *.pdf
rm -rf *.csv
rm -rf DiffExp_*
rm -rf output 

SCRIPT_DIR="$(pwd)"
cd ..

Level_UP=$(pwd)

mkdir ${Level_UP}/99_final_summary/output

cd ${SCRIPT_DIR}
echo "Script_Directory"
echo ${SCRIPT_DIR}
echo "Level_UP:"
echo ${Level_UP}

## remove bam/cram/bai/crai if they are not required
if [ ${KEEP_BAM_AFTER_RUN} -eq 0 ]; then
    extensions=( "*bam" "*bai" "*cram" "*crai" )
    for ext in ${extensions[@]}; do
	echo "removing '${ext}' files"
	set +eu
	find ${DATASET_DIR} -name ${ext} -delete
	set -eu
    done
fi

set +eu
# cp -rf  ${Level_UP}/04_TopHat_Paired_End/Job_Summary/TopHat2_Stats_BestMapped.txt  ./output/
cp -rf  ${Level_UP}/02_Read_Strandness/Job_Summary/Read_Strandness_Stats.txt ./output/
cp -rf  ${Level_UP}/03_FASTQC/Job_Summary/number_of_unique_reads.csv ./output/
cp -rf  ${Level_UP}/06_CollectMetrics/Job_Summary/CollectRnaSeqMetrics_Stats.txt  ./output/
cp -rf  ${Level_UP}/06_CollectMetrics/Job_Summary/CollectInsertSizeMetrics_Stats.txt  ./output/
cp -rf  ${Level_UP}/06_CollectMetrics/Job_Summary/CollectInsertSizeMetrics_Plots.pdf ./output/
cp -rf  ${Level_UP}/12_Venn/Job_Summary/*.pdf ./output/
cp -rf ${Level_UP}/00_Setup_Pipeline/*.xlsx ./output/

mkdir -p ./output/EnrichmentAnalysis
cp -rf ${Level_UP}/14_Enrichment/Job_Summary/* ./output/EnrichmentAnalysis/

mkdir -p ./output/Heatmaps
cp -rf ${Level_UP}/15_Heatmaps/Job_Summary/* ./output/Heatmaps/

# cp -rf  ${Level_UP}/08_Extract_Counts/Job_Summary/featureCounts_summary_LncRNA15k_ExonCollapsed_GTF.txt  ./output/
set -eu

mkdir -p output/Combined_PCA_and_Correlation_plots
find ${Level_UP}/13_Correlation/Job_Summary/ -name "*.pdf" | grep -E "Combined|SignificantOnly" | xargs -n1 -I{} cp {} ./output/Combined_PCA_and_Correlation_plots/

## SPECIFY FUNCTIONS

copy_feature(){
    local de_index=$1
    local feature=$2

    mkdir -p output/Segex_${de_index}/Segex${de_index}_${feature}

    # find all segex files associated with fpkm, edger, not located in
    # output, contained ${feature} in the name and copy them to the
    # required folder
    find ../${de_index}_DE_* -name "*SEGEX*" | grep -iv output | grep -i ${feature} | grep -i tpm | grep -i edger | xargs -n1 -I{} cp {} ./output/Segex_${de_index}/Segex${de_index}_${feature}/
}


# copy segex files from DE directories to step 99 output
copy_all_segex_files() {
    local de_index=$1

    ## copy exoncollapsed files to up level of each feature directory (deprecated, 2021-04-15)
    # mkdir -p output/Segex_${de_index}
    # find ../${de_index}_DE_* -name "*SEGEX*" | grep -iv output | grep -i fpkm | grep -i exoncollapsed | grep -i edger | xargs -n1 -I{} cp {} ./output/Segex_${de_index}/
    
    # copy features (ex. for index 9a -> ExonOnly, IntronOnly, ExonCollapsed, ...)
    features=$(find ../${de_index}* -iname "output*" | grep -Po "\K([a-zA-Z]*)(?=$)" | sort | uniq)
    
    # loop over each feature associated with de_index
    for feature in $features
    do
	echo "copy_feature: ${de_index} ${feature}"
	copy_feature ${de_index} ${feature}
    done
}


# make links for segex files by normalization method (fpkm/tpm)
link_all_segex_by_norm() {
    local de_index=$1
    local norm_method=$2 # fpkm/tpm
    local feature=$3

    find ../../${de_index}_DE_* -name "*SEGEX*" | grep -i ${feature} | grep -i edger |grep -i ${norm_method} | xargs -n1 -I{} ln -s {} ./
    
    # if [[ ${norm_method} == "TPM" ]]; then
    # 	find ../../${de_index}_DE_* -name "*SEGEX*" | grep -i ${feature} | grep -i edger |grep -i ${norm_method} | xargs -n1 -I{} ln -s {} ./	
    # else
    # 	find ../../${de_index}_DE_* -name "*SEGEX*" | grep -i ${feature} | grep -i edger |grep -iv "TPM" | xargs -n1 -I{} ln -s {} ./	
    # fi
    
}


# make links for segex files by the DE index
link_all_segex_by_ix() {
    local de_index=$1
    find ../../${de_index}_DE_* -name "*SEGEX*" | grep -i "TPM" | xargs -n1 -I{} ln -s {} ./
}


# summarize up and down genes
function updown_genes() {
    local de_index=$1
    
    # INPUT_UPDOWN_FILE="./gene_counts.csv"
    SAMPLES_FILE="../../00_Setup_Pipeline/Sample_Labels.txt"
    CMP_FILE="../../00_Setup_Pipeline/Comparisons.txt"
    OUTPUT_FILE="$(basename $(pwd)).txt"
    
    cp ../Scripts/updown_genes.R ./

    # two last parameters for updown_genes.R should be paired arrays
    # "c(2,1,1,4)" -- array with FC
    # "c(0.05,0.05,0.001,0.05)" -- array with FDR
    Rscript updown_genes.R ${CMP_FILE} ${SAMPLES_FILE} ${OUTPUT_FILE} ${DATASET_LABEL} "c(2,1,1,4)" "c(0.05,0.05,0.001,0.0001)"
    echo ${OUTPUT_FILE}
}

# make combined files for segex output
function segex_combined_files() {
    local DIR=$1

    output_header=""
    afnames=()
    
    for f in $(find . -name "*.txt" | sort -n); do
	fname=$(basename $f)
	# subdirname=$(basename $(dirname $DIR))
	subdirname=$(basename $DIR)
	
	# add file name to the array
	afnames+=($fname)
	
	# add file to the
	output_header+="$fname\t\t\t\t\t\t\t"
    done

    # combine all files
    output_fname="Combined_${subdirname}.txt"
    tmp_fname="tmp_fname"
    paste ${afnames[@]} > ${tmp_fname}

    # add header (only in the first 'cell' for each 8 cells )
    cat <(echo -ne "${output_header}\n") ${tmp_fname} > ${output_fname}
    rm ${tmp_fname}
    
    # return output filename
    echo ${output_fname}
}

# extract fastqc headers
function get_fastq_headers(){
    # INPUT: info about fq reads location
    # OUTPUT: fq headers for each read assocciated with SAMPLE_ID, CONDNAME
    
    # 1. get info for each sample from index file using using the pipeline cmd
    #    (01_Pipeline_Setup.py)
    # 2. extract 1st and 2nd location info and get header from these files
    # 3. save output to 'outname' file

    local outname=$1

    # remove outname if exists
    rm -rf ${outname}
    printf "%s,%s,%s,%s\n" "SAMPLE_ID" "CONDITION_NAME" "READ1_HEADER" "READ2_HEADER"> ${outname}
    
    # get all samples from Sample_labels.txt
    samples=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --samples))
    
    for ((i=0;i< ${#samples[@]} ;i+=2));
    do
	SAMPLE_ID=${samples[i]}
	CONDNAME=${samples[i+1]}

	# extract info about fastq files locations
	sample_info=($("${SETUP_PIPELINE_DIR}"/01_Pipeline_Setup.py --get_sample_info ${SAMPLE_ID}))
	READ1=${sample_info[1]}
	READ2=${sample_info[2]}
	
	# get header
	set +eu
	READ1_HEADER=$(zcat ${READ1} | head -n1)
	READ2_HEADER=$(zcat ${READ2} | head -n1)
	set -eu

	# save to file
	printf "%s,%s,%s,%s\n" "${SAMPLE_ID}" "${CONDNAME}" "${READ1_HEADER}" "${READ2_HEADER}" >> ${outname}
    done
}

function multiqc_report() {
    local local_dir=$1
    mkdir -p ${local_dir} && pushd ${local_dir}
    multiqc --no-ansi -c ${SCRIPT_DIR}/Scripts/multiqc_config.yaml ${DATASET_DIR} ${Level_UP}/03_FASTQC/ ${Level_UP}/06_CollectMetrics/ > multiqc_log.txt
    popd
}


## RUN FUNCTIONS

# find all 09a, 09b, 09c, ... directories
dedirs=$(find ../09* -iname "output*" | grep -Po '09[a-z]' | sort | uniq)
norm_methods=( "TPM" "FPKM" )

# for each 09a, 09b,...
for de_ix in $dedirs
do
    ## copy all segex files by DE index
    copy_all_segex_files ${de_ix}

    ## make combined files
    ## find all features associated with de index
    features=$(find ../${de_ix}* -iname "output*" | grep -Po "\K([a-zA-Z]*)(?=$)" | sort | uniq)
    
    for norm_method in "${norm_methods[@]}"; do
	for feature in $features; do
	    # create temporary dir
	    tmpdir=${de_ix}_${feature}_${norm_method}
	    rm -rf tmpdir && mkdir $tmpdir
	    pushd $tmpdir

	    # process files
	    link_all_segex_by_norm ${de_ix} ${norm_method} ${feature}
	    combined_fname=$(segex_combined_files $(pwd))

	    # create the output directory
	    combined_output_dir="../output/Segex_${de_ix}/Combined_${de_ix}/"
	    mkdir -p ${combined_output_dir}
	    
	    # copy output to correct dir
	    mv ${combined_fname} ${combined_output_dir}
	    popd && rm -rf $tmpdir
	done
    done

   
    ## postprocessing files for SEGEX, adding index
    pushd "./output/Segex_${de_ix}"
    cp "../../Scripts/add_segex_index.R" ./
    (set -x; Rscript add_segex_index.R ${SEGEX_EXPT_IX} ${SEGEX_FEATURE_SORT_IX})
    rm add_segex_index.R
    # TODO: rename "Upload_SEGEX_table.tsv" to dir name
    popd

    ## make summary for up/down genes
    tmpdir=${de_ix}_DE_Genes_counts
    rm -rf $tmpdir && mkdir $tmpdir && pushd ${tmpdir}

    # process files
    link_all_segex_by_ix ${de_ix}
    outname=$(updown_genes ${de_ix})

    # copy output to correct dir
    mv ${outname} ../output/
    popd &&  rm -rf $tmpdir

done

## headers for FASTQ files
echo "Extracting headers from corresponding FASTQ files"
get_fastq_headers "FASTQ_headers.csv"
mv "FASTQ_headers.csv" ./output

## multiqc reports
multiqc_report "${DATASET_LABEL}_multiqc_report"
cp -a "${DATASET_LABEL}_multiqc_report" ${VM_DIR_UCSC}/MULTIQC_REPORTS/
mv "${DATASET_LABEL}_multiqc_report" ./output

## create xlsx with summary statistics including counts/tpm/rpkm for each
## separate sample
## aggregate_counts.R <directory with steps path> <output filename>
(set -x; Rscript ./Scripts/aggregate_counts.R "../" "./output/${DATASET_LABEL}_counts_summary.xlsx")

## Create Volcano plots
echo "Creating Volcano plots"
set +eu
conda deactivate
set -eu
module load R/4.2.3

## for each segex directory with exoncollapsed files create volcano plots
cp ./Scripts/volcano_plots.R ./

## echo alldots volcano plots
for dpath in `find . -type d -name "Segex*Exon*"`; do
    dname=`basename $dpath`
    (set -x; Rscript volcano_plots.R\
		     --segex_files_path "${dpath}" \
		     --sample_labels "../00_Setup_Pipeline/Sample_Labels.txt" \
		     --comparisons "../00_Setup_Pipeline/Comparisons.txt" \
		     --output_prefix "./output/alldots_${dname}")
done

## echo 10 percent volcano plots (for fast view)
for dpath in `find . -type d -name "Segex*Exon*"`; do
    dname=`basename $dpath`
    (set -x; Rscript volcano_plots.R\
		     --segex_files_path "${dpath}" \
		     --sample_labels "../00_Setup_Pipeline/Sample_Labels.txt" \
		     --comparisons "../00_Setup_Pipeline/Comparisons.txt" \
		     --nonsignif_fraction 0.1 \
		     --output_prefix "./output/10percent_${dname}")
done

rm volcano_plots.R

module unload R
echo "All files copied. Done "


