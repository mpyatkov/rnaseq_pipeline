#!/usr/bin/bash -l
######################################################################################
# Kritika Karri, 08.04.2017
# Way to run script:
#Usage: ./Run_Jobs.sh
#Example:
#./Run_Jobs.sh
# The result fo this script copies all the Diffexp_v2_genebody (.txt) files from DE analysis (HTSeq Method) and copies in the current directory. These files are then used as input bu the Pearson_Script.R to generate pearson correlation plots and matrices. 

set -o errexit
set -o pipefail
set -o nounset

module load gcc/8.1.0
module load R/3.6.0

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

rm -rf *.txt
rm -rf *.pdf
rm -rf *.csv
rm -rf DiffExp_*
rm -rf output 

SCRIPT_DIR="$(pwd)"
cd ..

Level_UP=$(pwd)

mkdir ${Level_UP}/14_final_summary/output

cd ${SCRIPT_DIR}
echo "Script_Directory"
echo ${SCRIPT_DIR}
echo "Level_UP:"
echo ${Level_UP}

cp -rf  ${Level_UP}/04_TopHat_Paired_End/Job_Summary/TopHat2_Stats_BestMapped.txt  ./output/
cp -rf  ${Level_UP}/05_Read_Strandness/Job_Summary/Read_Strandness_Stats.txt ./output/
cp -rf  ${Level_UP}/06_CollectRnaSeqMetrics/Job_Summary/CollectRnaSeqMetrics_Stats.txt  ./output/
cp -rf  ${Level_UP}/07_CollectInsertSizeMetrics/Job_Summary/CollectInsertSizeMetrics_Plots.pdf ./output/
cp -rf  ${Level_UP}/08_Extract_Counts/Job_Summary/featureCounts_summary_LncRNA15k_ExonCollapsed_GTF.txt  ./output/
cp -rf ${Level_UP}/13_Correlation/Job_Summary/* ./output/
find ./output/13* -name "*.pdf" | grep -iv "combined" | xargs rm -rf
find ./output/13* \( -name "*.csv" -o -name "*.txt" \) | xargs rm -rf

# remove from 09abc/Output* all segex files
# TODO: refactor this, data should be removed in the DE step
set +eu
find ../09* -name "*_forSEGEXUpload*.txt" | grep Output | xargs -n1 rm -rf 
set -eu

copy_feature(){
    local de_index=$1
    local feature=$2
    mkdir -p output/Segex_${de_index}/Segex${de_index}_${feature}
    find ../${de_index}_DE_* -name "*SEGEX*" | grep -iv output | grep -i ${feature} | grep -i fpkm | grep -i edger | xargs -n1 -I{} cp {} ./output/Segex_${de_index}/Segex${de_index}_${feature}/
}

# copy all segex files in 
# 09a 
mkdir -p output/Segex_09a/
find ../09a_DE_* -name "*SEGEX*" | grep -iv output | grep -i fpkm | grep -i exoncollapsed | grep -i edger | xargs -n1 -I{} cp {} ./output/Segex_09a/

copy_feature 09a ExonCollapsed
copy_feature 09a IntronOnly
copy_feature 09a ExonOnly

# # 09b
mkdir output/Segex_09b
find ../09b_DE_* -name "*SEGEX*" | grep -iv output | grep -i fpkm | grep -i exoncollapsed | grep -i edger | xargs -n1 -I{} cp {} ./output/Segex_09b/

copy_feature 09b ExonCollapsed
copy_feature 09b IntronOnly
copy_feature 09b ExonOnly


# # 09c
mkdir output/Segex_09c
find ../09c_DE_* -name "*SEGEX*" | grep -iv output | grep -i fpkm | grep -i exoncoll | grep -i edger | xargs -n1 -I{} cp {} ./output/Segex_09c/

copy_feature 09c ExonCollapsed
copy_feature 09c IntronicOnly
copy_feature 09c ExonicOnly
copy_feature 09c FullGeneBody

# # 09d
mkdir output/Segex_09d
find ../09d_DE_* -name "*SEGEX*" | grep -iv output | grep -i fpkm | grep -i exoncoll | grep -i edger | xargs -n1 -I{} cp {} ./output/Segex_09d/

copy_feature 09d ExonCollapsed
copy_feature 09d FullGeneBody


# summarize up and down genes
function updown_genes() {
    local de_index=$1
    
    INPUT_UPDOWN_FILE="./gene_counts.csv"
    SAMPLES_FILE="../00_Setup_Pipeline/Sample_Labels.txt"
    CMP_FILE="../00_Setup_Pipeline/Comparisons.txt"
    
    OUTPUT_FILE="$(pwd)/${de_index}_DE_Genes_counts.csv"
    find ../"${de_index}"_DE_* -iname  "*DE_Gene_C*" | xargs -n1 -I{} cat {} | grep -v "Output_File" >> ${INPUT_UPDOWN_FILE}
    Rscript ./Scripts/updown_genes.R ${INPUT_UPDOWN_FILE} ${CMP_FILE} ${SAMPLES_FILE} ${OUTPUT_FILE} ${DATASET_LABEL}
    mv ${OUTPUT_FILE} ./output && rm ${INPUT_UPDOWN_FILE}

}

updown_genes 09a
updown_genes 09b
updown_genes 09c
updown_genes 09d

# mv *.txt ${Level_UP}/14_final_summary/output
# mv *.pdf ${Level_UP}/14_final_summary/output

echo "All files copied. Done "
