#! /bin/sh
######################################################################################
# Kritika Karri, 08.04.2017
# Way to run script:
#Usage: ./Run_Jobs.sh
#Example:
#./Run_Jobs.sh
# The result fo this script copies all the Diffexp_v2_genebody (.txt) files from DE analysis (HTSeq Method) and copies in the current directory. These files are then used as input bu the Pearson_Script.R to generate pearson correlation plots and matrices. 
# The . of this script is two folders: 1) Pearson_All (All genes were involved to calculate the pearson) 2) Pearson_Filtered (Genes filtered by rpkm  > 1 are used for calculating the pearson values)

set -o errexit
set -o pipefail
set -o nounset

module load gcc/8.1.0
#module load R/3.1.1
#module load R/3.2
module load R/3.6.0

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

rm -rf *.txt
rm -rf *.pdf
rm -rf *.csv
rm -rf DiffExp_*
rm -rf output 

SCRIPT_DIR="$(pwd)"

echo "DATASET_DIR:"
echo ${DATASET_DIR}
echo "DATASET_LABEL:"
echo ${DATASET_LABEL}
echo "SCRIPT_DIR:"
echo ${SCRIPT_DIR}
cd ..

Level_UP=$(pwd)

mkdir ${Level_UP}/14_final_summary/output

cd ${SCRIPT_DIR}
echo "Script_Directory"
echo ${SCRIPT_DIR}
echo "Level_UP:"
echo ${Level_UP}

cp -rf  ${Level_UP}/04_TopHat_Paired_End/Job_Summary/TopHat2_Stats_BestMapped.txt  .
cp -rf  ${Level_UP}/05_Read_Strandness/Job_Summary/Read_Strandness_Stats.txt .
cp -rf  ${Level_UP}/06_CollectRnaSeqMetrics/Job_Summary/CollectRnaSeqMetrics_Stats.txt  .
cp -rf  ${Level_UP}/07_CollectInsertSizeMetrics/Job_Summary/CollectInsertSizeMetrics_Plots.pdf .
cp -rf  ${Level_UP}/08b_Extract_Counts_featureCounts/Job_Summary/featureCounts_summary_Illumina_GTF.txt  .
cp -rf  ${Level_UP}/08b_Extract_Counts_featureCounts/Job_Summary/featureCounts_summary_LncRNA_Exon_Collapsed_GTF.txt  .
cp -rf ${Level_UP}/13b_Pearson_Correlation_featureCounts/PCA* .
cp -rf ${Level_UP}/13c_Pearson_Correlation_featureCounts_lncRNA/PCA* .

mv *.txt ${Level_UP}/14_final_summary/output
mv *.pdf ${Level_UP}/14_final_summary/output
mv *.csv ${Level_UP}/14_final_summary/output



#echo 'Setup_Pipeline_DIR :' $Setup_Pipeline_DIR
cc=$(find ../ -name "[0-9][0-9]*b_DE_*_FEATURECOUNTS" | wc -l);
echo $cc
start=1
for((i=$start; i <=$cc; i++))
   { 	
       cp -rf ${Level_UP}/09b_DE_${i}_FEATURECOUNTS/Summary_Differential_Expression/SEGEX_Upload_Files/Individual/*_GeneBody_forSEGEXUpload_TPM_EdgeR_featureCounts.txt .
       cp -rf ${Level_UP}/09b_DE_${i}_FEATURECOUNTS/Summary_Differential_Expression/SEGEX_Upload_Files/Individual/*_GeneBody_forSEGEXUpload_EdgeR_featureCounts.txt .
       cp -rf ${Level_UP}/09c_DE_${i}_LNCRNA/Summary_Differential_Expression/SEGEX_Upload_Files/Individual/*_GeneBody_forSEGEXUpload_TPM_EdgeR_featureCounts.txt .
       cp -rf ${Level_UP}/09c_DE_${i}_LNCRNA/Summary_Differential_Expression/SEGEX_Upload_Files/Individual/*_GeneBody_forSEGEXUpload_EdgeR_featureCounts.txt .
   }
   mv *.txt ${Level_UP}/14_final_summary/output	

   echo "All files copied  Done "