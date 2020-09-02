######################################################################################
# Kritika Karri, 08.04.2017
# Way to run script:
#Usage: ./Run_Jobs.sh
#Example:
#./Run_Jobs.sh
# The result fo this script copies all the Diffexp_v2_genebody (.txt) files from DE analysis (FeatureCount Method) and copies in the current directory. These files are then used as input bu the Pearson_Script.R to generate pearson correlation plots and matrices.
# The output of this script is two folders: 1) Pearson_All (All genes were involved to calculate the pearson) 2) Pearson_Filtered (Genes filtered by rpkm  > 1 are used for calculating the pearson values)

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
rm -rf featureCounts

rm -rf Pearson_FPKM_Filtered 
rm -rf Pearson_All
rm -rf PCA

mkdir -p Pearson_FPKM_Filtered 
mkdir -p Pearson_All
mkdir -p PCA

SCRIPT_DIR="$(pwd)"

echo "DATASET_DIR:"
echo ${DATASET_DIR}
echo "Dataset_Label:"
echo ${DATASET_LABEL}
echo "SCRIPT_DIR:"
echo ${SCRIPT_DIR}
echo "Scripts for DiffEXP:"
#echo ${Setup_Pipeline_DIR}

cd ..

Level_UP=$(pwd)

cd ${SCRIPT_DIR}
echo "Script_Directory"
echo ${SCRIPT_DIR}
echo "Level_UP:"
echo ${Level_UP}


dd=$(find ../09b_DE_* -name "DiffEx*GeneBody*txt")
for fname in $dd
do
    # extract comparison number (09d_1_... --> 1)
    prefix=$(echo $fname | grep -Po "_\K([0-9][0-9]?)(?=_)")
    
    # copy with prefix in the begining
    cp $fname "./${prefix}_$(basename $fname)"
done


# dd=$(find ../ -name "[0-9][0-9]*b_DE_*_FEATURECOUNTS" | wc -l);
# echo $dd
# start_dd=1
# for((j=$start_dd; j <=$dd; j++))
# do
#     cp -rf ${Level_UP}/09b_DE_${j}_FEATURECOUNTS/Output_DiffExp_${j}a_featureCounts_GeneBody/DiffExp_v2_GeneBody* .
# done

Rscript pearson_script.R $DATASET_LABEL

mv Pearson_Filtered_* Pearson_FPKM_Filtered/
mv Pearson_All_* Pearson_All/
mv *_PCA_* ./PCA/
rm -rf *DiffExp*
