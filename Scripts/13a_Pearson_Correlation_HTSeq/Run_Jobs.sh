######################################################################################
# Kritika Karri, 08.04.2017
# Way to run script:
#Usage: ./Run_Jobs.sh
#Example:
#./Run_Jobs.sh
# The result fo this script copies all the Diffexp_v2_genebody (.txt) files from DE analysis (HTSeq Method) and copies in the current directory. These files are then used as input bu the Pearson_Script.R to generate pearson correlation plots and matrices. 
# The output of this script is two folders: 1) Pearson_All (All genes were involved to calculate the pearson) 2) Pearson_Filtered (Genes filtered by rpkm  > 1 are used for calculating the pearson values)
module load gcc/8.1.0
#module load R/3.1.1
#module load R/3.2
module load R/3.6.0

set -o errexit
set -o pipefail
set -o nounset

# export all variables from Pipeline_Setup.conf
eval "$(../00_Setup_Pipeline/01_Pipeline_Setup.py --export)"

SCRIPT_DIR="$(pwd)"

rm -rf *.txt
rm -rf *.pdf
rm -rf *.csv
rm -rf DiffExp_*
rm -rf Pearson_All
rm -rf Pearson_Filtered

rm -rf  Pearson_All
rm -rf Pearson_Filtered
mkdir -p Pearson_All 
mkdir -p Pearson_Filtered
rm -rf PCA && mkdir -p PCA


echo "DATASET_DIR:"
echo ${DATASET_DIR}
echo "DATASET_LABEL:"
echo ${DATASET_LABEL}
echo "SCRIPT_DIR:"
echo ${SCRIPT_DIR}

cd ..

Level_UP=$(pwd)

cd ${SCRIPT_DIR}
echo "Script_Directory"
echo ${SCRIPT_DIR}
echo "Level_UP:"
echo ${Level_UP}

#echo 'Setup_Pipeline_DIR :' $Setup_Pipeline_DIR
cc=$(find ../ -name "[0-9][0-9]*a_DE_*_HTSEQ" | wc -l);
echo $cc
start=1
for((i=$start; i <=$cc; i++))
do 
    cp -rf ${Level_UP}/09a_DE_${i}_HTSEQ/Output_DiffExp_${i}a_HTSeq_GeneBody/DiffExp_v2_GeneBody* .
done
Rscript pearson_script.R

mv  Pearson_Filtered_* Pearson_Filtered/
mv  Pearson_All_* Pearson_All/ 
mv PCA_* ./PCA/
echo "HTSeq Pearson Analysis Done "

echo "#######################################################"

rm -rf DiffEx
