#!/usr/bin/bash

set -eu
# set -x

# ask user about removing the BAM files
make_cleanup=0
while true; do
    read -p "Do you want to clean up the SAMPLES directory? [yN]. Press Enter to N. " yn
    case $yn in
        [Yy]* ) make_cleanup=1; break;;
        [Nn]* ) break;;
	* ) break;;
    esac
done

# remove files if it is required
if [ ${make_cleanup} -eq 1 ]; then
    echo "Cleaning up SAMPLES directory"
    ./make_cleanup.sh
else
    echo "Skipping cleaning step"
fi

echo "Make compression..."

VERSION=$(ls -1 VERSIONS* | cut -d"_" -f2)
VERSION=${VERSION%.txt}
DATE=$(date +'%Y-%m-%d_%H%M')
DATASET_LABEL=$(cat ./Scripts/00_Setup_Pipeline/Pipeline_Setup.conf | grep -i dataset_label | head -n1 | cut -d "=" -f 2)
ARCHNAME="Scripts_${DATASET_LABEL}_${DATE}"

git mv Scripts ${ARCHNAME}
# zip -1 -qdgds 50m -r ${ARCHNAME}.zip ./${ARCHNAME}/* ./VERSIONS*
7z a ${ARCHNAME}.zip ./${ARCHNAME}/ ./VERSIONS*
git mv ${ARCHNAME} Scripts

echo "The Scripts directory was zipped in the ${ARCHNAME}.zip file"
