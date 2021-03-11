#!/usr/bin/bash

set -eu
# set -x

make_cleanup=0
while true; do
    read -p "Do you want to clean up the SAMPLES directory? [yN]. Press Enter to N. " yn
    case $yn in
        [Yy]* ) make_cleanup=1; break;;
        [Nn]* ) break;;
	* ) break;;
    esac
done

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
ARCHNAME="Scripts_${VERSION}_${DATE}"

git mv Scripts ${ARCHNAME}
zip -1 -qdgds 50m -r ${ARCHNAME}.zip ./${ARCHNAME}/09* ./${ARCHNAME}/13* ./${ARCHNAME}/12* ./${ARCHNAME}/14* ./${ARCHNAME}/00*
git mv ${ARCHNAME} Scripts

echo "The Scripts directory (steps 09*,12,13,14) was zipped in the ${ARCHNAME}.zip file"
