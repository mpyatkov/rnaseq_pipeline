#!/usr/bin/bash

set -eu
# set -x
VERSION=$(ls -1 VERSIONS* | cut -d"_" -f2)
VERSION=${VERSION%.txt}
DATE=$(date +'%Y-%m-%d_%H%M')
ARCHNAME="Scripts_${VERSION}_${DATE}"
# echo ${ARCHNAME}
git mv Scripts ${ARCHNAME}
zip -1 -qdgds 50m -r ${ARCHNAME}.zip ./${ARCHNAME}/09* ./${ARCHNAME}/13* ./${ARCHNAME}/14* ./${ARCHNAME}/00*
git mv ${ARCHNAME} Scripts
echo "The Scripts directory (steps 09*,12,13,14) was zipped in the --> ${ARCHNAME}.zip <-- file"
