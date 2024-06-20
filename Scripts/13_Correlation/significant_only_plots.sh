#!/bin/bash
# This script removes from Combined PCA,TSNE,RPCA plots
# non-significant and all genes plots. So only Significant plots here

set -eu

function remove_nonsignificant() {

    local fname=$1
    local number_of_comparisons=$2

    ## replace Combined to SignificantOnly in filename
    output_fname=${fname/Combined/SignificantOnly}

    stop=$(((number_of_comparisons + 1) * 2))

    pdfseparate ${fname} %d_PCA.pdf

    ## iterate over *.pdf files and remove unnecessary one
    for file in *_PCA.pdf; do
        number=$(echo "${file}" | grep -oE '^[0-9]+')
        if ((number % 2 != 0 || number > stop)); then
            rm -rf "${number}_PCA.pdf"
        fi
    done

    pdfunite $(ls *_PCA.pdf | sort -n | paste -s -d " ") ${output_fname}
    rm *_PCA.pdf
}

number_of_comparisons=$(ls -l ../ | grep 09d | wc -l)

for f in `find ./Job_Summary/ -name "*Combined*" | grep -vE "Spearman|Pearson|score|Tiled"`; do
    echo "Processing $f..."
    remove_nonsignificant $f ${number_of_comparisons}
done

