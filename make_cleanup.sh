#!/usr/bin/bash

# clean up BAM files in SAMPLES directory
set +eu
find ./SAMPLES -name "*.bam" | xargs rm -rf
set -eu
