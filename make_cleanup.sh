#!/usr/bin/bash

# clean up BAM files in SAMPLES directory
find ./SAMPLES -name "*.bam" | xargs rm -rf 
