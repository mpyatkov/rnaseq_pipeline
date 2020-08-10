#!/usr/bin/env bash

# prevent some unexpected behaviour
set -o errexit
set -o pipefail
set -o nounset

##########
# Run this file to make sure you don't have any error or warning messages
# If you have no messages, go to the next step
# Usage: ./02_Review_Pipeline_Parameters.sh
##########

# export and print all variables from Pipeline_Setup.conf 
eval $(./01_Pipeline_Setup.py --export)

# Print all exported variables from config files
echo "All exported variables:"
./01_Pipeline_Setup.py --export

# turn off hard exit when a command fails
# TODO: Before starting the pipeline, users must be sure that they have the necessary
# access to the directories, because in the next versions (set + x) the parameters
# will be removed.

set +o errexit

# checking existing VM_DIR_FASTQC 
if [[ -d ${VM_DIR_FASTQC} ]]
then
    echo "VM_DIR_FASTQC exists!"
else
    echo "Attempt to create ${VM_DIR_FASTQC} directory..."
    mkdir -p "${VM_DIR_FASTQC}"
    if [[ $? -ne 0 ]] ; then
        printf "\nWARNING: Ask somebody with rights to give you write access to $VM_DIR_FASTQC directory\n\n"
    else
        echo "Directory ${VM_DIR_FASTQC} created"
    fi
fi

# checking existing VM_DIR_UCSC
if [[ -d ${VM_DIR_UCSC} ]]
then
    echo "VM_DIR_UCSC_FASTQC exists..."
else
    echo "Need to create VM_DIR_UCSC!"
    mkdir -p ${VM_DIR_UCSC}
    if [[ $? -ne 0 ]] ; then
        printf "\nWARNING: Ask somebody with rights to give you write access to $VM_DIR_UCSC directory\n\n"
    else
        echo "Directory ${VM_DIR_UCSC} created."
    fi
fi
set -o errexit
