#!/bin/bash

# *********

## This script is for assembly a metagenomie sample
## from row reads. Ituses just metagahit 'cos we
## wanted to use few RAM, but maybe it better SPADES
## It is enough to change the path of the beginning
## and in case to uncomment the right extension e=

# *********

d=${1}  # dataset_name
s=${2}  # sample_name
e=".fastq.bz2" #extension of the reads

# *********
# Path declaration

root=path/to/assemblies/project
pathOut=${root}/${d}/assemblies_output/${s} #output folder
pathReads=/path/${d}/reads/${s}/${s} # Raw Reads folder and basename
py=/path/to/python #if necessary
pmega=path/to/megahit #if necessary

# *********

NOTh= #NumberOfThreads: insert reasonable number of cores for each parallel jobs according to hardware available
mkdir -p ${root}/${d}/assemblies_output/

# *********

if [ -f "${pathReads}_R1${e}" ]; then
    if [ -f "${pathReads}_UN${e}"  ]; then
        ${pmega} -1 ${pathReads}_R1${e} -2 ${pathReads}_R2${e} -r ${pathReads}_UN${e} -o ${pathOut} -t ${NOTh}
    else
        ${pmega} -1 ${pathReads}_R1${e} -2 ${pathReads}_R2${e} -o ${pathOut} -t ${NOTh}
    fi
else
    if [ -f "${pathReads}_UN${e}" ]; then
        ${pmega} -r ${pathReads}_UN${e} -o ${pathOut} -t ${NOTh}
    fi
fi

echo "No Problems for ${s}."
# *********
