#!/bin/bash

## This script is used to summarize, on the basis of
## the stats file from checkm, the bins that are
## >= 50 % completeness and < 5% contamination
## then it parses the outputtted file and copies on
## the genomes_comp50_cont05 folder the genomes
## following these caracteristics
## to be used giving the name of the dataset as argument

# *********

d=${1}  #dataset_name
sp='bins_dir/checkm/storage/'

# *********

for f in ${d}*;
do

    s=`basename ${f}`
    bins_place="${d}${s}/${sp}"
    cd ${bins_place}   ##   /storage/

    sed "s/\t{.*//g" bin_stats_ext.tsv > qa.tsv
    sed "s/.*'Completeness': //g" bin_stats_ext.tsv | sed "s/,.*//g" >> qa.tsv;
    sed "s/.*'Contamination': //g" bin_stats_ext.tsv | sed "s/,.*//g" >> qa.tsv
    sed "s/.*'Genome size': //g" bin_stats_ext.tsv | sed "s/,.*//g" >> qa.tsv
    sed "s/.*'N50 (contigs)': //g" bin_stats_ext.tsv | sed "s/,.*//g" >> qa.tsv
    sed "s/.*'marker lineage': //g" bin_stats_ext.tsv | sed "s/,.*//g" >> qa.tsv;

    for i in $(seq 1 $(($(cat qa.tsv | wc -l)/6)))
    do
        sed -n ${i}~$(($(cat qa.tsv | wc -l)/6))p qa.tsv > tmp.tsv
        tr "\n" "\t" < tmp.tsv >> qa2.tsv; echo "" >> qa2.tsv
    done

    rm qa.tsv tmp.tsv
    mv qa2.tsv qa.tsv

done

# *********

cd ${d}
grep "" */bins_dir/checkm/storage/qa.tsv | sed 's,/bins_dir/checkm/storage/qa.tsv:,\t,g' > summary_checkm.txt

mkdir -p genomes_comp50_cont05

cat summary_checkm.txt | awk '$3 >= 50' | awk '$4 <= 5' | awk '{print $1"/bins_dir/"$2".fa"}' | while read i; do cp ${i} genomes_comp50_cont05/${i/\/bins_dir\//_}; done
