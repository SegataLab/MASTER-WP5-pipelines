#!/bin/bash

phylophlan_metagenomic \
    -i input_metagenomic \
    -o output_metagenomic \
    --nproc 4 \
    -n 1 \
    -d SGB.Jan21 \
    --verbose 2>&1 | tee logs/phylophlan_metagenomic.log

phylophlan_draw_metagenomic \
    -i output_metagenomic.tsv \
    -o output_heatmap \
    --map bin2meta.tsv \
    --top 20 \
    --verbose 2>&1 | tee logs/phylophlan_draw_metagenomic.log
