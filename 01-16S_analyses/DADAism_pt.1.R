#!/usr/bin/env Rscript

# Author: Livio Antonielli (livio.antonielli@gmail.com)

#####################
# DADA2 workflow #1 #
#####################

# https://benjjneb.github.io/dada2/tutorial.html

# Data requirements accordin to DADA2 tutorial:
# This workflow assumes that your sequencing data meets certain criteria:
# 1.Samples have been demultiplexed, i.e. split into individual per-sample fastq files.
# 2.Non-biological nucleotides have been removed, e.g. primers, adapters, linkers, etc.
# 3.Paired-end sequencing data contain reads in matched order.

# This script takes demultiplexed, primer-free Illumina PE reads and generate diagnostic plots to check before going further with the next filtering step.

# Sourcing libraries
#if (!requireNamespace("BiocManager"))
  #install.packages("BiocManager")
#BiocManager::install(c("dada2"), update = FALSE, ask = FALSE)
if (!require("pacman")) install.packages("pacman", repos = "http://cran.rstudio.com/")
pacman::p_load(dada2, 
               gtools,
               ggplot2,
               install = TRUE)

# Reading arguments from main bash script:
args <- commandArgs(trailingOnly = TRUE)
current_path <- as.character(args[1])
seq_path <- as.character(args[2]) #directory containing the fastq files after uncompressing.
setwd(current_path)

# Listing the files:
list.files(seq_path)

fnFs <- sort(list.files(seq_path, pattern = "_R1.fastq", full.names = TRUE)) 
fnRs <- sort(list.files(seq_path, pattern = "_R2.fastq", full.names = TRUE))

# Examining quality profiles:
plotQualityProfile(fnFs, n = 1e+06, aggregate = TRUE)
ggsave("stripped_read_R1_qual_plot.png", device = "png", width = 30, height = 30, units = "cm")
plotQualityProfile(fnRs, n = 1e+06, aggregate = TRUE)
ggsave("stripped_read_R2_qual_plot.png", device = "png", width = 30, height = 30, units = "cm")
