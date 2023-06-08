#!/usr/bin/env Rscript

# Author: Livio Antonielli (livio.antonielli@gmail.com)

#####################
# DADA2 workflow #2 #
#####################

# https://benjjneb.github.io/dada2/tutorial.html

# Data requirements accordin to DADA2 tutorial:
# This workflow assumes that your sequencing data meets certain criteria:
# 1.Samples have been demultiplexed, i.e. split into individual per-sample fastq files.
# 2.Non-biological nucleotides have been removed, e.g. primers, adapters, linkers, etc.
# 3.Paired-end sequencing data contain reads in matched order.

# This script takes demultiplexed, primer-free Illumina PE reads together with the read length values at which you would like to trim R1 and R2. Reads will later be trimmed, filtered, denoised, merged and ASV (Amplicon Sequence Variants) inferred. After ASV de-novo chimera filtering, a table (similar to a good old classic OTU table) will be generated. FASTA sequences representative of each ASV will be also saved and used for further purposes.  

# Sourcing libraries
#if (!requireNamespace("BiocManager"))
  #install.packages("BiocManager")
#BiocManager::install(c("dada2"), update = FALSE, ask = FALSE)
if (!require("pacman")) install.packages("pacman", repos = "http://cran.rstudio.com/")
pacman::p_load(dada2, 
               gtools,
               ggplot2,
               install = TRUE)

# Reading arguments from bash:
args <- commandArgs(trailingOnly = TRUE)
current_path <- as.character(args[1])
seq_path <- as.character(args[2]) #directory containing the fastq files after uncompressing.
truncLenR1 <- as.numeric(args[3])
truncLenR2 <- as.numeric(args[4])
filt_path <- as.character(args[5])
setwd(current_path)

# Listing the files:
fnFs <- sort(list.files(seq_path, pattern = "_R1.fastq", full.names = TRUE)) 
fnRs <- sort(list.files(seq_path, pattern = "_R2.fastq", full.names = TRUE))

# Extracting the names of the files:
sample.names <- sapply(strsplit(basename(fnFs), "_"), function(x) x[1])

# Creating a new folder and assigning new file names:
filtFs <- file.path(filt_path, paste0(sample.names, "_R1_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R2_filt.fastq.gz"))

# Filtering the reads and save them into the newly created files:
cat("\n", "==================== Filtering reads... ====================", "\n", "\n")

filt_trim_out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(truncLenR1, truncLenR2), maxN=0, maxEE=c(5, 5), truncQ=0, rm.phix=FALSE, compress=TRUE, multithread=TRUE, verbose=TRUE)

# Learning the error rates:
cat("\n", "==================== Estimating the error rates... ====================", "\n", "\n")

errF <- learnErrors(filtFs, multithread=TRUE, randomize=TRUE, MAX_CONSIST=999, verbose=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE, randomize=TRUE, MAX_CONSIST=999, verbose=TRUE)

cat("\n", "==================== Please check the error diagnostic plots in stats folder... ====================", "\n", "\n")
plotErrors(errF, nominalQ=TRUE)
ggsave("filtered_read_R1_error_plot.png", device = "png", width = 30, height = 30, units = "cm") # save a diagnostic plot for R1 reads. 

err_R2_plot <- plotErrors(errR, nominalQ=TRUE)
ggsave("filtered_read_R2_error_plot.png", device = "png", width = 30, height = 30, units = "cm") # save a diagnostic plot for R2 reads. 

# Dereplicating:
cat("\n", "==================== Dereplicating reads... ====================", "\n", "\n") 
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

names(derepFs) <- sample.names # name the derep-class objects by the sample names.
names(derepRs) <- sample.names

# Inferring the sequence variants in each sample:
cat("\n", "==================== Inferring SV in each sample... ====================", "\n", "\n")
dadaFs <- dada(derepFs, err=errF, multithread=TRUE) # apply the core sequence-variant inference algorithm to the dereplicated data.
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Merging the denoised forward and reverse reads:
cat("\n", "==================== Merging denoised reads... ====================", "\n", "\n")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = 10, maxMismatch = 1, verbose=TRUE)

# Constructin the sequence table:
cat("\n", "==================== Bulding the sequence table... ====================", "\n", "\n")
seqtab <- makeSequenceTable(mergers)

cat("\n", "==================== Inspecting the table size (samples x ASVs)... ====================", "\n", "\n")
dim(seqtab) # inspecting the table size.

cat("\n", "==================== Inspecting the distribution of sequence lengths... ====================", "\n", "\n")
table(nchar(getSequences(seqtab))) # inspecting distribution of sequence lengths.

# Removing de-novo chimeras:
cat("\n", "==================== Removing chimeras... ====================", "\n", "\n")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

cat("\n", "==================== Inspecting the table size without chimeras (samples x ASVs)... ====================", "\n", "\n")
dim(seqtab.nochim) # inspecting distribution of sequence lengths after removing chimeras.

cat("\n", "==================== Chimeric sequence relative abundance... ====================", "\n", "\n")
1-sum(seqtab.nochim)/sum(seqtab) 

# Tracking reads through the pipeline steps:
cat("\n", "==================== Please check the read counts at each step of the pipeline in stats folder... ====================", "\n", "\n")
getN <- function(x) sum(getUniques(x))
track <- cbind(filt_trim_out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("stripped", "filtered", "denoised", "merged", "non-chimeric")
rownames(track) <- sample.names

write.table(track, file = "read.counts", sep = "\t")

# Generating a FASTA file with ASV seqs:
cat("\n", "==================== Generating an ASV FASTA file... ====================", "\n", "\n")
uniquesToFasta(getUniques(seqtab.nochim), fout="seqs.fasta", ids=paste0("ASV_", seq(length(getUniques(seqtab.nochim)))))

# Exporting a count matrix with sequences as headear and the same one but with simplified header, for later...
seqtab.nochim_mod <- seqtab.nochim
colnames(seqtab.nochim_mod) <- paste0("ASV_", seq(length(colnames(seqtab.nochim))))
write.table(seqtab.nochim, file="seqtab_head_seqs.txt", col.names = NA, sep = "\t")
write.table(seqtab.nochim_mod, file="seqtab_head_names.txt", col.names = NA, sep = "\t")
