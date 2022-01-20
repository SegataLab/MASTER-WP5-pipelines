#!/usr/bin/env Rscript

# Author: Livio Antonielli (livio.antonielli@gmail.com)

#####################
# DADA2 workflow #3 #
#####################

# Specific for fungal ITS amplicons (UNITE 8.2 db)

# https://benjjneb.github.io/dada2/tutorial.html

# Data requirements:
# to be written.

# Sourcing libraries
#if (!requireNamespace("BiocManager"))
  #install.packages("BiocManager")
#BiocManager::install(c("dada2"), update = FALSE, ask = FALSE)
if (!require("pacman")) install.packages("pacman", repos = "http://cran.rstudio.com/")
pacman::p_load(data.table,
		   seqinr,
		   dada2,
		   plyr,
		   dplyr,
		   vegan, 
               install = TRUE)

# Reading arguments from bash:
args <- commandArgs(trailingOnly = TRUE)
current_path <- as.character(args[1])
seqtab_mod_path <- as.character(args[2])
asv.fasta_path <- as.character(args[3])
unite.fasta_path <- as.character(args[4])
setwd(current_path)

# Importing data from previous script and select ASVs according to ITSx output and replace ASV sequences with targeted extracted sequences:
cat("\n", "==================== Importing data... ====================", "\n", "\n")
seqtab_mod <- fread(seqtab_mod_path, sep = "\t", skip = 0, showProgress = TRUE, data.table = FALSE) #way faster than read.table importing
rownames(seqtab_mod) <- seqtab_mod$V1
seqtab_mod <- subset(seqtab_mod, select = -V1)
filtered_fasta_seqs <- as.character(unlist(read.fasta(file = asv.fasta_path, as.string = TRUE, forceDNAtolower = FALSE, seqonly = FALSE))) #importing fasta targeted sequences 
filtered_fasta_names <- names(read.fasta(file = asv.fasta_path, as.string = TRUE, forceDNAtolower = FALSE, seqonly = FALSE)) #same as before but for sequence names
seqtab_mod_sel <- seqtab_mod[, filtered_fasta_names] #selecting targeted ASVs on the table
colnames(seqtab_mod_sel) <- mapvalues(colnames(seqtab_mod_sel), filtered_fasta_names, filtered_fasta_seqs) #changing the names of the ASVs with their actual sequence 
seqtab_filt <- as.matrix(seqtab_mod_sel)

# Assigning taxonomy:
cat("\n", "==================== Assigning taxonomy using RDP classifier against UNITE 8... ====================", "\n", "\n")
taxa <- assignTaxonomy(seqtab_filt, unite.fasta_path, minBoot=80, multithread=TRUE, tryRC=TRUE)

# Generating a taxon table:
taxon_table <- as.data.frame(taxa)
taxon_table$sequence <- rownames(taxon_table)
rownames(taxon_table) <- filtered_fasta_names
write.table(taxon_table, file="taxon_seq_table.txt", col.names = NA, sep = "\t")

# Removing sequence information:
rownames(taxa) <- NULL

# Formatting the taxonomy old style:
 taxonomy <- as.data.frame(taxa) %>% 
 			 	transmute(taxonomy = paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = ';'))

# Building an ASV table with taxonomy formatted like old style OTU tables:
cat("\n", "==================== Building an ASV table in old style OTU table like... ====================", "\n", "\n")
asv_table <- as.data.frame(seqtab_filt)
names(asv_table) <- filtered_fasta_names
asv_table <- cbind(t(asv_table), taxonomy)

# Removing taxonomic levels with NAs:
tax_NA.list <- list(";p__NA", ";c__NA", ";o__NA", ";f__NA", ";g__NA", ";s__.* NA", ";s__NA NA")
tax_NA.string <- paste(unlist(tax_NA.list), collapse = "|")
asv_table$taxonomy <- gsub(tax_NA.string, replacement = "", x = asv_table$taxonomy)

# Filtering out ASV: taxonomy based
cat("\n", "==================== Filtering out contaminants... ====================", "\n", "\n")
exclude_these <- c("Unclassified", "unclassified", "Unknown", "unknown")
asv_table_filt <- subset(asv_table, !(asv_table$taxonomy %in% exclude_these))

# Filtering out ASV: relative abundance based 
cat("\n", "==================== Filtering out rare ASV with max relative abundance < 0.1%... ====================", "\n", "\n") 
asv_table_filt_rel <- decostand(t(asv_table_filt[, 1:length(colnames(asv_table_filt))-1]), method = "total")
apply(asv_table_filt_rel, 1, sum) #just to verify

asv_table_filt_rel_0001 <- asv_table_filt_rel %>% 
						   	as.data.frame() %>% 
						   	select_if(funs(max(.) > 0.001)) #ASVs with max relative value > 0.1%

asv_table_filt_0001 <- asv_table_filt[colnames(asv_table_filt_rel_0001), ]

# Calculating read counts and saving the FASTA and filtered ASV table for later...
cat("\n", "==================== Calculating read counts and saving the filtered ASV table... ====================", "\n", "\n") 
write.table(asv_table_filt, file="asv_table.txt", col.names = NA, sep = "\t") 
write.table(asv_table_filt_0001, file="asv_table_filtered.txt", col.names = NA, sep = "\t")
ordered_read_count <- data.frame(count=sort(colSums(subset(asv_table_filt, select= -taxonomy))))
ordered_read_count_0001 <- data.frame(count=sort(colSums(subset(asv_table_filt_0001, select= -taxonomy))))
write.table(ordered_read_count, file="asv_table_reads.counts", col.names = NA, sep = "\t")
write.table(ordered_read_count_0001, file="asv_table_filtered_reads.counts", col.names = NA, sep = "\t")
