library(tidyverse)
#Probe targets are in a table in the thesis; rip those out into a .csv file
# and read into R. 

probe_info <- read.csv("data/probe_info.csv", header = TRUE)

#Read in the raw enriched RNA_seq data:

RNA_raw <- read.table("data/STAR_baits_liver_gene_counts.txt", 
                      sep = "\t", 
                      header = TRUE) 

#Unfortunately; the RNA_raw contains the Chromosome Accession #s, and the 
# probe_info contains teh RefSeq transcript IDs. I will need to map one of these
# to the other, then combine the dataframes if I want to select corresponding
# Geneids for the targets. 

#Grabbed the feature table from NCBI: GRCz12tu 2025-05-13 09:55

features <- read.table("data/GCF_049306965.1_GRCz12tu_feature_table.txt", 
                       sep = "\t",
                       header = TRUE) 

features <- features %>% 
  select(chromosome, genomic_accession, product_accession, symbol, name)


#Map probe_info to features

probe_feat <- probe_info %>% 
  left_join(features, by = c("target.name" = "product_accession"))

head(probe_feat)

length(probe_info$target.name)

#Strange; why are there 1075 probe targets?

#Tried refseqR but it doesnt seem to work for these predicted transcripts
library(refseqR)

targets <- probe_info$target.name

#gene_symbols <- sapply(targets, function(x) {
#  refseq_geneSymbol(x, db = "nuccore")
#})

#str(gene_symbols)

#gene_df <- data.frame(
#  target.name = names(gene_symbols),
#  gene_symbol = sapply(gene_symbols, function(x) if (is.null(x)) NA else x),
#  stringsAsFactors = FALSE
#)  

#refseq_geneSymbol(XM_009293954.3, db = "nuccore")

#Try rentrez

library(rentrez)

