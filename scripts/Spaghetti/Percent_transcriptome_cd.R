library(readr)
library(DESeq2)
library(dplyr)
library(tibble)
library(tidyr)
unenriched_data <- as.data.frame(read_delim("data/STAR_library_liver_gene_counts.txt", delim = "\t", skip = 1, show_col_types = FALSE))
genes <- read.table("data/defensome_genes_symbol_only.tsv", sep = "\t", header = FALSE)
colnames(genes) <- "Target_Gene"
unenriched_filt <- unenriched_data %>% 
  filter(Geneid %in% genes$Target_Gene)
total <- sum(unenriched_data$L1_S1_L001.bam)
cd <- sum(unenriched_filt$L1_S1_L001.bam)
cd/total
