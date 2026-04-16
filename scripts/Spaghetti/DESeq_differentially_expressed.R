library(readr)
library(DESeq2)
library(dplyr)
library(tibble)
library(tidyr)
library(ggvenn)

genes <- read.table("data/defensome_genes_symbol_only.tsv", sep = "\t", header = FALSE)
colnames(genes) <- "Target_Gene"
enriched_data <- as.data.frame(read_delim("data/STAR_baits_liver_gene_counts.txt", delim = "\t", skip = 1, show_col_types = FALSE))
unenriched_data <- as.data.frame(read_delim("data/STAR_library_liver_gene_counts.txt", delim = "\t", skip = 1, show_col_types = FALSE))

enriched_filt <- enriched_data %>% 
  filter(Geneid %in% genes$Target_Gene) %>% 
  column_to_rownames(var = "Geneid") %>% 
  select(-Chr,-Start,-End,-Strand,-Length) %>% 
  .[rowSums(.[, 1:6] > 3) >= 3, ] %>% 
  as.matrix()

unenriched_filt <- unenriched_data %>% 
  filter(Geneid %in% genes$Target_Gene) %>% 
  column_to_rownames(var = "Geneid") %>%  
  select(-Chr,-Start,-End,-Strand,-Length) %>% 
  .[rowSums(.[, 1:6] > 10) >= 3, ] %>% 
  as.matrix()

coldata_en <- data.frame(condition = c("male","male","male","female","female","female"))
coldata_un <- data.frame(condition = c("male","male","male","female","female","female"))

rownames(coldata_en) <- colnames(enriched_filt)
rownames(coldata_un) <- colnames(unenriched_filt)

dds_en <- DESeqDataSetFromMatrix(
  countData = enriched_filt,
  colData = coldata_en,
  design = ~ condition)
results_en<- results(DESeq(dds_en))

dds_sf_en <- estimateSizeFactors(dds_en)
enriched_counts <- as.data.frame(counts(dds_sf_en, normalized = TRUE))

dds_un <- DESeqDataSetFromMatrix(
  countData = unenriched_filt,
  colData = coldata_un,
  design = ~ condition)
results_un<- results(DESeq(dds_un))
dds_sf_un <- estimateSizeFactors(dds_un)
unenriched_counts <- as.data.frame(counts(dds_sf_un, normalized = TRUE))

sig_en <- as.data.frame(results_en[!is.na(results_en$padj) & results_en$padj < 0.05 & abs(results_en$log2FoldChange) > 1, ])
sig_un <- as.data.frame(results_un[!is.na(results_un$padj) & results_un$padj < 0.05 & abs(results_un$log2FoldChange) > 1, ])
sig_df <- genes

sig_df$sig_en <- (sig_df$Target_Gene %in% rownames(sig_en)) #189
sig_df$sig_un <- (sig_df$Target_Gene %in% rownames(sig_un)) #121 
sig_df

sig_in_en <- sig_df %>% filter(sig_df$sig_en == TRUE & sig_df$sig_un == FALSE) #87
sig_in_un <- sig_df %>% filter(sig_df$sig_en == FALSE & sig_df$sig_un == TRUE) #19

# 102 of genes are differentially expressed in both 
# 19 only seen in unenriched
# 87 only seen in enriched 

enriched <- rownames(sig_df)[sig_df$sig_en == TRUE]
unenriched <- rownames(sig_df)[sig_df$sig_un == TRUE]


venn_list <- list(
  Enriched = enriched,
  Unenriched = unenriched
)

ggvenn(venn_list)
