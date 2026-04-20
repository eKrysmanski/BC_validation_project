library(readr)
library(DESeq2)
library(dplyr)
library(tibble)
library(tidyr)
library(ggvenn)
library(ggplot2)
library(multcomp)
library(emmeans)

# Read in data
genes <- read.table("data/defensome_genes_symbol_only.tsv", sep = "\t", header = FALSE)
colnames(genes) <- "Target_Gene"
enriched_data <- as.data.frame(read_delim("data/STAR_baits_liver_gene_counts.txt", delim = "\t", skip = 1, show_col_types = FALSE))
unenriched_data <- as.data.frame(read_delim("data/STAR_library_liver_gene_counts.txt", delim = "\t", skip = 1, show_col_types = FALSE))

# Filter data to chemical defensome genes and store as matrix for DESeq
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

# Create column data for DESeq
coldata_en <- data.frame(sex = c("male","male","male","female","female","female"))
coldata_un <- data.frame(sex = c("male","male","male","female","female","female"))

rownames(coldata_en) <- colnames(enriched_filt)
rownames(coldata_un) <- colnames(unenriched_filt)

# Run DESeq for each both enriched and unenriched libraries storing the results
dds_en <- DESeqDataSetFromMatrix(
  countData = enriched_filt,
  colData = coldata_en,
  design = ~ sex)
results_en<- results(DESeq(dds_en))

dds_un <- DESeqDataSetFromMatrix(
  countData = unenriched_filt,
  colData = coldata_un,
  design = ~ sex)
results_un<- results(DESeq(dds_un))

# From the results store biologically relevant and statistically significant DE genes
sig_en <- as.data.frame(results_en[!is.na(results_en$padj) & results_en$padj < 0.05 & abs(results_en$log2FoldChange) > 1, ])
sig_un <- as.data.frame(results_un[!is.na(results_un$padj) & results_un$padj < 0.05 & abs(results_un$log2FoldChange) > 1, ])

# Make new data frame to send DE results to 
sig_df <- genes %>% rename(Geneid = Target_Gene)

# Create TRUE FALSE columns based on if the gene is DE in either library
sig_df$sig_en <- (sig_df$Geneid %in% rownames(sig_en)) 
sig_df$sig_un <- (sig_df$Geneid %in% rownames(sig_un))  

# Create lists of genes DE in each library
enriched <- rownames(sig_df)[sig_df$sig_en == TRUE]
unenriched <- rownames(sig_df)[sig_df$sig_un == TRUE]
venn_list <- list(Enriched = enriched, Unenriched = unenriched)

# Venn Diagram
ggvenn(venn_list, fill_color = c("green4","purple"))

# Create mean columns
unenriched_data_means <- unenriched_data %>%
  mutate(male_mean = rowMeans(select(., `L1_S1_L001.bam`, `L2_S2_L001.bam`, `L3_S3_L001.bam`), na.rm = TRUE), 
         female_mean = rowMeans(select(., `L4_S4_L001.bam`, `L5_S5_L001.bam`, `L6_S6_L001.bam`), na.rm = TRUE))

# Create DF for stats and plotting that joins the DE with the average counts and creates columns with
# categories based on where genes were found to be DE
plotting<- sig_df %>%
  inner_join(unenriched_data_means, by = "Geneid") %>%
  select(Geneid, male_mean, female_mean,sig_en,sig_un)%>%
  mutate(DEG = case_when(
    sig_en == TRUE & sig_un == TRUE ~ "Both",
    sig_en == FALSE & sig_un == FALSE ~ "Neither",
    sig_en == TRUE & sig_un == FALSE ~ "Enriched",
    sig_en == FALSE & sig_un == TRUE ~ "Unenriched",
  ))

# One-way anova combining male and female averages as genes could be more expressed in both sexes
# by the library they were found to be DE 
# Using +2 allows for lowly expressed genes not found in the unenriched but found in enriched to be included 
anova_res <- aov(log10((male_mean + female_mean+2)/2) ~ DEG, data = plotting)

# Diagnostic Plots of residuals vs fitted and Q-Q residuals 
plot(anova_res, 1)
plot(anova_res, 2)


# Makes letter groups based on the pairwise comparisons and join them with plotting DF by DEG 
emm <- emmeans(anova_res, ~ DEG )
pairs(emm)
cld_res <- cld(emm, Letters = letters)
plot_labels <- merge(plotting, cld_res, by = "DEG")

# makes the label position to be %5 higher than the maximum value in the category
label_pos <- plotting %>%
  group_by(DEG) %>%
  summarise(y = max(log10((male_mean + female_mean + 2) / 2), na.rm = TRUE) *1.05)

#joins labels groups and position to same df
labelling <- left_join(label_pos, plot_labels, by = "DEG")


# Graphs box plot with points over the boxes 
ggplot(data = plotting, aes(x = DEG, y = log10((male_mean + female_mean+2)/2), fill = DEG)) +
  geom_boxplot(color = NA, alpha =0.6)+
  geom_jitter(aes(fill = DEG), shape = 21, color = "black", width = 0.2, size = 2) +
  geom_boxplot(fill=NA) +
  geom_text(data = labelling, aes(x = DEG, y = y, label = .group), vjust = 0) +
  theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.margin = margin(20, 20, 20, 20),
    legend.position = "none",  
    text = element_text(size = 13, colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title.x = element_text(vjust = -1),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.tag = element_text(face = "bold", vjust = 1),
    axis.ticks = element_line(colour = "black"),
    axis.ticks.length = unit(4, "pt"))+
  labs(x = "Library where gene is detected as DE",y = "Log10 mean expression in unenriched \n (male + female)")






