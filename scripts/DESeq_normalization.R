library(readr)
library(DESeq2)
library(dplyr)
library(tibble)
library(ggplot2)
library(TOSTER)

read_data <- read_delim("data/STAR_all_liver_gene_counts.txt", delim = "\t", skip = 1, show_col_types = FALSE)
genes <- read.table("data/defensome_genes_symbol_only.tsv", sep = "\t", header = FALSE)
colnames(genes) <- "Target_Gene"
df_filtered <- read_data[rowSums(read_data[, 7:12] > 3) >= 3, ]

filtered_data <- df_filtered[rowSums(df_filtered[, 13:18] > 10) >= 3, ]
target_data <- filtered_data %>% filter(Geneid %in% genes$Target_Gene)
view(target_data)
counts <- target_data %>%
  column_to_rownames(var = "Geneid") %>%   
  select(-Length, -Strand, -Start, -End, -Chr) %>%                     
  as.matrix()


coldata <- data.frame(
  condition = c("enriched_m","enriched_m","enriched_m","enriched_f",
                "enriched_f","enriched_f",
                "unenriched_m","unenriched_m","unenriched_m",
                "unenriched_f","unenriched_f","unenriched_f")
)

rownames(coldata) <- colnames(counts)
counts

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)
dds_filtered <- dds[rowSums(counts(dds)) > 10, ]
dds_sf <- estimateSizeFactors(dds_filtered)
sizeFactors(dds_sf)
normalized_counts <- counts(dds_sf, normalized = TRUE)


normalized_counts <- counts(dds_sf, normalized = TRUE)
normalized_counts
normalized_df <- as.data.frame(normalized_counts)
View(normalized_df)

means_table <- normalized_df %>% 
  rowwise() %>% 
  mutate(enriched_mean = mean(c_across(2:7)), 
         enriched_sd = sd(c_across(2:7)),
         unenriched_mean = mean(c_across(8:13)),
         unenriched_sd = sd(c_across(8:13)),
         enriched_m_mean = mean(c_across(2:4)),
         enriched_f_mean = mean(c_across(5:7)),
         unenriched_m_mean = mean(c_across(8:10)),
         unenriched_f_mean = mean(c_across(11:13)),
         enriched_m_sd = sd(c_across(2:4)),
         enriched_f_sd = sd(c_across(5:7)),
         unenriched_m_sd = sd(c_across(8:10)),
         unenriched_f_sd = sd(c_across(11:13)))

view(means_table)
write.csv(means_table, file = "data/clean/DESeq2_normalized_liver_combined.csv")

# Male and Female Enriched vs Unenriched
ggplot(data = means_table, aes(x = unenriched_mean, y = enriched_mean)) +
  geom_point() +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000, 100000)) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000)) +
  geom_abline(intercept = 0, slope = 1, color = "indianred", linetype = "dashed") +
  geom_smooth(method = "lm",formula = y ~ x, se = TRUE, colour ="steelblue") +
  theme_bw()

# Male Enriched vs unenriched
ggplot(data = means_table, aes(x = unenriched_m_mean, y = enriched_m_mean)) +
  geom_point() +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000, 100000)) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000)) +
  geom_abline(intercept = 0, slope = 1, color = "indianred", linetype = "dashed") +
  geom_smooth(method = "lm",formula = y ~ x, se = TRUE, colour ="steelblue") +
  theme_bw()

# Female Enriched vs Unenriched
ggplot(data = means_table, aes(x = unenriched_f_mean, y = enriched_f_mean)) +
  geom_point() +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000, 100000)) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000)) +
  geom_abline(intercept = 0, slope = 1, color = "indianred", linetype = "dashed") +
  geom_smooth(method = "lm",formula = y ~ x, se = TRUE, colour ="steelblue") +
  theme_bw()

lm_combined <- lm(log10(enriched_mean) ~ log10(unenriched_mean), data = means_table)
lm_male <- lm(log10(enriched_m_mean) ~ log10(unenriched_m_mean), data = means_table)
lm_female <- lm(log10(enriched_f_mean) ~ log10(unenriched_f_mean), data = means_table)

value_table <- summary(lm_combined)$coefficients
value_table_m <- summary(lm_male)$coefficients
value_table_f <- summary(lm_female)$coefficients

summary(lm_combined)
# Run TOST
TOST_combined_slope <- tsum_TOST(m1 = value_table[2,"Estimate"], m2 = 1, 
                           sd1 = value_table[2,"Std. Error"], sd2 = 0, 
                           n1 = nrow(means_table), n2 = nrow(means_table), 
                           eqb = 0.05, eqbound_type = "raw", 
                           var.equal = FALSE)
#Save a df
TOST_combined <- data.frame(
  d_slope = TOST_combined_slope$effsize[1,1],
  c_int_low = TOST_combined_slope$effsize[1,3],
  c_int_high = TOST_combined_slope$effsize[1,4],
  eqb_high = TOST_combined_slope$eqb[1,2],
  eqb_low = TOST_combined_slope$eqb[1,3])

# TOST for Males
TOST_male_slope <- tsum_TOST(m1 = value_table_m[2,"Estimate"], m2 = 1, 
                                 sd1 = value_table_m[2,"Std. Error"], sd2 = 0, 
                                 n1 = nrow(means_table), n2 = nrow(means_table), 
                                 eqb = 0.05, eqbound_type = "raw", 
                                 var.equal = FALSE)
#Save a male df
TOST_male <- data.frame(
  d_slope = TOST_male_slope$effsize[1,1],
  c_int_low = TOST_male_slope$effsize[1,3],
  c_int_high = TOST_male_slope$effsize[1,4],
  eqb_high = TOST_male_slope$eqb[1,2],
  eqb_low = TOST_male_slope$eqb[1,3])

# TOST for females
TOST_female_slope <- tsum_TOST(m1 = value_table_f[2,"Estimate"], m2 = 1, 
                             sd1 = value_table_f[2,"Std. Error"], sd2 = 0, 
                             n1 = nrow(means_table), n2 = nrow(means_table), 
                             eqb = 0.05, eqbound_type = "raw", 
                             var.equal = FALSE)
#Save a female df
TOST_female <- data.frame(
  d_slope = TOST_female_slope$effsize[1,1],
  c_int_low = TOST_female_slope$effsize[1,3],
  c_int_high = TOST_female_slope$effsize[1,4],
  eqb_high = TOST_female_slope$eqb[1,2],
  eqb_low = TOST_female_slope$eqb[1,3])


# Visualization 
ggplot(TOST_female, aes(x = d_slope)) +
  geom_vline(aes(xintercept = eqb_high), linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = eqb_low),  linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black") +
  geom_point(aes(y = 1), size = 5, shape = 15) +
  geom_errorbar(aes(y = 1, 
                    xmin = c_int_low,
                    xmax = c_int_high), 
                width = 0.35, 
                size = 0.75) +
  labs(x = expression(Delta~"slope"), 
       title = "Difference in Slope to Zero (Full Data)") +
  scale_y_continuous(limits = c(0.5, 1.5)) +
  scale_x_continuous(limits = c(TOST_combined$eqb_low*2, TOST_combined$eqb_high*2)) +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.line.x = element_line(color = "black", width = 0.75), 
    axis.ticks.x = element_line(color = "black", width = 0.75),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(), 
    axis.title.x = element_text(size = 15, vjust = -1),
    panel.grid.major.y = element_blank(), 
    panel.grid.minor.y = element_blank(), 
    margins = margin(t = 15, r = 15, l = 15, b = 15, unit = "pt"), 
  )
