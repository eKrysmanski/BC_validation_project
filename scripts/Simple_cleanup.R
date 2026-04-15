#Data Cleanup before doing anything with it:

library(tidyverse)

######################      UNENRICHED DATA      ###############################

#Read in data
RNA_raw <- read.table("data/STAR_all_liver_gene_counts.txt", sep = "\t", header = TRUE)

#Filter to unenriched data, and geneIDs + Counts
RNA_raw_un <- RNA_raw %>% 
  select(Geneid, matches("^L[0-9]+_"))

#Check for missing data
anyNA(RNA_raw_un[3:8,])

#Calculate sequencing depths
seq_depth_un <- RNA_raw_un %>% 
  summarise(L1_depth = sum(L1_S1_L001.bam), 
            L2_depth = sum(L2_S2_L001.bam), 
            L3_depth = sum(L3_S3_L001.bam), 
            L4_depth = sum(L4_S4_L001.bam), 
            L5_depth = sum(L5_S5_L001.bam), 
            L6_depth = sum(L6_S6_L001.bam))

seq_depth_un

rowMeans(seq_depth_un)

#Normalize counts to CPM
# CPM = counts / seq_depth * 10^6
RNA_norm_un <- RNA_raw_un %>%
  mutate(L1_CPM = L1_S1_L001.bam/seq_depth_un$L1_depth*10^6, 
         L2_CPM = L2_S2_L001.bam/seq_depth_un$L2_depth*10^6, 
         L3_CPM = L3_S3_L001.bam/seq_depth_un$L3_depth*10^6, 
         L4_CPM = L4_S4_L001.bam/seq_depth_un$L4_depth*10^6, 
         L5_CPM = L5_S5_L001.bam/seq_depth_un$L5_depth*10^6, 
         L6_CPM = L6_S6_L001.bam/seq_depth_un$L6_depth*10^6
  ) %>% 
  select(Geneid, ends_with("_CPM"))

#Simple summary statistics

RNA_sum_un <- RNA_norm_un %>% 
  rowwise() %>% 
  mutate(CPM_mean_UN = mean(c_across(2:7)), 
         CPM_sd_UN = sd(c_across(2:7)))


head(RNA_sum_un)

#From the thesis: "Zero and low-level count data was filtered out based on counts per million 
# (CPM) to account for varying library sizes. Genes with > 0.33 CPM (~10 reads for a 30M library 
# size) from the unenriched library and > 0.5 CPM (~3 reads for a 6M library size) for the 
# enriched library in at least 3 of the 6 samples were retained.

#So, filter out low-level count data; for unenriched >0.33 CPM in at least 3/6 samples

RNA_sum_un <- RNA_sum_un %>% 
  filter(rowSums(across(ends_with("_CPM"), ~ .x > 0.33)) >= 3)

#Save the data for later use

write.csv(RNA_sum_un, file = "data/clean/CPM_liver_unenriched_full.csv")


######################        ENRICHED DATA      ###############################

#Filter to enriched geneIDs + Counts
RNA_raw_en <- RNA_raw %>% 
  select(Geneid, matches("^E[0-9]+_"))

#Check for missing counts
anyNA(RNA_raw_en[3:8,])

#Calculate Sequencing Depths
seq_depth_en <- RNA_raw_en %>% 
  summarise(E1_depth = sum(E1_S1_L001.bam), 
            E2_depth = sum(E2_S2_L001.bam), 
            E3_depth = sum(E3_S3_L001.bam), 
            E4_depth = sum(E4_S4_L001.bam), 
            E5_depth = sum(E5_S5_L001.bam), 
            E6_depth = sum(E6_S6_L001.bam))

rowMeans(seq_depth_en)

#Normalize all the counts (Counts per Million)
# CPM = (Counts / Seq_depth) x 10^6

RNA_norm_en <- RNA_raw_en %>%
  mutate(E1_cpm = E1_S1_L001.bam/seq_depth_en$E1_depth*10^6, 
         E2_cpm = E2_S2_L001.bam/seq_depth_en$E2_depth*10^6, 
         E3_cpm = E3_S3_L001.bam/seq_depth_en$E3_depth*10^6, 
         E4_cpm = E4_S4_L001.bam/seq_depth_en$E4_depth*10^6, 
         E5_cpm = E5_S5_L001.bam/seq_depth_en$E5_depth*10^6, 
         E6_cpm = E6_S6_L001.bam/seq_depth_en$E6_depth*10^6) %>% 
  select(Geneid, ends_with("_cpm"))

#Simple summary stats
RNA_sum_en <- RNA_norm_en %>% 
  rowwise() %>% 
  mutate(CPM_mean_UN = mean(c_across(2:7)), 
         CPM_sd_UN = sd(c_across(2:7)))


#Performing the filter from the
RNA_sum_en <- RNA_sum_en %>% 
  filter(rowSums(across(ends_with("_CPM"), ~ .x > 0.5)) >= 3)


#write the file
write.csv(RNA_sum, 
          file = "data/clean/CPM_liver_enriched_full.csv")

