#Cleanup and normalize unenriched liver RNA_seq data

library(tidyverse)

#Read in the data

RNA_raw <- read.table("data/STAR_all_liver_gene_counts.txt", 
                      sep = "\t", 
                      header = TRUE)

head(RNA_raw)
str(RNA_raw)

#Okay, what the heck is in each column:
#Self-explanatory Stuff:
#"Geneid"         "Chr"            "Start"          "End"            "Strand"        
#"Length"         
#Ambiguous Stuff:
# "E1_S1_L001.bam" "E2_S2_L001.bam" "E3_S3_L001.bam" "E4_S4_L001.bam" "E5_S5_L001.bam" "E6_S6_L001.bam" 
# "L1_S1_L001.bam" "L2_S2_L001.bam" "L3_S3_L001.bam" "L4_S4_L001.bam" "L5_S5_L001.bam" "L6_S6_L001.bam"

# If I had to guess, E might be enriched 1-6? and L1-6 might be the unenriched liver?
# These are confusing naming conventions; I thought the enriched was in a seperate file
# Check if they are the same, maybe they merged them into the same dataframe

RNA_raw_en <- read.table("data/STAR_baits_liver_gene_counts.txt", 
                      sep = "\t", 
                      header = TRUE)

table(RNA_raw_en$E1_S1_L001.bam == RNA_raw$E1_S1_L001.bam)

#Okay, that is exactly what that is; I'll just subset the data since i've already taken
# the time to clean up the enriched data

RNA_raw <- RNA_raw %>% 
  select(Geneid, matches("^L[0-9]+_"))

head(RNA_raw)

#Check for NAs before moving forwards:

anyNA(RNA_raw[3:8,])

#Calculate Seq_depths for each sample

seq_depth <- RNA_raw %>% 
  summarise(L1_depth = sum(L1_S1_L001.bam), 
            L2_depth = sum(L2_S2_L001.bam), 
            L3_depth = sum(L3_S3_L001.bam), 
            L4_depth = sum(L4_S4_L001.bam), 
            L5_depth = sum(L5_S5_L001.bam), 
            L6_depth = sum(L6_S6_L001.bam))

seq_depth

#Normalize counts to CPM
# CPM = counts / seq_depth * 10^6

RNA_norm <- RNA_raw %>% 
  rowwise() %>% 
  mutate(L1_CPM = L1_S1_L001.bam/seq_depth$L1_depth*10^6, 
         L2_CPM = L2_S2_L001.bam/seq_depth$L2_depth*10^6, 
         L3_CPM = L3_S3_L001.bam/seq_depth$L3_depth*10^6, 
         L4_CPM = L4_S4_L001.bam/seq_depth$L4_depth*10^6, 
         L5_CPM = L5_S5_L001.bam/seq_depth$L5_depth*10^6, 
         L6_CPM = L6_S6_L001.bam/seq_depth$L6_depth*10^6
         ) %>% 
  select(Geneid, ends_with("_CPM"))

RNA_norm

#Simple summary statistics

RNA_sum <- RNA_norm %>% 
  rowwise() %>% 
  mutate(mean_CPM = mean(c_across(2:7)), 
         std_CPM = sd(c_across(2:7)))


head(RNA_sum)

#subset to only expressed genes
RNA_sub <- RNA_sum %>% 
  filter(mean_CPM > 0)

length(RNA_sub$Geneid)

#Save the cleaned files as the full and only genes with counts

write.csv(RNA_sum, 
          file = "data/clean/CPM_liver_unenriched_full.csv")

write.csv(RNA_sub, 
          file = "data/clean/CPM_liver_unenriched_expressed.csv")
          