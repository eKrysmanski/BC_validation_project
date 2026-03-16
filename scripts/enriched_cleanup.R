##Cleanup Script for Liver_bait_RNA_seq
library(tidyverse)
#Begin by reading in the data
#I think it looks tab-seperated

RNA_raw <- read.table("data/STAR_baits_liver_gene_counts.txt", 
                      sep = "\t", 
                      header = TRUE) 

#Check for missing counts
anyNA(RNA_raw[,3:8])

#Check structure
head(RNA_raw)
str(RNA_raw)

#What do all these columns mean?
#   Geneid --------------- Name of the gene
#   Chr    --------------- Looks like genomic locus
#   Start  --------------- Looks like start position
#   End    --------------- Looks like end position
#   Strand --------------- Looks like + or - strand
#   Length --------------- Looks like transcript length
#   Ex_Sx_L001.bam
#                  Ex I think is enriched 1-6   (User sample identifier)
#                  Sx I think is sample 1-6     (illumina sample index)
#                  L001 I think is the lane number on sequencer
#                  .bam I think is because these colums were spliced into 
#                  this dataframe from individual .bam files

#Let's do a few simple calculations

#Total transcripts (i.e. sequencing depth)

seq_depth <- RNA_raw %>% 
  summarise(E1_depth = sum(E1_S1_L001.bam), 
            E2_depth = sum(E2_S2_L001.bam), 
            E3_depth = sum(E3_S3_L001.bam), 
            E4_depth = sum(E4_S4_L001.bam), 
            E5_depth = sum(E5_S5_L001.bam), 
            E6_depth = sum(E6_S6_L001.bam))

seq_depth

#Normalize all the counts (Counts per Million)
# CPM = (Counts / Seq_depth) x 10^6

RNA_norm <- RNA_raw %>% 
  mutate(E1_cpm = E1_S1_L001.bam/seq_depth$E1_depth*10^6, 
         E2_cpm = E2_S2_L001.bam/seq_depth$E2_depth*10^6, 
         E3_cpm = E3_S3_L001.bam/seq_depth$E3_depth*10^6, 
         E4_cpm = E4_S4_L001.bam/seq_depth$E4_depth*10^6, 
         E5_cpm = E5_S5_L001.bam/seq_depth$E5_depth*10^6, 
         E6_cpm = E6_S6_L001.bam/seq_depth$E6_depth*10^6) %>% 
  select(Geneid, Length, ends_with("_cpm"))

head(RNA_norm)

#Some summary stats
#NOTE: I think that half of these are male and half are female; but it's unclear
# which is which based on the data that we were given; for now just pooling it
# because worse-case scenario we ignore sex differences and treat it as pooled. 

RNA_sum <- RNA_norm %>%
  rowwise() %>%
  mutate(
    CPM_mean = mean(c_across(3:8)),
    CPM_stdev = sd(c_across(3:8))
  ) %>%
  ungroup()

head(RNA_sum)

#Check my largest values...
sort(RNA_sum$mean, decreasing = TRUE) %>% 
  head(n = 12)

#From the thesis there is 1 value above 100'000, and rest are below with like 10-15
# above 10'000; this seems consistent with what I've computed here so I'm confident
# with this cleanup. 

#These genes should be a representation of ~600 that were bait-captured
#There are definitely more than 600 genes in the list;

length(RNA_raw$Geneid)

#Filter out only genes with transcript counts; but that's nearly 6000 genes. 

RNA_sub <- RNA_sum %>% 
  filter(CPM_mean > 0)

length(RNA_sub$Geneid)

rm(RNA_sub)

# There is definately some noise, and this may need to be cleaned up.

# First thing I think I need to do is figure out what genes were probed for then
# I can subset the data to only probe targets to compare with the conventional
# RNA-seq data. Will do this work in a separate script, since I need to do some 
# mapping to get from the probe target to Geneids, and that can get messy 
# "scripts/probe_targets.R"

#After more thinking; I don't think it is possible to take these gene-level annotations
# in the RNA_seq data, and map that to the individual transcripts in the probes. 
# I can just accept that there is probably noise; or figure out a way to filter out
# things that are likely contaminants? I think the transcript level stuff was pooled
# to the gene level; and I don't have a good way to figure out identity of the ~600 
# genes that were targeted without getting a list. The probe list I found is 
# like 1075 targets, which I believe includes known and predicted splice variants. 

#Save the cleaned data;
#  Save a smaller file with only genes were counts > 0
#  Save a larger file with all genes and counts even if 0

#Full file
write.csv(RNA_sum, 
          file = "data/clean/RNA_sum_CPM_full.csv")

write.csv(RNA_sub, 
          file = "data/clean/RNA_sum_CPM_grtr_zero.csv")



