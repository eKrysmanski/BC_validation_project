#Genewise TOSTs
library(tidyverse)
library(TOSTER)

data <- read.csv("data/clean/DESeq2_normalized_liver_combined.csv")

data_compact <- data %>% 
  select(Geneid, ends_with("_L001.bam"), enriched_mean, enriched_sd, unenriched_mean, unenriched_sd)

####################
#Performing gene by gene TOSTs, and extracting p-values/decisions from results into data frame
#Note: the NHST results are not explicitely labeled, and are only found in a character string in the 
#       decisions part of the results; could extract it I'm sure, but it's not strictly necessary
e_bounds <- 0.1
n_val <- 6

data_tost <- data_compact %>% 
  rowwise() %>% 
  mutate(TOST_up_p = as.vector(tsum_TOST(m1 = enriched_mean, m2 = unenriched_mean, 
                             sd1 = enriched_sd, sd2 = unenriched_sd, 
                             n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$TOST["TOST Upper", "p.value"]), 
         TOST_low_p = as.vector(tsum_TOST(m1 = enriched_mean, m2 = unenriched_mean, 
                                         sd1 = enriched_sd, sd2 = unenriched_sd, 
                                         n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$TOST["TOST Lower", "p.value"]),
         NHST_p = as.vector(tsum_TOST(m1 = enriched_mean, m2 = unenriched_mean, 
                                          sd1 = enriched_sd, sd2 = unenriched_sd, 
                                          n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$TOST["t-test", "p.value"]),
         TOST_eff_size = as.vector(tsum_TOST(m1 = enriched_mean, m2 = unenriched_mean, 
                                             sd1 = enriched_sd, sd2 = unenriched_sd, 
                                             n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$effsize["Raw", "estimate"]), 
         TOST_ci_low = as.vector(tsum_TOST(m1 = enriched_mean, m2 = unenriched_mean, 
                                           sd1 = enriched_sd, sd2 = unenriched_sd, 
                                           n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$effsize["Raw", "lower.ci"]), 
         TOST_ci_up = as.vector(tsum_TOST(m1 = enriched_mean, m2 = unenriched_mean, 
                                          sd1 = enriched_sd, sd2 = unenriched_sd, 
                                          n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$effsize["Raw", "upper.ci"]))

#Check on how many have passed equivalence, and how many have failed equivalence
equivalent <- data_tost %>% 
  filter(TOST_up_p <= 0.05, 
         TOST_low_p <= 0.05)

#Welp; that's not great, looks like none of the genes are "practically equivalent"
# But that's not entirely surprising, there are differences between male and female tissue
# which may be driving up standard deviation, and as much as we try to normalize to the 
# same scale, the enriched had much higher counts and much lower sequencing depth, and
# FC differences were not very uniform between enriched and unenriched. 

#Try with the male, and female tissues independently... 

male_data <- data %>% 
  select(Geneid, 
         E1_S1_L001.bam, E2_S2_L001.bam, E3_S3_L001.bam, 
         L1_S1_L001.bam, L2_S2_L001.bam, L3_S3_L001.bam, 
         enriched_m_mean, enriched_m_sd, 
         unenriched_m_mean, unenriched_m_sd)

E_mean <-
E_sd <-
UN_mean <-
  

male_tost <- male_data %>% 
  rowwise() %>% 
  mutate(TOST_up_p = as.vector(tsum_TOST(m1 = enriched_mean, m2 = unenriched_mean, 
                                         sd1 = enriched_sd, sd2 = unenriched_sd, 
                                         n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$TOST["TOST Upper", "p.value"]), 
         TOST_low_p = as.vector(tsum_TOST(m1 = enriched_mean, m2 = unenriched_mean, 
                                          sd1 = enriched_sd, sd2 = unenriched_sd, 
                                          n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$TOST["TOST Lower", "p.value"]),
         NHST_p = as.vector(tsum_TOST(m1 = enriched_mean, m2 = unenriched_mean, 
                                      sd1 = enriched_sd, sd2 = unenriched_sd, 
                                      n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$TOST["t-test", "p.value"]),
         TOST_eff_size = as.vector(tsum_TOST(m1 = enriched_mean, m2 = unenriched_mean, 
                                             sd1 = enriched_sd, sd2 = unenriched_sd, 
                                             n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$effsize["Raw", "estimate"]), 
         TOST_ci_low = as.vector(tsum_TOST(m1 = enriched_mean, m2 = unenriched_mean, 
                                           sd1 = enriched_sd, sd2 = unenriched_sd, 
                                           n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$effsize["Raw", "lower.ci"]), 
         TOST_ci_up = as.vector(tsum_TOST(m1 = enriched_mean, m2 = unenriched_mean, 
                                          sd1 = enriched_sd, sd2 = unenriched_sd, 
                                          n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$effsize["Raw", "upper.ci"]))





#Visualize the data as...



#my "TOST-plot"



#Volcano Plot



#Regression Plot with Colour

