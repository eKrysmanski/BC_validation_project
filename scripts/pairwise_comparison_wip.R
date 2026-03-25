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
         unenriched_m_mean, unenriched_m_sd) %>% 
  filter(enriched_m_mean > 0, 
         unenriched_m_mean > 0)

#Calculating stats for male data
n_val <- 3
e_bounds <- 0.1

male_tost <- male_data %>% 
  rowwise() %>% 
  mutate(TOST_up_p = as.vector(tsum_TOST(m1 = enriched_m_mean, m2 = unenriched_m_mean, 
                                         sd1 = enriched_m_sd, sd2 = unenriched_m_sd, 
                                         n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$TOST["TOST Upper", "p.value"]), 
         TOST_low_p = as.vector(tsum_TOST(m1 = enriched_m_mean, m2 = unenriched_m_mean, 
                                          sd1 = enriched_m_sd, sd2 = unenriched_m_sd, 
                                          n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$TOST["TOST Lower", "p.value"]),
         NHST_p = as.vector(tsum_TOST(m1 = enriched_m_mean, m2 = unenriched_m_mean, 
                                      sd1 = enriched_m_sd, sd2 = unenriched_m_sd, 
                                      n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$TOST["t-test", "p.value"]),
         TOST_eff_size = as.vector(tsum_TOST(m1 = enriched_m_mean, m2 = unenriched_m_mean, 
                                             sd1 = enriched_m_sd, sd2 = unenriched_m_sd, 
                                             n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$effsize["Raw", "estimate"]), 
         TOST_ci_low = as.vector(tsum_TOST(m1 = enriched_m_mean, m2 = unenriched_m_mean, 
                                           sd1 = enriched_m_sd, sd2 = unenriched_m_sd, 
                                           n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$effsize["Raw", "lower.ci"]), 
         TOST_ci_up = as.vector(tsum_TOST(m1 = enriched_m_mean, m2 = unenriched_m_mean, 
                                          sd1 = enriched_m_sd, sd2 = unenriched_m_sd, 
                                          n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$effsize["Raw", "upper.ci"]))

#Check for equivalent genes...

equivalent_male <- male_tost %>% 
  filter(TOST_up_p <= 0.05, 
         TOST_low_p <= 0.05)

#Welp; let's make sure it's doing the calculations correctly...

subset <- male_tost %>% 
  head(n = 3) %>% 
  select(!contains(".bam")) %>% 
  print()

tsum_TOST(
  m1 = male_tost[male_tost$Geneid == "txnl4b", "enriched_m_mean"] %>%  pull(),
  sd1 = male_tost[male_tost$Geneid == "txnl4b", "enriched_m_sd"] %>%  pull(),
  n1 = 3,
  m2 = male_tost[male_tost$Geneid == "txnl4b", "unenriched_m_mean"] %>%  pull(),
  sd2 = male_tost[male_tost$Geneid == "txnl4b", "unenriched_m_sd"] %>%  pull(),
  n2 = 3,
  eqb = 0.05,
  var.equal = FALSE
)

#Try female data...

fem_data <- data %>% 
  select(Geneid, 
         E1_S1_L001.bam, E2_S2_L001.bam, E3_S3_L001.bam, 
         L1_S1_L001.bam, L2_S2_L001.bam, L3_S3_L001.bam, 
         enriched_f_mean, enriched_f_sd, 
         unenriched_f_mean, unenriched_f_sd) %>% 
  filter(enriched_f_mean > 0, 
         unenriched_f_mean > 0)

#Calculating stats for male data
n_val <- 3
e_bounds <- 0.1

fem_tost <- fem_data %>% 
  rowwise() %>% 
  mutate(TOST_up_p = as.vector(tsum_TOST(m1 = enriched_f_mean, m2 = unenriched_f_mean, 
                                         sd1 = enriched_f_sd, sd2 = unenriched_f_sd, 
                                         n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$TOST["TOST Upper", "p.value"]), 
         TOST_low_p = as.vector(tsum_TOST(m1 = enriched_f_mean, m2 = unenriched_f_mean, 
                                          sd1 = enriched_f_sd, sd2 = unenriched_f_sd, 
                                          n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$TOST["TOST Lower", "p.value"]),
         NHST_p = as.vector(tsum_TOST(m1 = enriched_f_mean, m2 = unenriched_f_mean, 
                                      sd1 = enriched_f_sd, sd2 = unenriched_f_sd, 
                                      n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$TOST["t-test", "p.value"]),
         TOST_eff_size = as.vector(tsum_TOST(m1 = enriched_f_mean, m2 = unenriched_f_mean, 
                                             sd1 = enriched_f_sd, sd2 = unenriched_f_sd, 
                                             n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$effsize["Raw", "estimate"]), 
         TOST_ci_low = as.vector(tsum_TOST(m1 = enriched_f_mean, m2 = unenriched_f_mean, 
                                           sd1 = enriched_f_sd, sd2 = unenriched_f_sd, 
                                           n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$effsize["Raw", "lower.ci"]), 
         TOST_ci_up = as.vector(tsum_TOST(m1 = enriched_f_mean, m2 = unenriched_f_mean, 
                                          sd1 = enriched_f_sd, sd2 = unenriched_f_sd, 
                                          n1 = n_val, n2 = n_val, eqb = e_bounds, var.equal = FALSE)$effsize["Raw", "upper.ci"]))


#Visualize the data as...

equivalent_female <- fem_tost %>% 
  filter(TOST_up_p <= 0.05, 
         TOST_low_p <= 0.05)

#Welp; let's make sure it's doing the calculations correctly...

subset <- fem_tost %>% 
  head(n = 3) %>% 
  select(!contains(".bam")) %>% 
  print()

tsum_TOST(
  m1 = fem_tost[fem_tost$Geneid == "txnl4b", "enriched_f_mean"] %>%  pull(),
  sd1 = fem_tost[fem_tost$Geneid == "txnl4b", "enriched_f_sd"] %>%  pull(),
  n1 = 3,
  m2 = fem_tost[fem_tost$Geneid == "txnl4b", "unenriched_f_mean"] %>%  pull(),
  sd2 = fem_tost[fem_tost$Geneid == "txnl4b", "unenriched_f_sd"] %>%  pull(),
  n2 = 3,
  eqb = 0.05,
  var.equal = FALSE
)

#my "TOST-plot"

ggplot(data_tost, aes(x = TOST_eff_size, y = Geneid)) +
  geom_vline(aes(xintercept = -.1), linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = 0.1),  linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black") +
  geom_point(size = 2, alpha = 0.5) +
  geom_errorbar(aes(y = Geneid, 
                    xmin = TOST_ci_low,
                    xmax = TOST_ci_up)) +
  labs(x = expression(Delta~"slope"), 
       title = "Difference in Slope to Zero (Full Data)") +
  scale_x_continuous(limits = c(-500, 500)) +
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

#Volcano Plot



#Regression Plot with Colour

