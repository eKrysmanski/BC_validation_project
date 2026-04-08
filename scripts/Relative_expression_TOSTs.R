#Looking into relative expression to a reference for comparing gene by gene
#Author = Evan C. Krysmanski
##Reading in data##  
library(tidyverse)
library(ggplot2)
library(ctrlGene)
library(TOSTER)

data <- read.csv(file = "data/clean/probe_targets_full.csv", header = TRUE)

data_sort_en <- data %>%
  select(Geneid, ends_with("_cpm", ignore.case = FALSE), CPM_mean_E, CPM_sd_E) %>% 
  mutate(CV = CPM_sd_E/CPM_mean_E) %>% 
  arrange(CV)

data_sort_un <- data %>% 
  select(Geneid, ends_with("_CPM", ignore.case = FALSE), CPM_mean_UN, CPM_sd_UN) %>% 
  mutate(CV = CPM_sd_UN/CPM_mean_UN) %>% 
  arrange(CV)

#Grab the genes with the lowest variation between the enriched and unenriched to 
# use as candidates for a reference gene.. 

candidates_full_en <- data_sort_en %>%
  mutate(rank_full_en = rank(CV, ties.method = "first"))

candidates_full_un <- data_sort_un %>%
  mutate(rank_full_un = rank(CV, ties.method = "first"))

candidates_full <- inner_join(candidates_full_en, candidates_full_un, by = "Geneid") %>%
  filter(rank_full_en <= 50,
         rank_full_un <= 50)

#Formating the dataframe for CtrlGene
candidates_full_en_CG <- data_sort_en %>%
  filter(Geneid %in% candidates_full$Geneid) %>% 
  select(Geneid, ends_with("_CPM")) %>% 
  column_to_rownames(var = "Geneid") %>% 
  t()

#Repeat with the unenriched data
candidates_full_un_CG <- data_sort_un %>%
  filter(Geneid %in% candidates_full$Geneid) %>% 
  select(Geneid, ends_with("_CPM")) %>% 
  column_to_rownames(var = "Geneid") %>% 
  t()

#Get M-values
results_full_en <- geNorm(candidates_full_en_CG, ctVal = FALSE)

results_full_un <- geNorm(candidates_full_un_CG, ctVal = FALSE)
#Simple index to decide on best choice of reference

rank_full_un <- data.frame(results_full_un) %>%
  separate_rows(Genes, sep = "-") %>% 
  mutate(rank_full_un = c((nrow(.)-1):1,1))

rank_full_en <- data.frame(results_full_en) %>%
  separate_rows(Genes, sep = "-") %>% 
  mutate(rank_full_en = c((nrow(.)-1):1,1))

rankings_full <- left_join(rank_full_en, rank_full_un, by = "Genes") %>% 
  mutate(total = rank_full_en + rank_full_un) %>% 
  rename(Gene_ID = "Genes", 
         M_en = "Avg.M.x", Rank_en = "rank_full_en",  
         M_un = "Avg.M.y", Rank_un = "rank_full_un", 
         Combined_rank = "total") %>% 
  print()
#Looks like best reference gene here is txn2; lowest M-value
################################################################################

#Create a dataframe with the reference gene counts per sample
ref_en <- data_sort_en %>%
  filter(Geneid == "txn2") %>%
  select(ends_with("_cpm", ignore.case = FALSE))

#Calculate log2(FC) for each gene relative to the appropriate reference count

en_fc <- data_sort_en %>%
  rowwise() %>%
  mutate(FC_E1 = (E1_cpm / ref_en$E1_cpm),
         FC_E2 = (E2_cpm / ref_en$E2_cpm),
         FC_E3 = (E3_cpm / ref_en$E3_cpm),
         FC_E4 = (E4_cpm / ref_en$E4_cpm),
         FC_E5 = (E5_cpm / ref_en$E5_cpm),
         FC_E6 = (E6_cpm / ref_en$E6_cpm)) %>%
  mutate(mean_FC_e = mean(c_across(starts_with("FC_E"))), 
         sd_FC_e = sd(c_across(starts_with("FC_E")))) %>% 
  ungroup()

#Do the same for the unenriched data.... 

ref_un <- data_sort_un %>% 
  filter(Geneid == "txn2") %>% 
  select(ends_with("_CPM", ignore.case = FALSE))

un_fc <- data_sort_un %>% 
  rowwise() %>% 
  mutate(FC_UN1 = (L1_CPM / ref_un$L1_CPM),
         FC_UN2 = (L2_CPM / ref_un$L2_CPM),
         FC_UN3 = (L3_CPM / ref_un$L3_CPM),
         FC_UN4 = (L4_CPM / ref_un$L4_CPM),
         FC_UN5 = (L5_CPM / ref_un$L5_CPM),
         FC_UN6 = (L6_CPM / ref_un$L6_CPM)) %>%
  mutate(mean_FC_un = mean(c_across(starts_with("FC_UN"))), 
         sd_FC_un = sd(c_across(starts_with("FC_UN")))) %>% 
  ungroup() 

#Join the dataframes

FC_data <- inner_join(en_fc, un_fc, by = "Geneid") %>% 
  select(Geneid, mean_FC_e, sd_FC_e, mean_FC_un, sd_FC_un) %>% 
  filter(is.finite(mean_FC_e), 
         is.finite(mean_FC_un))

ggplot(data = FC_data, aes(x = mean_FC_un, y = mean_FC_e)) +
  geom_point(alpha = 0.5) +
  scale_y_log10(breaks = 10^(-1:3),
                labels = scales::label_parse()(paste0("10^", -1:3))) +
  scale_x_log10(breaks = 10^(-1:3),
                labels = scales::label_parse()(paste0("10^", -1:3))) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "purple") +
  labs(x = "Mean(CPM/txn2)\n[Unenriched]", 
       y = "Mean(CPM/txn2)\n[Enriched]") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    margins = margin(r = 15, l = 15, t = 15, b = 15),
    axis.title.x = element_text(vjust = .5, size = 12),
    axis.title.y = element_text(hjust = .5, size = 12),
    axis.text = element_text(size = 10, colour = "black", font = "arial"),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.line.y.right = element_blank(),
    axis.line.x.top = element_blank()
  )

#Linear Model
full_lm <- lm(formula = log10(mean_FC_e) ~ log10(mean_FC_un), data = FC_data)
#plot(full_lm) #not the most normal thing in the world; but alright and lm's are pretty robust against escape from normality
full_lm_res <- summary(full_lm)

#Intercepy TOST
full_int_tost <- tsum_TOST(m1 = full_lm_res$coefficients["(Intercept)","Estimate"], 
                          sd1 = full_lm_res$coefficients["(Intercept)", "Std. Error"], 
                          n1 = nrow(FC_data), 
                          eqb = 0.1, eqbound_type = "raw", var.equal = FALSE)

#Slope TOST
full_slope_tost <- tsum_TOST(m1 = full_lm_res$coefficients["log10(mean_FC_un)","Estimate"] - 1, 
                            sd1 = full_lm_res$coefficients["log10(mean_FC_un)", "Std. Error"]*sqrt(nrow(FC_data)), 
                            n1 = nrow(FC_data),
                            eqb = 0.1, eqbound_type = "raw", var.equal = FALSE)

full_int_tost   #Intercept is clearly different from 0, and  practically equivalent to 0 with 10% bounds
full_slope_tost #Intercept is clearly different from 0, and practically equivalent to 1 within 10% bounds


#TOST
e_bounds <- 3
n_val <- 6

TOST_full <- FC_data %>% 
  filter(mean_FC_e > 0, 
         mean_FC_un > 0, 
         Geneid != "txn2") %>% 
  rowwise() %>% 
  mutate(TOST_up_p = as.vector(tsum_TOST(m1 = mean_FC_e, m2 = mean_FC_un, 
                                         sd1 = sd_FC_e, sd2 = sd_FC_un, 
                                         n1 = n_val, n2 = n_val, eqb = e_bounds, 
                                         eqbound_type = "raw", var.equal = FALSE)$TOST["TOST Upper", "p.value"]), 
         TOST_low_p = as.vector(tsum_TOST(m1 = mean_FC_e, m2 = mean_FC_un, 
                                          sd1 = sd_FC_e, sd2 = sd_FC_un, 
                                          n1 = n_val, n2 = n_val, eqb = e_bounds,
                                          eqbound_type = "raw", var.equal = FALSE)$TOST["TOST Lower", "p.value"]),
         NHST_p = as.vector(tsum_TOST(m1 = mean_FC_e, m2 = mean_FC_un, 
                                      sd1 = sd_FC_e, sd2 = sd_FC_un, 
                                      n1 = n_val, n2 = n_val, eqb = e_bounds, 
                                      eqbound_type = "raw", var.equal = FALSE)$TOST["t-test", "p.value"]),
         TOST_eff_size = as.vector(tsum_TOST(m1 = mean_FC_e, m2 = mean_FC_un, 
                                             sd1 = sd_FC_e, sd2 = sd_FC_un, 
                                             n1 = n_val, n2 = n_val, eqb = e_bounds,
                                             eqbound_type = "raw", var.equal = FALSE)$effsize["Raw", "estimate"]), 
         TOST_ci_low = as.vector(tsum_TOST(m1 = mean_FC_e, m2 = mean_FC_un, 
                                           sd1 = sd_FC_e, sd2 = sd_FC_un, 
                                           n1 = n_val, n2 = n_val, eqb = e_bounds, 
                                           eqbound_type = "raw", var.equal = FALSE)$effsize["Raw", "lower.ci"]), 
         TOST_ci_up = as.vector(tsum_TOST(m1 = mean_FC_e, m2 = mean_FC_un, 
                                          sd1 = sd_FC_e, sd2 = sd_FC_un, 
                                          n1 = n_val, n2 = n_val, eqb = e_bounds, 
                                          eqbound_type = "raw", var.equal = FALSE)$effsize["Raw", "upper.ci"]))

TOST_full <- TOST_full %>% 
  mutate(TOST_up_p.adjust = p.adjust(TOST_up_p, method = "BH", n = length(TOST_full$TOST_up_p)), 
         TOST_low_p.adjust = p.adjust(TOST_low_p, method = "BH", n = length(TOST_full$TOST_low_p)), 
         NHST_p.adjust = p.adjust(NHST_p, method = "BH", n = length(TOST_full$NHST_p)))

equivalent_full <- TOST_full %>% 
  filter(TOST_low_p < 0.05, 
         TOST_up_p < 0.05)

equivalent_full_adjust <- TOST_full %>% 
  filter(TOST_low_p.adjust < 0.05, 
         TOST_up_p.adjust < 0.05)

nrow(equivalent_full)
nrow(equivalent_full_adjust)

nrow(equivalent_full) - nrow(equivalent_full_adjust)

#Plot the 'equivalent'

ggplot(equivalent_full %>% filter(NHST_p > 0.05), aes(x = TOST_eff_size, y = Geneid)) +
  geom_vline(aes(xintercept = -1*e_bounds), linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = e_bounds),  linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black") +
  geom_point(size = 1.5, 
             alpha = 0.5,
             shape = 22, 
             colour = "purple4", fill = "purple3") +
  geom_errorbar(aes(y = Geneid, xmin = TOST_ci_low, xmax = TOST_ci_up), 
                colour = "purple2", 
                width = 0, 
                linewidth = 1, 
                alpha = 0.5) +
  labs(x = expression(Delta~"FC (trancript/txn2)"), 
       title = "Difference in Slope to Zero (Full Data)") +
  scale_x_continuous(limits = c(-5, 5)) +
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

#plot everything

ggplot(TOST_full, aes(x = TOST_eff_size, y = Geneid)) +
  geom_vline(aes(xintercept = -1*e_bounds), linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = e_bounds),  linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black") +
  geom_point(size = 1.5, 
             alpha = 0.6,
             shape = 22, 
             colour = "purple4", fill = "purple3") +
  geom_errorbar(aes(y = Geneid, xmin = TOST_ci_low, xmax = TOST_ci_up),
                colour = "purple2", 
                width = 0, 
                linewidth = 1, 
                alpha = 0.5) +
  labs(x = expression(Delta~"FC (trancript/txn2)"), 
       title = "Difference in Relative Expression (Full Data)") +
  scale_x_continuous(limits = c(-20, 20)) +
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

#Quite a few closer to eqbounds than I expected; a couple FC differences are even clearly different from zero. 

#I imagine that sex differences are potentially inflating sd; so let's try to look at male and female results
# seperately now.. 

############################       MALE      ###################################

data_male <- read.csv(file = "data/clean/probe_targets_male.csv", header = TRUE)

#Check if there is a better reference for male than txn2; txn2 was least variable 
# overall, but with just the males there might be a better option

data_male_en <- data_male %>%
  select(Geneid, starts_with("E", ignore.case = FALSE), CPM_mean_E13) %>% 
  rowwise() %>% 
  mutate(CPM_sd_E13 = sd(c_across(starts_with("E"))), 
         CV = CPM_sd_E13/CPM_mean_E13) %>% 
  ungroup() %>% 
  arrange(CV)

data_male_un <- data_male %>% 
  select(Geneid, starts_with("L", ignore.case = FALSE), CPM_mean_UN13) %>% 
  rowwise() %>% 
  mutate(CPM_sd_UN13 = sd(c_across(starts_with("L"))), 
         CV = CPM_sd_UN13/CPM_mean_UN13) %>% 
  ungroup() %>% 
  arrange(CV)

#Grab the genes with the lowest variation between the enriched and unenriched to 
# use as candidates for a reference gene.. 

candidates_male_en <- data_male_en %>%
  mutate(rank_male_en = rank(CV, ties.method = "first"))

candidates_male_un <- data_male_un %>%
  mutate(rank_male_un = rank(CV, ties.method = "first"))

candidates_male <- inner_join(candidates_male_en, candidates_male_un, by = "Geneid") %>%
  filter(rank_male_en <= 20,
         rank_male_un <= 20)

#Formating the dataframe for CtrlGene
candidates_male_en_CG <- data_male_en %>%
  filter(Geneid %in% candidates_male$Geneid) %>% 
  select(Geneid, ends_with("_CPM")) %>% 
  column_to_rownames(var = "Geneid") %>% 
  t()

#Repeat with the unenriched data
candidates_male_un_CG <- data_male_un %>%
  filter(Geneid %in% candidates_male$Geneid) %>% 
  select(Geneid, ends_with("_CPM")) %>% 
  column_to_rownames(var = "Geneid") %>% 
  t()

#Get M-values
results_male_en <- geNorm(candidates_male_en_CG, ctVal = FALSE)

results_male_un <- geNorm(candidates_male_un_CG, ctVal = FALSE)

#Simple index to decide on best choice of reference
rank_male_un <- data.frame(results_male_un) %>%
  separate_rows(Genes, sep = "-") %>% 
  mutate(rank_male_un = c(nrow(.):1))

rank_male_en <- data.frame(results_male_en) %>%
  separate_rows(Genes, sep = "-") %>% 
  mutate(rank_male_en = c(nrow(.):1))

rankings_male <- left_join(rank_male_en, rank_male_un, by = "Genes") %>% 
  mutate(total = rank_male_en + rank_male_un) %>% 
  print()

#Looks like m msgst1.2 the best choice for reference for males; let's
# go ahead and use that... 

ref_male <- data_male %>%
  filter(Geneid == "mgst1.2") %>%
  select(ends_with("_cpm"))

#Calculate FC for each gene relative to the appropriate reference count

en_fc_male <- data_male_en %>%
  rowwise() %>%
  mutate(FC_E1 = (E1_cpm / ref_male$E1_cpm),
         FC_E2 = (E2_cpm / ref_male$E2_cpm),
         FC_E3 = (E3_cpm / ref_male$E3_cpm)) %>%
  mutate(mean_FC_e = mean(c_across(starts_with("FC_E"))), 
         sd_FC_e = sd(c_across(starts_with("FC_E")))) %>% 
  ungroup()

#Do the same for the unenriched data....

un_fc_male <- data_male_un %>% 
  rowwise() %>% 
  mutate(FC_UN1 = (L1_CPM / ref_male$L1_CPM),
         FC_UN2 = (L2_CPM / ref_male$L2_CPM),
         FC_UN3 = (L3_CPM / ref_male$L3_CPM)) %>%
  mutate(mean_FC_un = mean(c_across(starts_with("FC_UN"))), 
         sd_FC_un = sd(c_across(starts_with("FC_UN")))) %>% 
  ungroup() 

#Join dataframes
FC_male <- inner_join(en_fc_male, un_fc_male, by = "Geneid") %>% 
  select(Geneid, mean_FC_e, sd_FC_e, mean_FC_un, sd_FC_un) %>% 
  filter(is.finite(mean_FC_e), 
         is.finite(mean_FC_un))

#plot to take a look; looks pretty good

ggplot(data = FC_male, aes(x = mean_FC_un, y = mean_FC_e)) +
  geom_point(alpha = 0.5) +
  scale_y_log10(breaks = 10^(-4:2),
                labels = scales::label_parse()(paste0("10^", -4:2))) +
  scale_x_log10(breaks = 10^(-4:2),
                labels = scales::label_parse()(paste0("10^", -4:2))) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "steelblue") +
  labs(x = "Mean(CPM/mgst1.2)\n[Unenriched]", 
       y = "Mean(CPM/mgst1.2)\n[Enriched]") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    margins = margin(r = 15, l = 15, t = 15, b = 15),
    axis.title.x = element_text(vjust = .5, size = 12),
    axis.title.y = element_text(hjust = .5, size = 12),
    axis.text = element_text(size = 10, colour = "black", font = "arial"),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.line.y.right = element_blank(),
    axis.line.x.top = element_blank()
  )

#TOST on line for intercept of 0, and slope of 1; claim equivalence for relative expression patterns

#Linear Model
male_lm <- lm(formula = log10(mean_FC_e) ~ log10(mean_FC_un), data = FC_male)
#plot(male_lm) #not the most normal thing in the world; but alright and lm's are pretty robust against escape from normality
male_lm_res <- summary(male_lm)

#Intercept TOST
male_int_tost <- tsum_TOST(m1 = male_lm_res$coefficients["(Intercept)","Estimate"], 
                       sd1 = male_lm_res$coefficients["(Intercept)", "Std. Error"], 
                       n1 = nrow(FC_male), 
                       eqb = 0.1, eqbound_type = "raw", var.equal = FALSE)

#Slope TOST
male_slope_tost <- tsum_TOST(m1 = male_lm_res$coefficients["log10(mean_FC_un)","Estimate"] - 1, 
                             sd1 = male_lm_res$coefficients["log10(mean_FC_un)", "Std. Error"], 
                             n1 = nrow(FC_male), 
                             mu = 1,
                             eqb = 0.1, eqbound_type = "raw", var.equal = FALSE)


male_int_tost #Intercept is clearly different from zero, but also not practically equivalent to zero with 10% bounds
male_slope_tost #Slope is clearly different from 1; but practically equivalent within 10% eqbounds


#Now perform equivalence tests gene by gene;
n_val <- 3
e_bounds <- 3

TOST_male <- FC_male %>% 
  filter(mean_FC_e > 0, 
         mean_FC_un > 0, 
         Geneid != "mgst1.2") %>% 
  rowwise() %>% 
  mutate(TOST_up_p = as.vector(tsum_TOST(m1 = mean_FC_e, m2 = mean_FC_un, 
                                         sd1 = sd_FC_e, sd2 = sd_FC_un, 
                                         n1 = n_val, n2 = n_val, eqb = e_bounds, 
                                         eqbound_type = "raw", var.equal = FALSE,)$TOST["TOST Upper", "p.value"]), 
         TOST_low_p = as.vector(tsum_TOST(m1 = mean_FC_e, m2 = mean_FC_un, 
                                          sd1 = sd_FC_e, sd2 = sd_FC_un, 
                                          n1 = n_val, n2 = n_val, eqb = e_bounds,
                                          eqbound_type = "raw", var.equal = FALSE)$TOST["TOST Lower", "p.value"]),
         NHST_p = as.vector(tsum_TOST(m1 = mean_FC_e, m2 = mean_FC_un, 
                                      sd1 = sd_FC_e, sd2 = sd_FC_un, 
                                      n1 = n_val, n2 = n_val, eqb = e_bounds, 
                                      eqbound_type = "raw", var.equal = FALSE)$TOST["t-test", "p.value"]),
         TOST_eff_size = as.vector(tsum_TOST(m1 = mean_FC_e, m2 = mean_FC_un, 
                                             sd1 = sd_FC_e, sd2 = sd_FC_un, 
                                             n1 = n_val, n2 = n_val, eqb = e_bounds,
                                             eqbound_type = "raw", var.equal = FALSE)$effsize["Raw", "estimate"]), 
         TOST_ci_low = as.vector(tsum_TOST(m1 = mean_FC_e, m2 = mean_FC_un, 
                                           sd1 = sd_FC_e, sd2 = sd_FC_un, 
                                           n1 = n_val, n2 = n_val, eqb = e_bounds, 
                                           eqbound_type = "raw", var.equal = FALSE)$effsize["Raw", "lower.ci"]), 
         TOST_ci_up = as.vector(tsum_TOST(m1 = mean_FC_e, m2 = mean_FC_un, 
                                          sd1 = sd_FC_e, sd2 = sd_FC_un, 
                                          n1 = n_val, n2 = n_val, eqb = e_bounds, 
                                          eqbound_type = "raw", var.equal = FALSE)$effsize["Raw", "upper.ci"]))
#adjust p-vals
TOST_male <- TOST_male %>% 
  mutate(TOST_up_p.adjust = p.adjust(TOST_up_p, method = "BH", n = length(TOST_male$TOST_up_p)), 
         TOST_low_p.adjust = p.adjust(TOST_low_p, method = "BH", n = length(TOST_male$TOST_low_p)), 
         NHST_p.adjust = p.adjust(NHST_p, method = "BH", n = length(TOST_male$NHST_p)))

#Subset to equivalence
equivalent_male <- TOST_male %>% 
  filter(TOST_low_p < 0.05, 
         TOST_up_p < 0.05)

nrow(equivalent_male)

#Subset to equivalence with adjusted p-vals
equivalent_male_adjust <- TOST_male %>% 
  filter(TOST_low_p.adjust < 0.05, 
         TOST_up_p.adjust < 0.05)

nrow(equivalent_male_adjust)

#How many were lost (54 apparently)
nrow(equivalent_male) - nrow(equivalent_male_adjust)

#Simple plotting of the equivalent values
ggplot(equivalent_male %>% filter(NHST_p > 0.05), aes(x = TOST_eff_size, y = Geneid)) +
  geom_vline(aes(xintercept = -1*e_bounds), linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = e_bounds),  linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black") +
  geom_point(size = 1.5, 
             alpha = 0.7, 
             colour = "navy", fill = "darkblue", 
             shape = 22) +
  geom_errorbar(aes(y = Geneid, 
                    xmin = TOST_ci_low,
                    xmax = TOST_ci_up), 
                linewidth = 1, 
                alpha = 0.5, 
                colour = "steelblue") +
  labs(x = expression(Delta~"FC (trancript/mgst1.2)"), 
       title = "Difference in Relative Expression (Male Data)") +
  scale_x_continuous(limits = c(-5, 5)) +
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

#Look at unequivalent/full males

ggplot(TOST_male, aes(x = TOST_eff_size, y = Geneid)) +
  geom_vline(aes(xintercept = -1*e_bounds), linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = e_bounds),  linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black") +
  geom_point(size = 1.5, 
             alpha = 0.7, 
             colour = "navy", fill = "darkblue", 
             shape = 22) +
  coord_cartesian(clip = "off") +
  geom_errorbar(aes(y = Geneid, xmin = TOST_ci_low, xmax = TOST_ci_up), 
                linewidth = 1, 
                width = 0,
                alpha = 0.5, 
                colour = "steelblue") +
  labs(x = expression(Delta~"FC (transcript/mgst1.2)"), 
       title = "Difference in Relative Expresion (Male Data)") +
  scale_x_continuous(limits = c(-25, 20)) +
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

#honestly not that bad;

#######################        FEMALE         ##################################

data_fem <- read.csv(file = "data/clean/probe_targets_female.csv", header = TRUE)

#Check if there is a better reference for male than txn2; txn2 was least variable 
# overall, but with just the males there might be a better option

data_fem_en <- data_fem %>%
  select(Geneid, starts_with("E", ignore.case = FALSE), CPM_mean_E46) %>% 
  rowwise() %>% 
  mutate(CPM_sd_E46 = sd(c_across(starts_with("E"))), 
         CV = CPM_sd_E46/CPM_mean_E46) %>% 
  ungroup() %>% 
  arrange(CV)

data_fem_un <- data_fem %>% 
  select(Geneid, starts_with("L", ignore.case = FALSE), CPM_mean_UN46) %>% 
  rowwise() %>% 
  mutate(CPM_sd_UN46 = sd(c_across(starts_with("L"))), 
         CV = CPM_sd_UN46/CPM_mean_UN46) %>% 
  ungroup() %>% 
  arrange(CV)


#Grab the genes with the lowest variation between the enriched and unenriched to 
# use as candidates for a reference gene.. 

candidates_fem_en <- data_fem_en %>%
  mutate(rank_fem_en = rank(CV, ties.method = "first"))

candidates_fem_un <- data_fem_un %>%
  mutate(rank_fem_un = rank(CV, ties.method = "first"))

candidates_fem <- inner_join(candidates_fem_en, candidates_fem_un, by = "Geneid") %>%
  filter(rank_fem_en <= 50,
         rank_fem_un <= 50)

#Formating the dataframe for CtrlGene
candidates_fem_en_CG <- data_fem_en %>%
  filter(Geneid %in% candidates_fem$Geneid) %>% 
  select(Geneid, ends_with("_CPM")) %>% 
  column_to_rownames(var = "Geneid") %>% 
  t()

#Repeat with the unenriched data
candidates_fem_un_CG <- data_fem_un %>%
  filter(Geneid %in% candidates_fem$Geneid) %>% 
  select(Geneid, ends_with("_CPM")) %>% 
  column_to_rownames(var = "Geneid") %>% 
  t()

#Get M-values
results_fem_en <- geNorm(candidates_fem_en_CG, ctVal = FALSE)

results_fem_un <- geNorm(candidates_fem_un_CG, ctVal = FALSE)

#Simple index to decide on best choice of reference

rank_fem_un <- data.frame(results_fem_un) %>%
  separate_rows(Genes, sep = "-") %>% 
  mutate(rank_fem_un = c(nrow(.):1))

rank_fem_en <- data.frame(results_fem_en) %>%
  separate_rows(Genes, sep = "-") %>% 
  mutate(rank_fem_en = c(nrow(.):1))

rankings_fem <- left_join(rank_fem_en, rank_fem_un, by = "Genes") %>% 
  mutate(total = rank_fem_en + rank_fem_un) %>% 
  print()

#Looks like sod2 is the best choice for reference for males; let's
# go ahead and use that... 

ref_fem <- data_fem %>%
  filter(Geneid == "sod2") %>%
  select(ends_with("_cpm"))

#Calculate FC for each gene relative to the appropriate reference count

en_fc_fem <- data_fem_en %>%
  rowwise() %>%
  mutate(FC_E4 = (E4_cpm / ref_fem$E4_cpm),
         FC_E5 = (E5_cpm / ref_fem$E5_cpm),
         FC_E6 = (E6_cpm / ref_fem$E6_cpm)) %>%
  mutate(mean_FC_e = mean(c_across(starts_with("FC_E"))), 
         sd_FC_e = sd(c_across(starts_with("FC_E")))) %>% 
  ungroup()

#Do the same for the unenriched data....

un_fc_fem <- data_fem_un %>% 
  rowwise() %>% 
  mutate(FC_UN4 = (L4_CPM / ref_fem$L4_CPM),
         FC_UN5 = (L5_CPM / ref_fem$L5_CPM),
         FC_UN6 = (L6_CPM / ref_fem$L6_CPM)) %>%
  mutate(mean_FC_un = mean(c_across(starts_with("FC_UN"))), 
         sd_FC_un = sd(c_across(starts_with("FC_UN")))) %>% 
  ungroup() 

#Join dataframes
FC_fem <- inner_join(en_fc_fem, un_fc_fem, by = "Geneid") %>% 
  select(Geneid, mean_FC_e, sd_FC_e, mean_FC_un, sd_FC_un) %>% 
  filter(is.finite(mean_FC_e), 
         is.finite(mean_FC_un))

#plot to take a look; looks pretty good
ggplot(data = FC_fem, aes(x = mean_FC_un, y = mean_FC_e)) +
  geom_point(alpha = 0.5) +
  scale_y_log10(breaks = 10^(-2:3),
                labels = scales::label_parse()(paste0("10^", -2:3))) +
  scale_x_log10(breaks = 10^(-2:3),
                labels = scales::label_parse()(paste0("10^", -2:3))) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "indianred") +
  labs(x = "Mean(CPM/sod2)\n[Unenriched]", 
       y = "Mean(CPM/sod2)\n[Enriched]") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    margins = margin(r = 15, l = 15, t = 15, b = 15),
    axis.title.x = element_text(vjust = .5, size = 12),
    axis.title.y = element_text(hjust = .5, size = 12),
    axis.text = element_text(size = 10, colour = "black", font = "arial"),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.line.y.right = element_blank(),
    axis.line.x.top = element_blank()
  )

#TOST on line for intercept of 0, and slope of 1; claim equivalence for relative expression patterns

#Linear Model
fem_lm <- lm(formula = log10(mean_FC_e) ~ log10(mean_FC_un), data = FC_fem)
#plot(fem_lm) #not the most normal thing in the world; but alright and lm's are pretty robust against escape from normality
fem_lm_res <- summary(fem_lm)

#Intercepy TOST
fem_int_tost <- tsum_TOST(m1 = fem_lm_res$coefficients["(Intercept)","Estimate"], 
                       sd1 = fem_lm_res$coefficients["(Intercept)", "Std. Error"], 
                       n1 = nrow(FC_fem), 
                       eqb = 0.1, eqbound_type = "raw", var.equal = FALSE)

#Slope TOST
fem_slope_tost <- tsum_TOST(m1 = fem_lm_res$coefficients["log10(mean_FC_un)","Estimate"] - 1, 
                             sd1 = fem_lm_res$coefficients["log10(mean_FC_un)", "Std. Error"], 
                             n1 = nrow(FC_fem), 
                             mu = 1,
                             eqb = 0.1, eqbound_type = "raw", var.equal = FALSE)

fem_int_tost   #Intercept is clearly different from 0, and not practically equivlent to 0 with 10% bounds
fem_slope_tost #Intercept is clearly different from 0, and practically equivalent to 1 within 10% bounds

#Small fold-changes appear to be smaller in magnitude in the enriched data, relative to the unenriched
#Perhaps this is because BC enrichment is giving a closer estimate of true transcript numbers; but why would I think that?

#TOSTs for gene by gene...
n_val <- 3
e_bounds <- 3

TOST_fem <- FC_fem %>% 
  filter(mean_FC_e > 0, 
         mean_FC_un > 0, 
         Geneid != "sod2") %>% 
  rowwise() %>% 
  mutate(TOST_up_p = as.vector(tsum_TOST(m1 = mean_FC_e, m2 = mean_FC_un, 
                                         sd1 = sd_FC_e, sd2 = sd_FC_un, 
                                         n1 = n_val, n2 = n_val, eqb = e_bounds, 
                                         eqbound_type = "raw", var.equal = FALSE,)$TOST["TOST Upper", "p.value"]), 
         TOST_low_p = as.vector(tsum_TOST(m1 = mean_FC_e, m2 = mean_FC_un, 
                                          sd1 = sd_FC_e, sd2 = sd_FC_un, 
                                          n1 = n_val, n2 = n_val, eqb = e_bounds,
                                          eqbound_type = "raw", var.equal = FALSE)$TOST["TOST Lower", "p.value"]),
         NHST_p = as.vector(tsum_TOST(m1 = mean_FC_e, m2 = mean_FC_un, 
                                      sd1 = sd_FC_e, sd2 = sd_FC_un, 
                                      n1 = n_val, n2 = n_val, eqb = e_bounds, 
                                      eqbound_type = "raw", var.equal = FALSE)$TOST["t-test", "p.value"]),
         TOST_eff_size = as.vector(tsum_TOST(m1 = mean_FC_e, m2 = mean_FC_un, 
                                             sd1 = sd_FC_e, sd2 = sd_FC_un, 
                                             n1 = n_val, n2 = n_val, eqb = e_bounds,
                                             eqbound_type = "raw", var.equal = FALSE)$effsize["Raw", "estimate"]), 
         TOST_ci_low = as.vector(tsum_TOST(m1 = mean_FC_e, m2 = mean_FC_un, 
                                           sd1 = sd_FC_e, sd2 = sd_FC_un, 
                                           n1 = n_val, n2 = n_val, eqb = e_bounds, 
                                           eqbound_type = "raw", var.equal = FALSE)$effsize["Raw", "lower.ci"]), 
         TOST_ci_up = as.vector(tsum_TOST(m1 = mean_FC_e, m2 = mean_FC_un, 
                                          sd1 = sd_FC_e, sd2 = sd_FC_un, 
                                          n1 = n_val, n2 = n_val, eqb = e_bounds, 
                                          eqbound_type = "raw", var.equal = FALSE)$effsize["Raw", "upper.ci"]))
#adjust p-vals
TOST_fem <- TOST_fem %>% 
  mutate(TOST_up_p.adjust = p.adjust(TOST_up_p, method = "BH", n = length(TOST_fem$TOST_up_p)), 
         TOST_low_p.adjust = p.adjust(TOST_low_p, method = "BH", n = length(TOST_fem$TOST_low_p)), 
         NHST_p.adjust = p.adjust(NHST_p, method = "BH", n = length(TOST_fem$NHST_p)))

#Subset to equivalence
equivalent_fem <- TOST_fem %>% 
  filter(TOST_low_p < 0.05, 
         TOST_up_p < 0.05)

nrow(equivalent_fem)

#Subset to equivalence with adjusted p-vals
equivalent_fem_adjust <- TOST_fem %>% 
  filter(TOST_low_p.adjust < 0.05, 
         TOST_up_p.adjust < 0.05)

nrow(equivalent_fem_adjust)

#Visualize


#Just the equivalent ones:

ggplot(equivalent_fem %>% filter(NHST_p > 0.05), aes(x = TOST_eff_size, y = Geneid)) +
  geom_vline(aes(xintercept = -1*e_bounds), linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = e_bounds),  linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black") +
  geom_point(size = 1, 
             alpha = 0.7, 
             colour = "red4", fill = "red3", 
             shape = 22) +
  coord_cartesian(clip = "off") +
  geom_errorbar(aes(y = Geneid, xmin = TOST_ci_low, xmax = TOST_ci_up), 
                linewidth = 1.5, 
                width = 0,
                alpha = 0.5, 
                colour = "indianred") +
  labs(x = expression(Delta~"FC (transcript/sod2)"), 
       title = "Difference in Relative Expression (Female Data)") +
  scale_x_continuous(limits = c(-5, 5)) +
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

#Everything

ggplot(TOST_fem, aes(x = TOST_eff_size, y = Geneid)) +
  geom_vline(aes(xintercept = -1*e_bounds), linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = e_bounds),  linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black") +
  geom_point(size = 1.5, 
             alpha = 0.7, 
             colour = "red4", fill = "red3", 
             shape = 22) +
  coord_cartesian(clip = "off") +
  geom_errorbar(aes(y = Geneid, xmin = TOST_ci_low, xmax = TOST_ci_up), 
                linewidth = 1, 
                width = 0,
                alpha = 0.5, 
                colour = "indianred") +
  labs(x = expression(Delta~"FC (transcript/sod2)"), 
       title = "Difference in Relative Expression (Female Data)") +
  scale_x_continuous(limits = c(-20, 20)) +
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


#Is there meaninful inflation of the enriched FCs relative to unenriched

male_FC_data <- data.frame(mean_FC_e = mean(FC_male$mean_FC_e), 
                           sd_FC_e = sd(FC_male$mean_FC_e), 
                           mean_FC_un = mean(FC_male$mean_FC_un), 
                           sd_FC_un = sd(FC_male$mean_FC_un))

ggplot(data = male_FC_data) +
  geom_point(aes(x = 1, y = mean_FC_e), shape = 22, size = 5, fill = "black") +
  geom_point(aes(x = 2, y = mean_FC_un), shape = 22, size = 5, fill = "black") +
  geom_errorbar(aes(x = 1, 
                    ymin = mean_FC_e - sd_FC_e/sqrt(436), 
                    ymax = mean_FC_e + sd_FC_e/sqrt(436)), 
                width = 0.5) +
  geom_errorbar(aes(x = 2, 
                    ymin = mean_FC_un - sd_FC_un/sqrt(436), 
                    ymax = mean_FC_un + sd_FC_un/sqrt(436)), 
                width = 0.5) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()

#Not really any difference with respect to FC to a reference; but maybe not the 
# best way to look at it; let's try with a boxplot

ggplot(data = FC_male) +
  geom_violin(aes(x = 1, y = mean_FC_e), alpha = 0.25, colour = "forestgreen") +
  ggbeeswarm::geom_beeswarm(aes(x = 1, y = mean_FC_e), cex = 2.35, method = "hex", colour = "forestgreen", alpha = 0.6) +
  geom_violin(aes(x = 2, y = mean_FC_un), alpha = 0.25, colour = "steelblue") +
  ggbeeswarm::geom_beeswarm(aes(x = 2, y = mean_FC_un), cex = 2.35, method = "hex", colour = "steelblue", alpha = .6) +
  scale_y_log10() +
  xlab("") +
  ylab("Mean FC") +
  theme_bw()


#Distribution of FC is quite similar for both; not surprising since ploting the points against
# another produces a line with a slope of ~1 and intercept of nearly 0;


#Save some of the dataframes to use in a different script
FC_data_save <- inner_join(en_fc, un_fc, by = "Geneid") %>% 
  select(Geneid, starts_with("FC_"), mean_FC_e, sd_FC_e, mean_FC_un, sd_FC_un)

FC_male_save <- inner_join(en_fc_male, un_fc_male, by = "Geneid") %>% 
  select(Geneid, starts_with("FC_"), mean_FC_e, sd_FC_e, mean_FC_un, sd_FC_un)

FC_fem_save <- inner_join(en_fc_fem, un_fc_fem, by = "Geneid") %>% 
  select(Geneid, starts_with("FC_"), mean_FC_e, sd_FC_e, mean_FC_un, sd_FC_un)


write.csv(FC_data_save, file = "data/clean/RE_data_full.csv", col.names = TRUE)
write.csv(FC_male_save, file = "data/clean/RE_data_male.csv", col.names = TRUE)
write.csv(FC_fem_save, file = "data/clean/RE_data_female.csv", col.names = TRUE)
