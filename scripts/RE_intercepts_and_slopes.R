#Calculate Relative Expression for each dataset... 
library(tidyverse)
library(ggplot2)

################################# FULL DATA ####################################

#read in data
data <- read.csv(file = "data/clean/probe_targets_full.csv", header = TRUE)

#filter to reference for full dataset
ref_full <- data %>%
  filter(Geneid == "akr1a1b") %>%
  select(ends_with("_cpm", ignore.case = TRUE))

RE_full <- data %>%
  rowwise() %>%
  mutate(FC_E1 = (E1_cpm / ref_full$E1_cpm),
         FC_E2 = (E2_cpm / ref_full$E2_cpm),
         FC_E3 = (E3_cpm / ref_full$E3_cpm),
         FC_E4 = (E4_cpm / ref_full$E4_cpm),
         FC_E5 = (E5_cpm / ref_full$E5_cpm),
         FC_E6 = (E6_cpm / ref_full$E6_cpm), 
         FC_UN1 = (L1_CPM / ref_full$L1_CPM),
         FC_UN2 = (L2_CPM / ref_full$L2_CPM),
         FC_UN3 = (L3_CPM / ref_full$L3_CPM),
         FC_UN4 = (L4_CPM / ref_full$L4_CPM),
         FC_UN5 = (L5_CPM / ref_full$L5_CPM),
         FC_UN6 = (L6_CPM / ref_full$L6_CPM)) %>%
  mutate(mean_FC_un = mean(c_across(starts_with("FC_UN"))), 
         sd_FC_un = sd(c_across(starts_with("FC_UN"))), 
         mean_FC_e = mean(c_across(starts_with("FC_E"))), 
         sd_FC_e = sd(c_across(starts_with("FC_E")))) %>%
  filter(Geneid != "akr1a1b") %>%
  ungroup()

#Plot the results

full <- ggplot(data = RE_full, aes(x = mean_FC_un, y = mean_FC_e)) +
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

#Check the linear model
RE_full_lm <- lm(formula = log10(mean_FC_e) ~ log10(mean_FC_un), data = RE_full)
#plot(RE_full_lm)    
RE_full_lm_res <- summary(RE_full_lm)

#Intercept TOST
full_int_tost <- tsum_TOST(m1 = RE_full_lm_res$coefficients["(Intercept)","Estimate"], 
                           sd1 = RE_full_lm_res$coefficients["(Intercept)", "Std. Error"]*sqrt(nrow(RE_full)), 
                           n1 = nrow(RE_full), 
                           eqb = 0.1, eqbound_type = "raw", var.equal = FALSE)

#Slope TOST
full_slope_tost <- tsum_TOST(m1 = RE_full_lm_res$coefficients["log10(mean_FC_un)","Estimate"] - 1, 
                             sd1 = RE_full_lm_res$coefficients["log10(mean_FC_un)", "Std. Error"]*sqrt(nrow(RE_full)), 
                             n1 = nrow(RE_full),
                             eqb = 0.1, eqbound_type = "raw", var.equal = FALSE)

full_int_tost   #Intercept is clearly different from 0, and  practically equivalent to 0 with 10% bounds
full_slope_tost #Intercept is clearly different from 0, and practically equivalent to 1 within 10% bounds

############################    MALE DATA    ###################################

data_male <- read.csv(file = "data/clean/probe_targets_male.csv", header = TRUE)

#filter to reference for male dataset

ref_male <- data_male %>%
  filter(Geneid == "mgst1.2") %>%
  select(ends_with("_cpm"))

#calculate relative expression

RE_male <- data_male %>%
  rowwise() %>%
  mutate(FC_E1 = (E1_cpm / ref_male$E1_cpm),
         FC_E2 = (E2_cpm / ref_male$E2_cpm),
         FC_E3 = (E3_cpm / ref_male$E3_cpm), 
         FC_UN1 = (L1_CPM / ref_male$L1_CPM),
         FC_UN2 = (L2_CPM / ref_male$L2_CPM),
         FC_UN3 = (L3_CPM / ref_male$L3_CPM)) %>%
  mutate(mean_FC_un = mean(c_across(starts_with("FC_UN"))), 
         sd_FC_un = sd(c_across(starts_with("FC_UN"))), 
         mean_FC_e = mean(c_across(starts_with("FC_E"))), 
         sd_FC_e = sd(c_across(starts_with("FC_E")))) %>%
  filter(Geneid != "mgst1.2") %>%
  ungroup()

#Plot the results
male <- ggplot(data = RE_male, aes(x = mean_FC_un, y = mean_FC_e)) +
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

#Check the linear model
male_lm <- lm(formula = log10(mean_FC_e) ~ log10(mean_FC_un), data = RE_male)
#plot(male_lm) #not the most normal thing in the world; but alright and lm's are pretty robust against escape from normality
male_lm_res <- summary(male_lm)

#Intercept TOST
male_int_tost <- tsum_TOST(m1 = male_lm_res$coefficients["(Intercept)","Estimate"], 
                           sd1 = male_lm_res$coefficients["(Intercept)", "Std. Error"]*sqrt(nrow(RE_male)), 
                           n1 = nrow(RE_male), 
                           eqb = 0.1, eqbound_type = "raw", var.equal = FALSE)

#Slope TOST
male_slope_tost <- tsum_TOST(m1 = male_lm_res$coefficients["log10(mean_FC_un)","Estimate"] - 1, 
                             sd1 = male_lm_res$coefficients["log10(mean_FC_un)", "Std. Error"]*sqrt(nrow(RE_male)), 
                             n1 = nrow(RE_male), 
                             mu = 1,
                             eqb = 0.1, eqbound_type = "raw", var.equal = FALSE)


male_int_tost #Intercept is clearly different from zero, but also not practically equivalent to zero with 10% bounds
male_slope_tost #Slope is clearly different from 1; but practically equivalent within 10% eqbounds

########################## FEMALE DATA #########################################

data_fem <- read.csv(file = "data/clean/probe_targets_female.csv", header = TRUE)

#Grab the reference from the dataframe
ref_fem <- data_fem %>%
  filter(Geneid == "sod2") %>%
  select(ends_with("_cpm"))

#Calculate relative expression
RE_fem <- data_fem %>%
  rowwise() %>%
  mutate(FC_E4 = (E4_cpm / ref_fem$E4_cpm),
         FC_E5 = (E5_cpm / ref_fem$E5_cpm),
         FC_E6 = (E6_cpm / ref_fem$E6_cpm), 
         FC_UN4 = (L4_CPM / ref_fem$L4_CPM),
         FC_UN5 = (L5_CPM / ref_fem$L5_CPM),
         FC_UN6 = (L6_CPM / ref_fem$L6_CPM)) %>%
  mutate(mean_FC_un = mean(c_across(starts_with("FC_UN"))), 
         sd_FC_un = sd(c_across(starts_with("FC_UN"))), 
         mean_FC_e = mean(c_across(starts_with("FC_E"))), 
         sd_FC_e = sd(c_across(starts_with("FC_E")))) %>%
  filter(Geneid != "sod2") %>%
  ungroup()

#Plot the results
fem <- ggplot(data = RE_fem, aes(x = mean_FC_un, y = mean_FC_e)) +
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

#Check the linear model
fem_lm <- lm(formula = log10(mean_FC_e) ~ log10(mean_FC_un), data = RE_fem)
#plot(fem_lm) #not the most normal thing in the world; but alright and lm's are pretty robust against escape from normality
fem_lm_res <- summary(fem_lm)

#Intercepy TOST
fem_int_tost <- tsum_TOST(m1 = fem_lm_res$coefficients["(Intercept)","Estimate"], 
                          sd1 = fem_lm_res$coefficients["(Intercept)", "Std. Error"]*sqrt(nrow(RE_fem)), 
                          n1 = nrow(RE_fem), 
                          eqb = 0.1, eqbound_type = "raw", var.equal = FALSE)

#Slope TOST
fem_slope_tost <- tsum_TOST(m1 = fem_lm_res$coefficients["log10(mean_FC_un)","Estimate"] - 1, 
                            sd1 = fem_lm_res$coefficients["log10(mean_FC_un)", "Std. Error"]*sqrt(nrow(RE_fem)), 
                            n1 = nrow(RE_fem), 
                            mu = 1,
                            eqb = 0.1, eqbound_type = "raw", var.equal = FALSE)

fem_int_tost   #Intercept is clearly different from 0, and not practically equivlent to 0 with 10% bounds
fem_slope_tost #Intercept is clearly different from 0, and practically equivalent to 1 within 10% bounds

######################## Storing new RE dataframes #############################

write.csv(RE_full, file = "data/clean/RE_data_full.csv", col.names = TRUE)
write.csv(RE_male, file = "data/clean/RE_data_male.csv", col.names = TRUE)
write.csv(RE_fem, file = "data/clean/RE_data_female.csv", col.names = TRUE)

