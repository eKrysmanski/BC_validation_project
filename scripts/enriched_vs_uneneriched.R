#Simple plot to see if the data looks like it's supposed to
library(tidyverse)

enriched <- read.csv(file = "data/clean/CPM_liver_enriched_full.csv", header = TRUE)
unenriched <- read.csv(file = "data/clean/CPM_liver_unenriched_full.csv", header = TRUE)

#Try something; filter to the top like 600 expressed enriched

data_full <-inner_join(enriched, unenriched, by = "Geneid")

head(data_full)


ggplot(data = data_full, aes(x = CPM_mean_UN, y = CPM_mean_E)) +
  geom_point() +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000, 100000)) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  theme_bw()



#I'm going to make an assumption that samples 1-3 and 4-6 are male or female

#Do calculations for the enriched data...
#for E1-E3
enriched_13 <- enriched %>%
  select(Geneid, E1_cpm, E2_cpm, E3_cpm)

enriched_13_sum <- enriched_13 %>%
  rowwise () %>%
  mutate(CPM_mean_E13 = mean(c_across(2:4))) %>%
  filter(CPM_mean_E13 > 0)

#for E4-E6
enriched_46 <- enriched %>%
  select(Geneid, E4_cpm, E5_cpm, E6_cpm)

enriched_46_sum <- enriched_46 %>%
  rowwise () %>%
  mutate(CPM_mean_E46 = mean(c_across(2:4))) %>%
  filter(CPM_mean_E46 > 0)

#Do calculations for the unenriched data... 
#subset data
unenriched_13 <- unenriched %>%
  select(Geneid, L1_CPM, L2_CPM, L3_CPM)

#calculate mean and filter
unenriched_13_sum <- unenriched_13 %>%
  rowwise () %>%
  mutate(CPM_mean_UN13 = mean(c_across(2:4))) %>%
  filter(CPM_mean_UN13 > 0)

#Subset data
unenriched_46 <- unenriched %>%
  select(Geneid, L4_CPM, L5_CPM, L6_CPM)

#calculate mean and filter
unenriched_46_sum <- unenriched_46 %>%
  rowwise () %>%
  mutate(CPM_mean_UN46 = mean(c_across(2:4)))


#Join male un/en dataframes

full_13 <- inner_join(enriched_13_sum, unenriched_13_sum, by = "Geneid")
  

#Plot "male" or "female" against other

ggplot(data = full_13, aes(x = CPM_mean_UN13, y = CPM_mean_E13)) +
  geom_point() +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000, 100000)) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  theme_bw()


#Join female un/en dataframes

full_46 <- inner_join(enriched_46_sum, unenriched_46_sum, by = "Geneid")

#Plot 46 against other

ggplot(data = full_46, aes(x = CPM_mean_UN46, y = CPM_mean_E46)) +
  geom_point() +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000, 100000)) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  theme_bw()

#So 13 is almost certainly the males, and 46 are almost certainly the females... 
