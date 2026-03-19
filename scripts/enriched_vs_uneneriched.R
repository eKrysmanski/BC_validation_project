#Simple plot to see if the data looks like it's supposed to
library(tidyverse)

#Read in the data
enriched <- read.csv(file = "data/clean/CPM_liver_enriched_full.csv", header = TRUE)
unenriched <- read.csv(file = "data/clean/CPM_liver_unenriched_full.csv", header = TRUE)

#Read in probe list (gene symbols)
targets <- read.table("data/defensome_genes_symbol_only.tsv", sep = "\t", header = FALSE)
colnames(targets) <- "gene_symbol"

#Combine data (males and females) into single dataframe
data_full <-inner_join(enriched, unenriched, by = "Geneid") %>% 
  filter(Geneid %in% targets$gene_symbol)

#Plot meanCPM enriched against unenriched to take a peek at the data
ggplot(data = data_full, aes(x = CPM_mean_UN, y = CPM_mean_E)) +
  geom_point() +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000, 100000)) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000)) +
  geom_abline(intercept = 0, slope = 1, color = "indianred", linetype = "dashed") +
  geom_smooth(method = "lm",formula = y ~ x, se = TRUE, colour ="steelblue") +
  theme_bw()

full_lm <- lm(log10(CPM_mean_E) ~ log10(CPM_mean_UN), data = data_full)
plot(full_lm) #Looks alright; maybe a slight neg trend for scale location but data in general looks randomly scattered
summary(full_lm)    #Probably a few influential points causing this trend

#Equivalence Test against imaginary, perfectly correlated data
library(TOSTER)

TOST_res_full <- tsum_TOST(m1 = 1.0108, m2 = 1, 
          sd1 = (0.01491*sqrt(437)), sd2 = 0, 
          n1 = 437, n2 = 437, 
          eqb = 0.05, eqbound_type = "raw", 
          var.equal = FALSE)

#Passes equivalence test within 5% difference; chosen arbitrarily but is probably negligable differences
#Let's visualize the results nicely

TOST_full <- data.frame(
  d_slope = TOST_res_full$effsize[1,1],
  c_int_low = TOST_res_full$effsize[1,3],
  c_int_high = TOST_res_full$effsize[1,4],
  eqb_high = TOST_res_full$eqb[1,2],
  eqb_low = TOST_res_full$eqb[1,3])

library(ggplot2)

ggplot(TOST_full, aes(x = d_slope)) +
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
  scale_x_continuous(limits = c(TOST_full$eqb_low*2, TOST_full$eqb_high*2)) +
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

#Looks good; could probably work on equivalence tests from here using the full 
# dataset (males and females pooled). Going to go ahead and look at males/females
# individually anyways to make sure it all matches up properly. 

#####################Subsetting data and filtering to targets###################
#Samples E1-E3 are males, and samples E4-E6 are females. 
#Var names are as follows:
#      enriched_13 -> enriched samples E1 to E3 filtered to probe targets
#      enriched_13_sum -> enriched samples E1 to E3, with mean CPM
#      same convention for unenriched

#Subset the data and filter to targets
enriched_13 <- enriched %>%
  select(Geneid, E1_cpm, E2_cpm, E3_cpm) %>% 
  filter(Geneid %in% targets$gene_symbol)
#calculate mean CPM for subsetted data
enriched_13_sum <- enriched_13 %>%
  rowwise () %>%
  mutate(CPM_mean_E13 = mean(c_across(2:4)))

#Subset the data and filter to targets
enriched_46 <- enriched %>%
  select(Geneid, E4_cpm, E5_cpm, E6_cpm) %>% 
  filter(Geneid %in% targets$gene_symbol)
#calculate the mean CPM for subsetted data
enriched_46_sum <- enriched_46 %>%
  rowwise () %>%
  mutate(CPM_mean_E46 = mean(c_across(2:4))) %>%
  filter(CPM_mean_E46 > 0)

#Subset data and filter to targets
unenriched_13 <- unenriched %>%
  select(Geneid, L1_CPM, L2_CPM, L3_CPM) %>% 
  filter(Geneid %in% targets$gene_symbol)
#calculate mean CPM and remove zero count data
unenriched_13_sum <- unenriched_13 %>%
  rowwise () %>%
  mutate(CPM_mean_UN13 = mean(c_across(2:4))) %>%
  filter(CPM_mean_UN13 > 0)

#Subset data and filter to targets
unenriched_46 <- unenriched %>%
  select(Geneid, L4_CPM, L5_CPM, L6_CPM) %>% 
  filter(Geneid %in% targets$gene_symbol)
#calculate mean CPM and remove zero count data
unenriched_46_sum <- unenriched_46 %>%
  rowwise () %>%
  mutate(CPM_mean_UN46 = mean(c_across(2:4)))

###################  Linear Models and Equivalence Testing #####################


#Male Data______________________________________________________________________

male_data <- inner_join(enriched_13_sum, unenriched_13_sum, by = "Geneid")
  
#Plot "male" or "female" en vs. un 

ggplot(data = male_data, aes(x = CPM_mean_UN13, y = CPM_mean_E13)) +
  geom_point() +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000, 100000)) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "steelblue") +
  theme_bw()

#Fit the lm

male_lm <- lm(log10(CPM_mean_E13) ~ log10(CPM_mean_UN13), data = male_data)
plot(male_lm) #look ok, scale location line suggests a small negative trend; but data looks randomly scattered
performance::check_model(male_lm) #Basically showing the same thing; data is a bit non-normal but should be fine for lm
summary(male_lm)

#Equivalence test against hypothetical perfectly correlated data
TOST_res_male <- tsum_TOST(m1 = 1.0523, sd1 = (0.1808*sqrt(436)), n1 = 436, 
                           m2 = 1, sd2 = 0, n2 = 436, 
                           eqb = 0.2, eqbound_type = "raw", var.equal = FALSE)

#Let's visualize these results
TOST_male <- data.frame(
  d_slope = TOST_res_male$effsize[1,1],
  c_int_low = TOST_res_male$effsize[1,3],
  c_int_high = TOST_res_male$effsize[1,4],
  eqb_high = TOST_res_male$eqb[1,2],
  eqb_low = TOST_res_male$eqb[1,3])

#ggplot
ggplot(TOST_male, aes(x = d_slope)) +
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
       title = "Difference in Slope to Zero (Male Data)") +
  scale_y_continuous(limits = c(0.5, 1.5)) +
  scale_x_continuous(limits = c(TOST_male$eqb_low*2, TOST_male$eqb_high*2)) +
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


#Female Data____________________________________________________________________

female_data <- inner_join(enriched_46_sum, unenriched_46_sum, by = "Geneid")

#Plot 46 against other

ggplot(data = female_data, aes(x = CPM_mean_UN46, y = CPM_mean_E46)) +
  geom_point(alpha = 0.6, 
             color = "magenta") +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000, 100000)) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000)) +
  geom_abline(intercept = 0, slope = 1, color = "indianred", linetype = "dashed") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "steelblue") +
  theme_bw()

fem_lm <- lm(log10(CPM_mean_E46) ~ log10(CPM_mean_UN46), data = female_data)
performance::check_model(fem_lm)    #not sure I really like performance; maybe I'm just used to base R plots
plot(fem_lm)

summary(fem_lm)
dotwhisker::dwplot(fem_lm)

#I'm pretty sure that dotwhisker plots a 95% CI
#Might be better to visualize the slope by just plotting the raw slope, with 
# raw equivalence bounds rather than the d_slope; worth considering. 
#I will definately like to plot the total, male, and female on the same plot
# at some point to show that the data is not distorted for male, female, or 
# globally (both male and fem)


#Comparing slope of the linear reg. line for female liver from enriched v unenriched
# to a hypothetical slope of 1, where there is perfect correlation. 
#Using 5% difference for eqb for now just to look at it; will try to decide on better 
# eqb to use later, but 5% is probably pretty strict... 
TOST_res_fem <- tsum_TOST(m1 = 1.01474, m2 = 1,
                          sd1 = 0.01259*sqrt(437), sd2 = 0,
                          n1 = 437, n2 = 437,
                          var.equal = FALSE,
                          mu = 0, eqb = 0.05, 
                          eqbound_type = "raw")

#From this test... 
# NHST:   non-significant; therefore it is unclear if the difference is not equal to zero
# TOST:   significant; therefore it is clear that the effect (difference in slope) is 
#                      within predefined equivalence bounds

#Let's visualize these results
TOST_fem <- data.frame(
  d_slope = TOST_res_fem$effsize[1,1],
  c_int_low = TOST_res_fem$effsize[1,3],
  c_int_high = TOST_res_fem$effsize[1,4],
  eqb_high = TOST_res_fem$eqb[1,2],
  eqb_low = TOST_res_fem$eqb[1,3])

#ggplot
ggplot(TOST_fem, aes(x = d_slope)) +
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
       title = "Difference in Slope to Zero (Female Data)") +
  scale_y_continuous(limits = c(0.5, 1.5)) +
  scale_x_continuous(limits = c(TOST_fem$eqb_low*2, TOST_fem$eqb_high*2)) +
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

#Playing around with comparisons; try and figure out an acceptable eqbound empirically... 

#Try plotting M/F for conventional RNA_seq; see how much variance there is in 
# similar tissues; despite there probably being sex differences

un_MvF <- inner_join(unenriched_13_sum, unenriched_46_sum, by = "Geneid")

ggplot(data = un_MvF, aes(x = log2(CPM_mean_UN46), y = log2(CPM_mean_UN13))) +
  geom_point(alpha = 0.6, 
             color = "black") +
  geom_abline(intercept = 0, slope = 1, color = "indianred", linetype = "dashed") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "steelblue") +
  geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed", color = "forestgreen") +
  theme_bw()

MvF_lm <- lm(log10(CPM_mean_UN13) ~ log10(CPM_mean_UN46), data = un_MvF)
plot(MvF_lm)
summary(MvF_lm)

#Plots look alright; the lm lines up almost perfectly with a 1:1 line, but not quite
# This natural varaiation is what we want to estimate for deciding on equivalence bounds
# The standard error is 0.03046 for the slope 

#Slope of the fitted lm is 0.95970, with SE of 0.03046; this is a much greater
# difference than observed plotting enriched v. unenriched female; I think 5% is 
# probably more than reasonable to say the slope difference is practically negligable. 
#I think we could, and maybe should relax it to 7.5% or 10%



#Messing around with normalization methods....

#Going to try a different way of normalizing, log2(CPM)

log2_un <- unenriched %>% 
  mutate(raw_UN1 = L1_CPM/10E6 + 1, #multiply by
         raw_UN2 = L2_CPM/10E6 + 1, 
         raw_UN3 = L3_CPM/10E6 + 1,
         log2_UN1 = log2(raw_UN1*10E6), 
         log2_UN2 = log2(raw_UN2*10E6), 
         log2_UN3 = log2(raw_UN3*10E6)) %>% 
  select(Geneid, starts_with("log2_")) %>% 
  rowwise() %>% 
  mutate(log2_mean = mean(c_across(2:4)))

log2_en <- enriched %>% 
  mutate(raw_en1 = E1_cpm/10E6 + 1,
         raw_en2 = E2_cpm/10E6 + 1, 
         raw_en3 = E3_cpm/10E6 + 1,
         log2_EN1 = log2(raw_en1*10E6), 
         log2_EN2 = log2(raw_en2*10E6), 
         log2_EN3 = log2(raw_en3*10E6)) %>% 
  select(Geneid, starts_with("log2_")) %>% 
  rowwise() %>% 
  mutate(log2_mean = mean(c_across(2:4)))

log2_male <- inner_join(log2_en, log2_un, by = "Geneid") %>% 
  filter(Geneid %in% targets$gene_symbol)

ggplot(data = log2_male, aes(x = log2_mean.x, y = log2_mean.y)) +
  geom_point(alpha = 0.6, 
             color = "black") +
  geom_abline(intercept = 0, slope = 1, color = "indianred", linetype = "dashed") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "steelblue") +
  theme_bw()


#Try another way: mean centering log2(CPM)

cent_un
