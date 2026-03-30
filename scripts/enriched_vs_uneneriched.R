#Simple plot to see if the data looks like it's supposed to
#Author = Evan C. Krysmanski
library(tidyverse)
library(ggplot2)
library(TOSTER)
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
  geom_point(alpha = 0.5) +
  scale_y_log10(breaks = 10^(-1:6),
                labels = scales::label_parse()(paste0("10^", -1:6))) +
  scale_x_log10(breaks = 10^(-1:5),
                labels = scales::label_parse()(paste0("10^", -1:5))) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "purple") +
  labs(x = "Mean(CPM)\n[Unenriched]", 
       y = "Mean(CPM)\n[Enriched]") +
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

#linear model for full dataset
full_lm <- lm(log10(CPM_mean_E) ~ log10(CPM_mean_UN), data = data_full)
plot(full_lm) #Looks alright; maybe a slight neg trend for scale location but data in general looks randomly scattered
summary(full_lm)    #Probably a few influential points causing this trend

#Equivalence Test against imaginary, perfectly correlated data
TOST_res_full <- tsum_TOST(m1 = 1.0108 - 1, 
                           sd1 = (0.01491*sqrt(437)), 
                           n1 = 437, 
                           eqb = 0.1,
                           eqbound_type = "raw", 
                           var.equal = FALSE)

#Passes equivalence test within 5% difference, but I found a paper that recommended;
# using 10% eqbounds for testing a slope = 1, so I'll use that.

#Let's visualize the results nicely

#Create dataframe pulling from the TOST results directly, in case anything is ever changed
TOST_full <- data.frame(
  d_slope = TOST_res_full$effsize[1,1],
  c_int_low = TOST_res_full$effsize[1,3],
  c_int_high = TOST_res_full$effsize[1,4],
  eqb_low = TOST_res_full$eqb[1,2],
  eqb_high = TOST_res_full$eqb[1,3])

#ggplot
ggplot(TOST_full, aes(x = d_slope)) +
  geom_vline(aes(xintercept = eqb_high), linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = eqb_low),  linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black") +
  geom_point(aes(y = 1), size = 5, shape = 15) +
  geom_errorbar(aes(y = 1, 
                    xmin =  c_int_low,
                    xmax =  c_int_high), 
                width = 0.35, 
                size = 0.75) +
  labs(x = expression(Delta~"slope"), 
       title = "Difference in Slope to Zero (Full Data)") +
  scale_y_continuous(limits = c(0.5, 1.5)) +
  scale_x_continuous(limits = c((TOST_full$eqb_low)*2, 
                                (TOST_full$eqb_high)*2)) +
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

#Looks good; could probably make it more pretty but it's fine. Let's move on to 
# individual comparisons of male un/en and female un/en; perhaps using the full 
# data is clouding some distortion that may exist

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
  
#Plot "male" en vs. un 

ggplot(data = male_data, aes(x = CPM_mean_UN13, y = CPM_mean_E13)) +
  geom_point(alpha = 0.5) +
  scale_y_log10(breaks = 10^(-1:6),
                labels = scales::label_parse()(paste0("10^", -1:6))) +
  scale_x_log10(breaks = 10^(-1:5),
                labels = scales::label_parse()(paste0("10^", -1:5))) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "steelblue") +
  labs(x = "Mean(CPM)\n[Unenriched]", 
       y = "Mean(CPM)\n[Enriched]") +
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

#I can see the lm has a slightly greater slope than a theoretical 1:1 line;
#This might be driven by some of those high count transcripts, or perhaps some 
# of the low count transcripts. Thinking of the system the likelihood of capturing
# the high abundance is probably much higher than the low abundance; if there is similar
# amount of probe fishing for either; might explain what I'm seeing. 

#Fit the lm

male_lm <- lm(log10(CPM_mean_E13) ~ log10(CPM_mean_UN13), data = male_data)
plot(male_lm)                     #look ok, scale location line suggests a small negative trend; but data looks randomly scattered
performance::check_model(male_lm) #Basically showing the same thing; data is a bit non-normal but should be fine for lm
summary(male_lm)

#Extract relevant info from the model summary
male_lm_sum <- summary(male_lm)

slope_male     <- male_lm_sum$coefficients[2, "Estimate"]
se_male  <- male_lm_sum$coefficients[2, "Std. Error"]
n_male         <- male_lm$df.residual + length(coef(male_lm))

#One-sample Equivalence test for a slope

TOST_res_male <- tsum_TOST(m1 = slope_male - 1, 
                           sd1 = (se_male*sqrt(n_male)), 
                           n1 = n_male, 
                           eqb = 0.1,
                           eqbound_type = "raw", var.equal = FALSE)

#Let's visualize these results
TOST_male <- data.frame(
  d_slope = 1 + TOST_res_male$effsize[1,1],
  c_int_low = TOST_res_male$effsize[1,3],
  c_int_high = TOST_res_male$effsize[1,4],
  eqb_high = TOST_res_male$eqb[1,2],
  eqb_low = TOST_res_male$eqb[1,3])

#ggplot
ggplot(TOST_male, aes(x = d_slope - 1)) +
  geom_vline(aes(xintercept = eqb_high), linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = eqb_low),  linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black") +
  geom_point(aes(y = 0.75), 
             size = 5, 
             shape = 15, 
             colour = "navy", 
             fill = "darkblue") +
  geom_errorbar(aes(y = 0.75, 
                    xmin = c_int_low,
                    xmax = c_int_high), 
                width = 0.1, 
                size = 0.75, 
                colour = "navy") +
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

#So equivalence within 10% eqb; but the difference in slope from zero is clear; 
# So there is a small effect, but it's negligable based on the bounds we set. 
# Positive slope suggests high count genes are more enriched than lower count genes
# This is not super concerning; and potentially reflects the true counts, but that's
# difficult to conclude since at best these measurements are estimates of transcript 
# numbers. 

#Female Data____________________________________________________________________

female_data <- inner_join(enriched_46_sum, unenriched_46_sum, by = "Geneid")

#Plot 46 against other

ggplot(data = female_data, aes(x = CPM_mean_UN46, y = CPM_mean_E46)) +
  geom_point(alpha = 0.5) +
  scale_y_log10(breaks = 10^(-1:6),
                labels = scales::label_parse()(paste0("10^", -1:6))) +
  scale_x_log10(breaks = 10^(-1:5),
                labels = scales::label_parse()(paste0("10^", -1:5))) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "indianred") +
  labs(x = "Mean(CPM)\n[Unenriched]", 
       y = "Mean(CPM)\n[Enriched]") +
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

#Constructing linear model
fem_lm <- lm(log10(CPM_mean_E46) ~ log10(CPM_mean_UN46), data = female_data)
#performance::check_model(fem_lm)    #not sure I really like performance; maybe I'm just used to base R plots
#plot(fem_lm)

summary(fem_lm)

#Pulling out relevant numbers
fem_lm_sum <- summary(fem_lm)

slope_fem     <- fem_lm_sum$coefficients[2, "Estimate"]
se_fem  <- fem_lm_sum$coefficients[2, "Std. Error"]
n_fem         <- fem_lm$df.residual + length(coef(fem_lm))

#TOST Test for Female Tissues
TOST_res_fem <- tsum_TOST(m1 = slope_fem - 1,
                          sd1 = se_fem*sqrt(n_fem),
                          n1 = n_fem,
                          eqb = 0.1,
                          var.equal = FALSE, 
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
  geom_point(aes(y = 1), 
             size = 5, 
             shape = 15, 
             colour = "#E573A0", 
             fill = "#E573A0") +
  geom_errorbar(aes(y = 1, 
                    xmin = c_int_low,
                    xmax = c_int_high), 
                width = 0.1, 
                size = 0.75, 
                colour = "#E573A0") +
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
# The standard error is 0.03046 for the slope; and d_slope is ~0.4. 
#Based on the reccomendations of that article, and this; I think that 10% is reasonable
# albeit a little less strict that I would have gone with; maybe a happy medium of 7.5?

#Slope of the fitted lm is 0.95970, with SE of 0.03046; this is a much greater
# difference than observed plotting enriched v. unenriched female; I think 5% is 
# probably more than reasonable to say the slope difference is practically negligable. 
#I think we could, and maybe should relax it to 7.5% or 10%

#Some of these dataframes are actually quite useful; going to save some of them

write.csv(data_full, file = "data/clean/probe_targets_full.csv")
write.csv(male_data, file = "data/clean/probe_targets_male.csv")
write.csv(female_data, file = "data/clean/probe_targets_female.csv")


#Try combining all three plots;
ggplot(TOST_fem, aes(x = d_slope)) +
  geom_vline(aes(xintercept = eqb_high), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(aes(xintercept = eqb_low),  linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "indianred", alpha = 0.5) +
  geom_point(aes(y = 1), 
             size = 5, 
             shape = 23, 
             colour = "indianred", 
             fill = "indianred") +
  geom_errorbar(aes(y = 1, 
                    xmin = c_int_low,
                    xmax = c_int_high), 
                width = 0, 
                size = 1.5,
                alpha = 0.5,
                colour = "indianred") +
  geom_text(aes(x = d_slope, y = 1, label = "Female Data"), 
            size = 4, colour = "indianred", 
            vjust = -2.5) +
  geom_point(data = TOST_male,
             aes(y = 0.75, x = d_slope - 1), 
             size = 5, 
             shape = 23, 
             colour = "navy", 
             fill = "darkblue") +
  geom_errorbar(data = TOST_male, 
                aes(y = 0.75, x = d_slope, 
                    xmin = c_int_low,
                    xmax = c_int_high), 
                width = 0, 
                size = 1.5,
                alpha = 0.5,
                colour = "navy") + 
  geom_text(data = TOST_male, 
            aes(x = d_slope - 1, y = 0.75, label = "Male Data"), 
            size = 4, colour = "navy", 
            vjust = -2.5) +
  geom_point(data = TOST_full,
             aes(y = 1.25), 
             size = 5, 
             shape = 23, 
             colour = "black", 
             fill = "black") +
  geom_errorbar(data = TOST_full, 
                aes(y = 1.25, x = d_slope, 
                    xmin = c_int_low,
                    xmax = c_int_high), 
                width = 0,
                alpha = 0.5,
                size = 1.5, 
                colour = "black") +
  geom_text(data = TOST_full,
    aes(x = d_slope, y = 1.25, label = "Full Data"), 
            size = 4, colour = "black", 
            vjust = -2.5) +
  labs(x = expression(Delta ~ "slope")) +
  scale_y_continuous(limits = c(0.6, 1.4)) +
  scale_x_continuous(limits = c(TOST_fem$eqb_low*1.5, TOST_fem$eqb_high*1.5)) +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.line.x = element_line(color = "black", width = 0.75), 
    axis.ticks.x = element_line(color = "black", width = 0.75),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(), 
    axis.title.x = element_text(size = 15, vjust = -1),
    panel.grid.major.y = element_blank(), 
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(),
    margins = margin(t = 25, r = 15, l = 15, b = 15, unit = "pt"), 
  )

