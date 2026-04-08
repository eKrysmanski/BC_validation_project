#Authour: Evan Krysmanski
#Looking at evidence for "non-equivalence"
#Objective is to use two one-tail t-tests with a margin to determine if difference
# are clearly greater than or less than the predefined equivalence bounds
# Practically; we don't really care if genes are not clearly equivalent within strict
# bounds, but may be concerned if they are clearly different from our predefined 
# boundries

library(tidyverse)
library(ggplot2)

#read in the data; saved from Relative_expression_TOSTs________________________

#start with the male data:
male_data <- read.csv(file = "data/clean/RE_data_male.csv", header = TRUE) %>% 
  select(Geneid, starts_with("FC_")) %>% 
  filter(Geneid != "mgst1.2")

male_data_tests <- male_data %>% 
  rowwise() %>% 
  mutate(lower_test_p = (t.test(x = c(FC_E1, FC_E2, FC_E3), 
                               y = c(FC_UN1, FC_UN2, FC_UN3), 
                               mu = -3, 
                               alternative = "less"))$p.value, 
         upper_test_p = (t.test(x = c(FC_E1, FC_E2, FC_E3), 
                                y = c(FC_UN1, FC_UN2, FC_UN3), 
                                mu = 3, 
                                alternative = "greater"))$p.value, 
         lower_p_adj = p.adjust(lower_test_p, method = "BH", n = length(male_data$Geneid)), 
         upper_p_adj = p.adjust(upper_test_p, method = "BH", n = length(male_data$Geneid)), 
         mean_diff = mean(c_across(c(FC_E1, FC_E2, FC_E3))) - mean(c_across(c(FC_UN1, FC_UN2, FC_UN3)))) %>% 
  ungroup()


#Looks like nothing is clearly not practically equivalent; tfa is the worst, but variation seems to be making 
# it unclear; try now to create a volcano plot with this data... 


#Female data now... 
fem_data <- read.csv(file = "data/clean/RE_data_female.csv", header = TRUE) %>% 
  select(Geneid, starts_with("FC_")) %>% 
  filter(Geneid != "sod2")

fem_data_tests <- fem_data %>% 
  rowwise() %>% 
  mutate(lower_test_p = (t.test(x = c_across(starts_with("FC_E")), 
                                y = c_across(starts_with("FC_UN")), 
                                mu = -3, 
                                alternative = "less", 
                                paired = TRUE))$p.value, 
         upper_test_p = (t.test(x = c_across(starts_with("FC_E")),
                                y = c_across(starts_with("FC_UN")), 
                                mu = 3, 
                                alternative = "greater", 
                                paired = TRUE))$p.value, 
         lower_p_adj = p.adjust(lower_test_p, method = "BH", n = length(fem_data$Geneid)), 
         upper_p_adj = p.adjust(upper_test_p, method = "BH", n = length(fem_data$Geneid)), 
         mean_diff = mean(c_across(c(FC_E4, FC_E5, FC_E6))) - mean(c_across(c(FC_UN4, FC_UN5, FC_UN6)))) %>% 
  ungroup()


#Male and female combined... 

full_data <- read.csv(file = "data/clean/RE_data_full.csv", header = TRUE) %>% 
  select(Geneid, starts_with("FC_")) %>% 
  filter(Geneid != "txn2")

full_data_tests <- full_data %>% 
  rowwise() %>% 
  mutate(lower_test_p = (t.test(x = c_across(starts_with("FC_E")), 
                                y = c_across(starts_with("FC_UN")), 
                                mu = -3, 
                                alternative = "less", 
                                paired = TRUE))$p.value, 
         upper_test_p = (t.test(x = c_across(starts_with("FC_E")), 
                                y = c_across(starts_with("FC_UN")), 
                                mu = 3, 
                                alternative = "greater", 
                                paired = TRUE))$p.value, 
         lower_p_adj = p.adjust(lower_test_p, method = "BH", n = length(full_data$Geneid)), 
         upper_p_adj = p.adjust(upper_test_p, method = "BH", n = length(full_data$Geneid)), 
         mean_diff = mean(c_across(starts_with("FC_E"))) - mean(c_across(starts_with("FC_UN")))) %>% 
  ungroup()

#After adjusting p-values seems like nothing falls clearly outside of equivalence bounds
# This is surprising since CIs on the TOST-plot for males at least, suggested that 
# there were several

#I suppose this means we can claim there is evidence that x/y genes RE values are 
#  clearly practically equivalent, and there is insufficient evidence to claim any 
#  of the values are clearly outside the equivalence bounds. Taken together, this 
#  suggests that the bait-capture enrichement data is not distorting relative 
#  expression of the counts, and therefore we believe bait-capture enrichement is 
#  suitable to use for further investigation of the targets in zebrafish. 




#Realizing these samples are paired... re do the TOST
#Try full data:
e_bounds <- 3   # your equivalence bounds
alpha <- 0.05   # for 90% CI
n_val <- 6      # number of paired samples

TOST_full_paired <- full_data %>% 
  rowwise() %>% 
  mutate(Tost_low_p = (t.test(x = c_across(starts_with("FC_E")), 
                              y = c_across(starts_with("FC_UN")),
                              mu = -3,
                              paired = TRUE, 
                              alternative = "greater"))$p.value, 
         Tost_up_p = (t.test(x = c_across(starts_with("FC_E")),
                             y = c_across(starts_with("FC_UN")),
                             mu = 3,
                             paired = TRUE,
                             alternative = "less"))$p.value, 
         NHST_p = (t.test(x = c_across(starts_with("FC_E")), 
                          y = c_across(starts_with("FC_UN")),
                          mu = 0,
                          paired = TRUE, 
                          alternative = "two.sided"))$p.value, 
         mean_diff = mean(c_across(starts_with("FC_E")) - c_across(starts_with("FC_UN"))),
         sd_diff = sd(c_across(starts_with("FC_E")) - c_across(starts_with("FC_UN"))),
         TOST_ci_low = mean_diff - qt(1 - alpha, df = n_val - 1) * sd_diff / sqrt(n_val),
         TOST_ci_up  = mean_diff + qt(1 - alpha, df = n_val - 1) * sd_diff / sqrt(n_val), 
         TOST_ci_low_adj = mean_diff - qt(1 - 2*(alpha/length(full_data$Geneid)), df = n_val - 1) * sd_diff / sqrt(n_val), 
         TOST_ci_up_adj = mean_diff + qt(1 - 2*(alpha/length(full_data$Geneid)), df = n_val - 1) * sd_diff / sqrt(n_val))

TOST_full_paired <- TOST_full_paired %>% 
  mutate(Tost_low_p_adjust = p.adjust(Tost_low_p, method = "BH", n = length(TOST_full_paired$Geneid)), 
         Tost_up_p_adjust = p.adjust(Tost_up_p, method = "BH", n = length(TOST_full_paired$Geneid)), 
         Equivalence_p = Tost_low_p < 0.05 && Tost_up_p < 0.05,
         Equivalence_CI = TOST_ci_up_adj < 3 && TOST_ci_low_adj > -3,  
         max_p_TOST = pmax(Tost_low_p_adjust, Tost_up_p_adjust), 
         log2_FC = log2(mean(c_across(starts_with("FC_E"))) / mean(c_across(starts_with("FC_UN")))), 
         half_CI = (abs(TOST_ci_low_adj) + abs(TOST_ci_up_adj))/2)

TOST_eq_full <- TOST_full_paired %>% 
  filter(Tost_low_p < 0.05, 
         Tost_up_p < 0.05)

TOST_eq_adjust_full <- TOST_male_paired %>% 
  filter(Tost_low_p_adjust < 0.05, 
         Tost_up_p_adjust < 0.05)


#Plot the results:

ggplot(TOST_full_paired, aes(x = mean_diff, y = Geneid)) +
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
  coord_cartesian(xlim = c(-15, 15), 
                  clip = "off") +
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

#Volcano Style Plot
full <- ggplot(TOST_full_paired,
       aes(x = sign(mean_diff) * log10(abs(mean_diff) + 1), 
           y = -log10(max_p_TOST),
           fill = Equivalence_CI)) +
  geom_errorbar(
    aes(xmin = sign(TOST_ci_low_adj) * log10(abs(TOST_ci_low_adj) + 1),
        xmax = sign(TOST_ci_up_adj) * log10(abs(TOST_ci_up_adj) + 1),
        colour = Equivalence_CI),
    width = 0,
    linewidth = 1,
    alpha = 0.5) +
  geom_point(shape = 21, 
             size = 2.5, 
             alpha = 0.8, 
             colour = "black") +
  scale_fill_manual(values = c("TRUE" = "purple3","FALSE" = "grey60"), 
                    breaks = c("TRUE", "FALSE"), 
                    labels = c("TRUE" = "Equivalent (n = 304)", "FALSE" = "Unclear (n = 132)")) +
  scale_colour_manual(values = c("TRUE" = "purple3", "FALSE" = "grey60"), 
                      breaks = c(TRUE, FALSE),
                      labels = c("TRUE" = "Equivalent (n = 304)", "FALSE" = "Unclear (n = 132)")) +
  geom_vline(xintercept = c(log10(3+1)*-1, log10(3+1)),
             linetype = "dashed", colour = "red") +
  geom_vline(xintercept = 0,
             linetype = "dashed", colour = "black") +
  labs(x = expression(Log[10](Delta~"Relative Expression")), 
       y = expression(-Log[10](italic(p)[TOST])),
       tag = "A.",
       fill = "Practical\nEquivalence",
       colour = "Practical\nEquivalence") +
  coord_cartesian(xlim = c(-2.5, 2.5)) +
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.margin = margin(20, 20, 20, 20),
    legend.position = c(0.85, 0.45),
    legend.key.size = unit(10, "pt"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.background = element_rect(
      fill = alpha("white", 0.8),
      colour = "black",
      linewidth = 0.3),
    text = element_text(size = 13, colour = "black"), 
    axis.text = element_text(colour = "black"), 
    axis.title.x = element_text(vjust = -1), 
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.tag = element_text(face = "bold", vjust = 1),
    axis.ticks = element_line(colour = "black"),
    axis.ticks.length = unit(4, "pt"))

#Try with Male Samples
e_bounds <- 3   # your equivalence bounds
alpha <- 0.05   # for 90% CI
n_val <- 3      # number of paired samples

TOST_male_paired <- male_data %>% 
  rowwise() %>% 
  mutate(Tost_low_p = (t.test(x = c_across(starts_with("FC_E")), 
                              y = c_across(starts_with("FC_UN")),
                              mu = -3,
                              paired = TRUE, 
                              alternative = "greater"))$p.value, 
         Tost_up_p = (t.test(x = c_across(starts_with("FC_E")),
                             y = c_across(starts_with("FC_UN")),
                             mu = 3,
                             paired = TRUE,
                             alternative = "less"))$p.value, 
         NHST_p = (t.test(x = c_across(starts_with("FC_E")), 
                          y = c_across(starts_with("FC_UN")),
                          mu = 0,
                          paired = TRUE, 
                          alternative = "two.sided"))$p.value, 
         mean_diff = mean(c_across(starts_with("FC_E")) - c_across(starts_with("FC_UN"))),
         sd_diff = sd(c_across(starts_with("FC_E")) - c_across(starts_with("FC_UN"))),
         TOST_ci_low = mean_diff - qt(1 - alpha, df = n_val - 1) * sd_diff / sqrt(n_val),
         TOST_ci_up  = mean_diff + qt(1 - alpha, df = n_val - 1) * sd_diff / sqrt(n_val), 
         TOST_ci_low_adj = mean_diff - qt(1 - 2*(alpha/length(male_data$Geneid)), df = n_val - 1) * sd_diff / sqrt(n_val), 
         TOST_ci_up_adj = mean_diff + qt(1 - 2*(alpha/length(male_data$Geneid)), df = n_val - 1) * sd_diff / sqrt(n_val))

TOST_male_paired <- TOST_male_paired %>% 
  mutate(Tost_low_p_adjust = p.adjust(Tost_low_p, method = "BH", n = length(TOST_male_paired$Geneid)), 
         Tost_up_p_adjust = p.adjust(Tost_up_p, method = "BH", n = length(TOST_male_paired$Geneid)), 
         Equivalence_p = Tost_low_p < 0.05 && Tost_up_p < 0.05,
         Equivalence_CI = TOST_ci_up_adj < 3 && TOST_ci_low_adj > -3,  
         max_p_TOST = pmax(Tost_low_p_adjust, Tost_up_p_adjust))

TOST_eq <- TOST_male_paired %>% 
  filter(Tost_low_p < 0.05, 
         Tost_up_p < 0.05)

TOST_eq_adjust <- TOST_male_paired %>% 
  filter(Tost_low_p_adjust < 0.05, 
         Tost_up_p_adjust < 0.05)

sum(TOST_male_paired$Equivalence_CI == TRUE)

#Plot the results:

ggplot(TOST_male_paired, aes(x = mean_diff, y = Geneid)) +
  geom_vline(aes(xintercept = -1*e_bounds), linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = e_bounds),  linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black") +
  geom_point(size = 1.5, 
             alpha = 0.5,
             shape = 22, 
             colour = "navy", fill = "blue4") +
  geom_errorbar(aes(y = Geneid, xmin = TOST_ci_low, xmax = TOST_ci_up), 
                colour = "steelblue", 
                width = 0, 
                linewidth = 1, 
                alpha = 0.5) +
  labs(x = expression(Delta~"FC (trancript/txn2)"), 
       title = "Difference in Slope to Zero (Full Data)") +
  coord_cartesian(xlim = c(-22, 5), 
                  clip = "off") +
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

#Volcano-style Plot

male <- ggplot(TOST_male_paired,
       aes(x = sign(mean_diff) * log10(abs(mean_diff) + 1), 
           y = -log10(max_p_TOST),
           fill = Equivalence_CI)) +
  geom_errorbar(
    aes(xmin = sign(TOST_ci_low_adj) * log10(abs(TOST_ci_low_adj) + 1),
        xmax = sign(TOST_ci_up_adj) * log10(abs(TOST_ci_up_adj) + 1),
        colour = Equivalence_CI),
    width = 0,
    linewidth = 1,
    alpha = 0.5) +
  geom_point(shape = 21, 
             size = 2.5, 
             alpha = 0.8, 
             colour = "black") +
  scale_fill_manual(values = c("TRUE" = "steelblue","FALSE" = "grey60"), 
                    breaks = c("TRUE", "FALSE"), 
                    labels = c("TRUE" = "Equivalent (n = 403)", "FALSE" = "Unclear (n = 32)")) +
  scale_colour_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey60"), 
                      breaks = c(TRUE, FALSE),
                      labels = c("TRUE" = "Equivalent (n = 403)", "FALSE" = "Unclear (n = 32)")) +
  geom_vline(xintercept = c(log10(3+1)*-1, log10(3+1)),
             linetype = "dashed", colour = "red") +
  geom_vline(xintercept = 0,
             linetype = "dashed", colour = "black") +
  labs(x = expression(Log[10](Delta~"Relative Expression")), 
       y = expression(-Log[10](italic(p)[TOST])),
       tag = "B.",
       fill = "Practical\nEquivalence",
       colour = "Practical\nEquivalence") +
  coord_cartesian(xlim = c(-2.5, 2.5)) +
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.margin = margin(20, 20, 20, 20),
    legend.position = c(0.85, 0.45),
    legend.key.size = unit(10, "pt"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.background = element_rect(
      fill = alpha("white", 0.8),
      colour = "black",
      linewidth = 0.3),
    text = element_text(size = 13, colour = "black"), 
    axis.text = element_text(colour = "black"), 
    axis.title.x = element_text(vjust = -1), 
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.tag = element_text(face = "bold", vjust = 1),
    axis.ticks = element_line(colour = "black"),
    axis.ticks.length = unit(4, "pt"))

#Time for female samples
#Try full data:
e_bounds <- 3   # your equivalence bounds
alpha <- 0.05   # for 90% CI
n_val <- 3      # number of paired samples

TOST_fem_paired <- fem_data %>% 
  rowwise() %>% 
  mutate(Tost_low_p = (t.test(x = c_across(starts_with("FC_E")), 
                              y = c_across(starts_with("FC_UN")),
                              mu = -3,
                              paired = TRUE, 
                              alternative = "greater"))$p.value, 
         Tost_up_p = (t.test(x = c_across(starts_with("FC_E")),
                             y = c_across(starts_with("FC_UN")),
                             mu = 3,
                             paired = TRUE,
                             alternative = "less"))$p.value, 
         NHST_p = (t.test(x = c_across(starts_with("FC_E")), 
                          y = c_across(starts_with("FC_UN")),
                          mu = 0,
                          paired = TRUE, 
                          alternative = "two.sided"))$p.value, 
         mean_diff = mean(c_across(starts_with("FC_E")) - c_across(starts_with("FC_UN"))),
         sd_diff = sd(c_across(starts_with("FC_E")) - c_across(starts_with("FC_UN"))),
         TOST_ci_low = mean_diff - qt(1 - alpha, df = n_val - 1) * sd_diff / sqrt(n_val),
         TOST_ci_up  = mean_diff + qt(1 - alpha, df = n_val - 1) * sd_diff / sqrt(n_val), 
         TOST_ci_low_adj = mean_diff - qt(1 - 2*(alpha/length(fem_data$Geneid)), df = n_val - 1) * sd_diff / sqrt(n_val), 
         TOST_ci_up_adj = mean_diff + qt(1 - 2*(alpha/length(fem_data$Geneid)), df = n_val - 1) * sd_diff / sqrt(n_val))

TOST_fem_paired <- TOST_fem_paired %>% 
  mutate(Tost_low_p_adjust = p.adjust(Tost_low_p, method = "BH", n = length(TOST_fem_paired$Geneid)), 
         Tost_up_p_adjust = p.adjust(Tost_up_p, method = "BH", n = length(TOST_fem_paired$Geneid)), 
         Equivalence = Tost_low_p_adjust < 0.05 && Tost_up_p_adjust < 0.05, 
         Equivalence_CI = TOST_ci_up_adj < 3 && TOST_ci_low_adj > -3, 
         max_p_TOST = pmax(Tost_low_p_adjust, Tost_up_p_adjust), 
         log2_FC = log2(mean(c_across(starts_with("FC_E"))) / mean(c_across(starts_with("FC_UN")))), 
         half_CI = (abs(TOST_ci_low_adj) + abs(TOST_ci_up_adj))/2)

TOST_eq_fem <- TOST_fem_paired %>% 
  filter(Tost_low_p < 0.05, 
         Tost_up_p < 0.05)

TOST_eq_adjust_fem <- TOST_fem_paired %>% 
  filter(Tost_low_p_adjust < 0.05, 
         Tost_up_p_adjust < 0.05)


#Plot the results:

ggplot(TOST_fem_paired, aes(x = mean_diff, y = Geneid, fill = Equivalence)) +
  geom_vline(aes(xintercept = -1*e_bounds), linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = e_bounds),  linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black") +
  geom_point(size = 1.5, 
             alpha = 0.5,
             shape = 22 ) +
       #      colour = "red2", fill = "indianred") +
  geom_errorbar(aes(y = Geneid, xmin = TOST_ci_low_adj, xmax = TOST_ci_up_adj), 
                colour = Equivalence, 
                width = 0, 
                linewidth = 1, 
                alpha = 0.5) +
  labs(x = expression(Delta~"FC (transcript/sod2)"), 
       title = "Difference in Slope to Zero (Female Data)") +
  coord_cartesian(xlim = c(-50, 50), 
                  clip = "off") +
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

#Volcano Style Plot

female <- ggplot(TOST_fem_paired,
       aes(x = sign(mean_diff) * log10(abs(mean_diff) + 1), 
           y = -log10(max_p_TOST),
           fill = Equivalence_CI)) +
  geom_errorbar(
    aes(xmin = sign(TOST_ci_low_adj) * log10(abs(TOST_ci_low_adj) + 1),
        xmax = sign(TOST_ci_up_adj) * log10(abs(TOST_ci_up_adj) + 1),
        colour = Equivalence_CI),
    width = 0,
    linewidth = 1,
    alpha = 0.5) +
  geom_point(shape = 21, 
             size = 2.5, 
             alpha = 0.8, 
             colour = "black") +
  scale_fill_manual(values = c("TRUE" = "indianred","FALSE" = "grey60"), 
                    breaks = c("TRUE", "FALSE"), 
                    labels = c("TRUE" = "Equivalent (n = 337)", "FALSE" = "Unclear (n = 99)")) +
  scale_colour_manual(values = c("TRUE" = "indianred", "FALSE" = "grey60"), 
                      breaks = c(TRUE, FALSE),
                      labels = c("TRUE" = "Equivalent (n = 337)", "FALSE" = "Unclear (n = 99)")) +
  geom_vline(xintercept = c(log10(3+1)*-1, log10(3+1)),
             linetype = "dashed", colour = "red") +
  geom_vline(xintercept = 0,
             linetype = "dashed", colour = "black") +
  labs(x = expression(Log[10](Delta~"Relative Expression")), 
       y = expression(-Log[10](italic(p)[TOST])),
       fill = "Practical\nEquivalence",
       colour = "Practical\nEquivalence", 
       tag = "C.") +
  coord_cartesian(xlim = c(-2.5, 2.5)) +
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.margin = margin(20, 20, 20, 20),
    legend.position = c(0.85, 0.45),
    legend.key.size = unit(10, "pt"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.background = element_rect(
      fill = alpha("white", 0.8),
      colour = "black",
      linewidth = 0.3),
    text = element_text(size = 13, colour = "black"), 
    axis.text = element_text(colour = "black"), 
    axis.title.x = element_text(vjust = -1), 
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.tag = element_text(face = "bold", vjust = 1),
    axis.ticks = element_line(colour = "black"),
    axis.ticks.length = unit(4, "pt"))

female
    
library(patchwork)

(full | male | female)
