#Author = Evan C. Krysmanski
library(tidyverse)

#Let's take a look at the male data, why is it not practically equivalent?

male_data <- read.csv(file = "data/clean/probe_targets_male.csv", header = TRUE)

head(male_data)

#Let's look at a distribution of the counts

ggplot(data = male_data, aes(x = CPM_mean_E)) +
  geom_histogram(bins = 30, color = "white", fill = "steelblue") +
  scale_x_log10()

#Looks pretty normal; definately a slight right skew

#Before I try and subset; let's just see how much slope changes if I remove that
# super high count value. First figure out what the transcript is, then filter it
# out

male_sort <- male_data %>%
  arrange(desc(CPM_mean_E))

head(male_sort)

male_sort <- male_sort %>%
  filter(Geneid != "tfa")

head(male_sort)

#Fit linear model

male_lm <- lm(log10(CPM_mean_E) ~ log10(CPM_mean_UN), data = male_sort)
plot(male_lm)

ggplot(data = male_sort, aes(y = CPM_mean_E, x = CPM_mean_UN)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "steelblue") +
  geom_abline(intercept = 0, slope = 1, colour = "indianred", linetype = "dashed")

#Looks a bit better; 

lm_sum <- summary(male_lm)

#Equivalence test

TOST_res_male <- tsum_TOST(m1 = lm_sum$coefficients[2,1], m2 = 1, 
                           sd1 = lm_sum$coefficients[2,2]*sqrt(nrow(male_sort)), sd2 = 0,
                           n1 = nrow(male_sort), n2 = nrow(male_sort), 
                           eqb = 0.05, eqbound_type = "raw", 
                           var.equal = FALSE)

#Ok, so it was just that one point throwing off the slope, which is a bit less
# exciting/interesting, but also a simple adjustment to make... 

#Why is that point so crazy high?
#Why is it not an issue with the female data?

#Plot the results 
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
