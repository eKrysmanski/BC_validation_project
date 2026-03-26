#Mean Centering Fun
library(tidyverse)
library(TOSTER)

#Read in the data

data_full <- read.csv(file = "data/clean/probe_targets_full.csv", header = TRUE)

unenriched_mean <- sum(data_full$CPM_mean_UN)/nrow(data_full)

enriched_mean <- sum(data_full$CPM_mean_E)/nrow(data_full)

centered <- data_full %>%
  mutate(cent_mean_E = CPM_mean_E - enriched_mean, 
         cent_mean_UN = CPM_mean_UN - unenriched_mean) %>% 
  select(Geneid, cent_mean_E, cent_mean_UN, CPM_mean_E, CPM_mean_UN)

head(centered)

##Plot

ggplot(data = centered, aes(y = cent_mean_E, x = cent_mean_UN)) +
  geom_point() + 
  scale_y_continuous(limits = c(0, 10000)) +
  scale_x_continuous(limits = c(0, 10000)) +
  geom_abline(intercept = 0, slope = 1, color = "indianred", linetype = "dashed") 
#  geom_smooth(method = "lm",formula = y ~ x, se = TRUE, colour ="steelblue") +
#  theme_bw()

#Yea that doesn't really work; the issue somewhat arises from sequencing depth so
# trying to get CPMs to be equal is an issue

#By what factor is the enriched scaled to the unenriched?

data_full <- centered %>%
  mutate(FC = CPM_mean_E/CPM_mean_UN)

mean(data_full$FC)  

ggplot(data = data_full, aes(x = Geneid, y = FC)) +
  geom_point()

#Seems to be all over the place from 2-40 fold, mean FC is ~18.4...

ggplot(data = data_full, aes(x = Geneid, y = log2(FC))) +
  geom_point()

#That's more consistent

ggplot(data = data_full, aes(x = log2(FC))) +
  geom_histogram()

mean(log2(data_full$FC))
median(log2(data_full$FC))

#Not great.... Don't think mean-centering the CPM will work, and I dont think it would
# make sense to try with raw counts since count numbers need to be normalized to the 
# sequencing depth for the run to avoid distortion. The chapter reported FC enriched vs un
# to be an average of ~13-14, let's check to make sure it's correct

data_male <- read.csv("data/clean/probe_targets_male.csv")

male_data <- data_male %>% 
  mutate(FC = CPM_mean_E13/CPM_mean_UN13)

mean(male_data$FC)
sd(male_data$FC)

ggplot(data = male_data, aes(x = FC)) +
  geom_histogram()

#Is there any pattern; see if FC correlates with counts

ggplot(data = male_data, aes(x = CPM_mean_E13, y = FC)) +
  geom_point() +
  scale_x_log10()

cor(male_data$FC, male_data$CPM_mean_E13, method = "pearson")

#Not really; suggests FC is not being inflated from already highly expressed genes
#Might be a slight trend where higher expression = higher fold change; but it's
# a weak relationship, and looks completely stochastic with the CPM_mean_UN13. 



#Try Median of Ratios normalization

data <- read.table(file = "data/STAR_all_liver_gene_counts.txt", sep = "\t", header = TRUE) %>% 
  select(Geneid, ends_with(".bam"))

rownames(data) <- data$Geneid
data <- data %>% 
  select(!Geneid)

head(data)

#Got the dataframe into the correct format

#Take the log of all values

log_data <- log(data)
head(log_data)

#Take the average of each row

log_data <- log_data %>% 
  rownames_to_column('gene') %>% 
  mutate(pseudo_reference = rowMeans(log_data))

head(log_data)

#Filter out all genes with -inf as their average

filtered_log_data <- log_data %>% 
  filter(pseudo_reference != "-Inf")

head(filtered_log_data)

#substract the gene pseudoreferences from log counts
# since its log of hte original data, it's taking a ratio

ratio_data = sweep(filtered_log_data[,2:9], 1, filtered_log_data[,10], "-")
head(ratio_data)

#Find median of the ratios for each sample

sample_medians <- apply(ratio_data, 2, median)

#Convert medians to scaling factors

scaling_factors = exp(sample_medians)
scaling_factors

#Divide original counts by the scaling factors

norm_data <- sweep(data, 2, scaling_factors, "/")

#Plot normalized data against another
targets <- read.table("data/defensome_genes_symbol_only.tsv", sep = "\t")

norm_data <- norm_data %>% 
  rownames_to_column("Geneid") %>% 
  filter(Geneid %in% targets$V1)

norm_data_sum <- norm_data %>%
  rowwise() %>%
  mutate(
    mean_E  = mean(c_across(starts_with("E")), na.rm = TRUE),
    mean_UN = mean(c_across(starts_with("L")), na.rm = TRUE)
  ) %>%
  ungroup() %>% 
  filter(mean_E > 0, 
         mean_UN > 0) %>% 
  select(Geneid, mean_E, mean_UN)

head(norm_data)

ggplot(data = norm_data_sum, aes(x = mean_UN, y = mean_E)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  geom_smooth(method = "lm", formula = y ~ x, color = "steelblue") +
  geom_abline(intercept = 0, slope = 1, colour = "red")


lm_data <- lm(log10(mean_E) ~ log10(mean_UN), data = norm_data_sum)
plot(lm_data)

summary(lm_data)

library(TOSTER)

tsum_TOST(m1 = 0.97590, sd1 = 0.01113*sqrt(559), n1 = 559, 
          m2 = 1, sd2 = 0, n2 = 559, 
          eqb = 0.05, eqbound_type = "raw", 
          var.equal = FALSE)


head(norm_data_sum)

ratios <- norm_data_sum %>% 
  mutate(ratios = mean_E/mean_UN)

ratios

ggplot(data = ratios, aes(x = ratios)) +
  geom_histogram()

mean(ratios$ratios)



#Try raw counts, as much as I don't want to mess with these it may be the only way to 
# get to what I want

log2_count <- log2(data[,2:12] + 1)

log2_count_sum <- log2_count %>%
  rownames_to_column("Geneid") %>% 
  rowwise() %>%
  mutate(
    mean_E  = mean(c_across(starts_with("E")), na.rm = TRUE),
    mean_UN = mean(c_across(starts_with("L")), na.rm = TRUE)
  ) %>%
  ungroup() %>% 
  filter(mean_E > 0, 
         mean_UN > 0) %>% 
  select(Geneid, mean_E, mean_UN)

log2_centered <- log2_count_sum %>% 
  mutate(cent_E = mean_E - mean(mean_E), 
         cent_UN = mean_UN - mean(mean_UN)) %>% 
  filter(Geneid %in% targets$V1)

ggplot(data = log2_centered, aes(x = cent_UN, y = cent_E)) +
  geom_point() +
#  scale_y_log10() +
#  scale_x_log10() +
  geom_smooth(method = "lm", formula = y ~ x, color = "steelblue") +
  geom_abline(intercept = 0, slope = 1, colour = "red")




ratios_cent <- data_cent_sum %>% 
  mutate(ratios = mean_E/mean_UN)

ggplot(data = ratios_cent, aes(y = ratios, x = Geneid)) +
  geom_point() +
  geom_hline(yintercept = 1)


#Moving on to another method of comparing individual genes;
#So I wanted to be able to somehow put enriched and unenriched on the same scale so that
# I could then compare individual genes, but so fair attempts of htis have failed. 
# This may be because enrichement is not uniform for every gene which seems true based on 
# counts like I initially predicted. 

#Going back to what I want to show for this data, I want to show that the expression patterns
# are not being distorted in the enriched vs. unenriched data. One way to do this would be to
# look at relative expression, using a reference from each data set. If the enriched data has  
# similar relative expression, then this can be used as evidence that expression patterns are 
# not being distorted. This way, I don't need to try and center the data to get them on the same 
# scale. 

#Choosing a reference gene:
# I want to choose a reference gene, that varies the least in the liver tissue for 
# both enriched and unenriched male and female tissues, so that even if there are
# sex differences, it won't distort the relative expression calculation. To choose
# what this gene will be, simply look for a sample with the lowest CV across both 
# male and female tissue for each dataset; this would suggest little change in count
# number despite differences in sex. 

data <- read.csv(file = "data/clean/probe_targets_full.csv", header = TRUE)

data_sort_en <- data %>%
  select(Geneid, ends_with("_cpm", ignore.case = FALSE), CPM_mean_E, CPM_sd_E) %>% 
  mutate(CV = CPM_sd_E/CPM_mean_E) %>% 
  arrange(CV)

data_sort_un <- data %>% 
  select(Geneid, ends_with("_CPM", ignore.case = FALSE), CPM_mean_UN, CPM_sd_UN) %>% 
  mutate(CV = CPM_sd_UN/CPM_mean_UN) %>% 
  arrange(CV)

#Check for candidates; genes with the lowest CV in the enriched, that are also in 
# the lowest 50 CV for unenriched. 
candidates <- data_sort_en %>%
  filter(Geneid %in% intersect(
    head(data_sort_en$Geneid, 50),
    head(data_sort_un$Geneid, 50)
  )) %>% 
  print()

#Both akr1a1b and slc22a15 seem to be reasonable choices; will go with akr1a1b
# since the differences between sexes appears to be less than the others. 

#For both the enriched and unenriched data; seperately...
# calculate log2(FC) for each gene in each sample, relative to akr1a1b;
# calcualate the mean log2(FC), and sd(log2(FC)). 

#Create a dataframe with the reference gene counts per sample
ref_en <- data_sort_en %>%
  filter(Geneid == "akr1a1b") %>%
  select(ends_with("_cpm", ignore.case = FALSE))

#Calculate log2(FC) for each gene relative to the appropriate reference count

en_fc <- data_sort_en %>%
  rowwise() %>%
  mutate(log2FC_E1 = log2(E1_cpm / ref_en$E1_cpm),
         log2FC_E2 = log2(E2_cpm / ref_en$E2_cpm),
         log2FC_E3 = log2(E3_cpm / ref_en$E3_cpm),
         log2FC_E4 = log2(E4_cpm / ref_en$E4_cpm),
         log2FC_E5 = log2(E5_cpm / ref_en$E5_cpm),
         log2FC_E6 = log2(E6_cpm / ref_en$E6_cpm)) %>%
  mutate(mean_FC_e = mean(c_across(starts_with("log2FC_E"))), 
         sd_FC_e = sd(c_across(starts_with("log2FC_E")))) %>% 
  ungroup()

#Do the same for the unenriched data.... 

ref_un <- data_sort_un %>% 
  filter(Geneid == "akr1a1b") %>% 
  select(ends_with("_CPM", ignore.case = FALSE))

un_fc <- data_sort_un %>% 
  rowwise() %>% 
  mutate(log2FC_UN1 = log2(L1_CPM / ref_un$L1_CPM),
         log2FC_UN2 = log2(L2_CPM / ref_un$L2_CPM),
         log2FC_UN3 = log2(L3_CPM / ref_un$L3_CPM),
         log2FC_UN4 = log2(L4_CPM / ref_un$L4_CPM),
         log2FC_UN5 = log2(L5_CPM / ref_un$L5_CPM),
         log2FC_UN6 = log2(L6_CPM / ref_un$L6_CPM)) %>%
  mutate(mean_FC_un = mean(c_across(starts_with("log2FC_UN"))), 
         sd_FC_un = sd(c_across(starts_with("log2FC_UN")))) %>% 
  ungroup() 

#Join the dataframes

FC_data <- inner_join(en_fc, un_fc, by = "Geneid") %>% 
  select(Geneid, mean_FC_e, sd_FC_e, mean_FC_un, sd_FC_un) %>% 
  filter(is.finite(mean_FC_e), 
         is.finite(mean_FC_un))

ggplot(data = FC_data, aes(x = mean_FC_un, y = mean_FC_e)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, colour = "indianred", linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
  geom_errorbar(aes(ymin = mean_FC_e - sd_FC_e/sqrt(6), 
                ymax = mean_FC_e + sd_FC_e/sqrt(6)), 
                alpha = 0.5) +
  geom_errorbarh(aes(xmin = mean_FC_un - sd_FC_un/sqrt(6), 
                xmax = mean_FC_un + sd_FC_un/sqrt(6)), 
                alpha = 0.5) +
  geom_smooth(method = "lm", formulat = y ~ x, colour = "steelblue")


#just take a peak at a quick lm; probably breaks the normality assumption but it's the 
# least important and the rest of the plots look alright.. 

lm_FC <- lm(mean_FC_e ~ mean_FC_un, data = FC_data)
plot(lm_FC)
summary(lm_FC)

#Tie me up, cover me in honey and throw me in a bear cage, the points are along a 
# line with an slope of 1 AND intercept of zero!!!!! And honestly, looks quite 
# similar to just counts vs. counts;

#There are probably some genes that we might be skeptical of using for determining 
# relative expression, assuming conventional RNA-seq is closer to the true counts
# than bait-enriched is. There is probably some information that we can infer from
# this plot; but need to think carefully about what negative vs. positive values 
# mean for a log2(FC). 

#Took a minute to add in error bars (standard error) and they're not too bad, I 
# imagine there are going to be some genes that are not considered practically equivalent
# between the datasets. 

#Okay, let's review what we did to make sure we didnt do a stupid...

#Choice of candidate gene:
#  I chose a reference by choosing a gene that had the lowest variation after scaling for
#  the count level. 

#Apparently there is a better way to decide on a candidate to use; the geNorm M-value... 
# What I did was probably good enough, but if there's a better way to handle it I 
#  should probably do that so it is defensible. 

#Going to choose my set of candidates from the CV values; I think that is a good 
# starting point for genes that have relatively little variation between the datasets


#geNorm stuff; got to get the dataframe as genes x samples; t() flips rows and columns

candidate_ref <- candidates %>% 
  select(Geneid, ends_with("_cpm"))

rownames(candidate_ref) <- candidates$Geneid
candidate_ref <- candidate_ref %>% 
  select(!Geneid) %>% 
  t()

candidate_ref

#Load in the package to calculate the M-vals
#BiocManager::install("NormqPCR")

library(NormqPCR)

results <- geNorm(candidate_ref)

results$ranking

#Seems like pretty much all of these are acceptable reference genes; the highest
# ranking were akr1a1b; and txn2


#I wonder if this also holds true if i used the unenriched CPMs

candidates_un <- data_sort_un %>%
  filter(Geneid %in% intersect(
    head(data_sort_en$Geneid, 50),
    head(data_sort_un$Geneid, 50)
  )) %>% 
  print()

candidate_ref_un <- candidates_un %>% 
  dplyr::select(Geneid, ends_with("_CPM"))

rownames(candidate_ref_un) <- candidates_un$Geneid

candidate_ref_un <- candidate_ref_un %>% 
  dplyr::select(!Geneid) %>% 
  t()

results_un <- geNorm(candidate_ref_un)

results_un$ranking

#Interesting; looks like akr1a1b is rank 5
# txn2 is rank 3 in unenriched; and rank 2 in enriched; might be a better choice
# Try a simple index:


rank_un <-data.frame(results_un$ranking) %>% 
  mutate(rank_un = 1:12) %>% 
  rename(Geneid = results_un.ranking)

rank_en <-(data.frame(results$ranking)) %>% 
  mutate(rank_en = 1:12) %>% 
  rename(Geneid = results.ranking)

rankings <- left_join(rank_en, rank_un, by = "Geneid") %>% 
  mutate(rank_total = rank_en + rank_un) %>% 
  print()

#txnl1 looks slightly better than akr1a1b; probably no clear difference in using one
# over the other. 