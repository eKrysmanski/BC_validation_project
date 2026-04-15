#Just a little script to get a table for the writeup to explain process for choosing references

library(tidyverse)
library(ctrlGene)

##############################    FULL DATA  ###################################
data <- read.csv(file = "data/clean/probe_targets_full.csv", header = TRUE)


#Sort enriched and unenriched by CV
data_sort_en <- data %>%
  select(Geneid, ends_with("_cpm", ignore.case = FALSE), CPM_mean_E, CPM_sd_E) %>% 
  mutate(CV = CPM_sd_E/CPM_mean_E) %>% 
  arrange(CV)

data_sort_un <- data %>% 
  select(Geneid, ends_with("_CPM", ignore.case = FALSE), CPM_mean_UN, CPM_sd_UN) %>% 
  mutate(CV = CPM_sd_UN/CPM_mean_UN) %>% 
  arrange(CV)

#Sort the dataframes by lowest CV
candidates_full_en <- data_sort_en %>%
  mutate(rank_full_en = rank(CV, ties.method = "first"))

candidates_full_un <- data_sort_un %>%
  mutate(rank_full_un = rank(CV, ties.method = "first"))

#Pick 10-15 transcripts that appear in top 50 lowest CV for both enriched and unenriched
candidates_full <- inner_join(candidates_full_en, candidates_full_un, by = "Geneid") %>%
  filter(rank_full_en <= 50,
         rank_full_un <= 50) %>% 
  select(Geneid, ends_with("_CPM"))

#Subset the enriched and unenriched to this candidate list and format for geNorm ranking
candidates_full_en_CG <- data_sort_en %>%
  filter(Geneid %in% candidates_full$Geneid) %>% 
  select(Geneid, ends_with("_CPM")) %>% 
  column_to_rownames(var = "Geneid") %>% 
  t()

candidates_full_un_CG <- data_sort_un %>%
  filter(Geneid %in% candidates_full$Geneid) %>% 
  select(Geneid, ends_with("_CPM")) %>% 
  column_to_rownames(var = "Geneid") %>% 
  t()

#Get M-values
results_full_en <- geNorm(candidates_full_en_CG, ctVal = FALSE)

results_full_un <- geNorm(candidates_full_un_CG, ctVal = FALSE)

#Rank by m_score to help decide on best choice of reference
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
         Combined_rank_M = "total") %>% 
  print()

#Introduce the CV ranks

CV_rank_en <- candidates_full_en %>% 
  select(Geneid, CV, rank_full_en)

CV_rank_un <- candidates_full_un %>% 
  select(Geneid, CV, rank_full_un)

candidates_full_CV <- inner_join(CV_rank_en, CV_rank_un, by = "Geneid") %>%
  filter(rank_full_en <= 50,
         rank_full_un <= 50) %>% 
  select(Geneid, CV.x, CV.y, rank_full_en, rank_full_un) %>% 
  rename(Gene_ID = Geneid)

candidates_full_CV$rank_full_en <- rank(candidates_full_CV$rank_full_en, ties.method = "first")
candidates_full_CV$rank_full_un <- rank(candidates_full_CV$rank_full_un, ties.method = "first")

candidates_full_CV <- candidates_full_CV %>% 
  mutate(Combined_Rank_CV = rank_full_en + rank_full_un)

Full_list <- inner_join(rankings_full, candidates_full_CV, by = "Gene_ID") %>% 
  select(Gene_ID, Rank_en, Rank_un, Combined_rank_M, rank_full_en, rank_full_un, Combined_Rank_CV) %>% 
  rename("geNorm Rank (Enriched)" = Rank_en, 
         "geNorm Rank (Unenriched)" = Rank_un, 
         "geNorm Rank (Combined)" = Combined_rank_M, 
         "CV Rank (Enriched)" = rank_full_en, 
         "CV Rank (Unenriched)" = rank_full_un, 
         "CV Rank (Combined)" = Combined_Rank_CV)

Full_list

############################## MALE DATA  ######################################

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
  mutate(Combined_rank_M = rank_male_en + rank_male_un) %>% 
  rename("Gene_ID" = "Genes", 
         "rank_male_en_M" = "rank_male_en", 
         "rank_male_un_M" = "rank_male_un")
print()

#Introduce CV rankings

CV_rank_male_en <- candidates_male_en %>% 
  select(Geneid, CV, rank_male_en)

CV_rank_male_un <- candidates_male_un %>% 
  select(Geneid, CV, rank_male_un)

candidates_male_CV <- inner_join(CV_rank_male_en, CV_rank_male_un, by = "Geneid") %>%
  filter(rank_male_en <= 20,
         rank_male_un <= 20) %>% 
  select(Geneid, CV.x, CV.y, rank_male_en, rank_male_un) %>% 
  rename(Gene_ID = Geneid)

candidates_male_CV$rank_male_en <- rank(candidates_male_CV$rank_male_en, ties.method = "first")
candidates_male_CV$rank_male_un <- rank(candidates_male_CV$rank_male_un, ties.method = "first")

candidates_male_CV <- candidates_male_CV %>% 
  mutate(Combined_Rank_CV = rank_male_en + rank_male_un)

male_list <- inner_join(rankings_male, candidates_male_CV, by = "Gene_ID") %>% 
  select(Gene_ID, rank_male_en_M, rank_male_un_M, Combined_rank_M, rank_male_en, rank_male_un, Combined_Rank_CV) %>% 
  rename("geNorm Rank (Enriched)" = rank_male_en_M, 
         "geNorm Rank (Unenriched)" = rank_male_un_M, 
         "geNorm Rank (Combined)" = Combined_rank_M, 
         "CV Rank (Enriched)" = rank_male_en, 
         "CV Rank (Unenriched)" = rank_male_un, 
         "CV Rank (Combined)" = Combined_Rank_CV)

male_list

#############################    fem Data    ################################

data_fem <- read.csv(file = "data/clean/probe_targets_fem.csv", header = TRUE)

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
  mutate(Combined_rank_M = rank_fem_en + rank_fem_un) %>% 
  rename("Gene_ID" = "Genes", 
         "rank_fem_en_M" = "rank_fem_en", 
         "rank_fem_un_M" = "rank_fem_un") 

print(rankings_fem)

#Introduce CV rankings

CV_rank_fem_en <- candidates_fem_en %>% 
  select(Geneid, CV, rank_fem_en)

CV_rank_fem_un <- candidates_fem_un %>% 
  select(Geneid, CV, rank_fem_un)

candidates_fem_CV <- inner_join(CV_rank_fem_en, CV_rank_fem_un, by = "Geneid") %>%
  filter(rank_fem_en <= 50,
         rank_fem_un <= 50) %>% 
  select(Geneid, CV.x, CV.y, rank_fem_en, rank_fem_un) %>% 
  rename(Gene_ID = Geneid)

candidates_fem_CV$rank_fem_en <- rank(candidates_fem_CV$rank_fem_en, ties.method = "first")
candidates_fem_CV$rank_fem_un <- rank(candidates_fem_CV$rank_fem_un, ties.method = "first")

candidates_fem_CV <- candidates_fem_CV %>% 
  mutate(Combined_Rank_CV = rank_fem_en + rank_fem_un)

fem_list <- inner_join(rankings_fem, candidates_fem_CV, by = "Gene_ID") %>% 
  select(Gene_ID, rank_fem_en_M, rank_fem_un_M, Combined_rank_M, rank_fem_en, rank_fem_un, Combined_Rank_CV) %>% 
  rename("geNorm Rank (Enriched)" = rank_fem_en_M, 
         "geNorm Rank (Unenriched)" = rank_fem_un_M, 
         "geNorm Rank (Combined)" = Combined_rank_M, 
         "CV Rank (Enriched)" = rank_fem_en, 
         "CV Rank (Unenriched)" = rank_fem_un, 
         "CV Rank (Combined)" = Combined_Rank_CV)

fem_list
