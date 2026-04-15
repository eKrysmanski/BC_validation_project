normalized_df <- read.csv("data/clean/DESeq2_normalized_for_individual_lm.csv")
normalized_df
# Making table for individual gene regressions

normal_longer <- pivot_longer(normalized_df,
                              cols = -Geneid,              
                              names_to = "sample",
                              values_to = "expression")

sample_info <- data.frame(
  sample = colnames(normalized_df)[2:13],
  library_type = c("enriched_m","enriched_m","enriched_m","enriched_f", "enriched_f","enriched_f", "unenriched_m", "unenriched_m","unenriched_m","unenriched_f","unenriched_f","unenriched_f"))

linear_ready <- normal_longer %>%
  left_join(sample_info, by = "sample")

models <- linear_ready %>%
  group_by(Geneid) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(log10(expression+1) ~ library_type, data = .)))

results <- models %>%
  mutate(slope = map_dbl(model, ~ summary(.x)$coefficients[2,1]),sderror= map_dbl(model,~summary(.x)$coefficients[2,2]))
view(results)

TOST_results_c <- lapply(1:nrow(results), function(i) {
  tsum_TOST(
    m1 = results$slope[i],
    m2 = 0,
    sd1 = results$sderror[i],
    sd2 = 0,
    n1 = 6,
    n2 = 6,
    eqb = 0.05,
    eqbound_type = "raw",
    var.equal = FALSE
  )
})

names(TOST_results_c) <- results$Geneid

TOST_comb <- map_dfr(TOST_results_c, ~ data.frame(
  estimate = .x$effsize[1,1],
  c_int_low = .x$effsize[1,3],
  c_int_high = .x$effsize[1,4],
  t_test = .x$TOST[1,4],
  p_lower = .x$TOST[2,4],
  p_upper = .x$TOST[3,4]
), .id = "Geneid")
head(TOST_comb)
head(TOST_results_c)
view(TOST_comb)
filtered_tost <- filter(TOST_comb, p_upper <0.05 & p_lower <0.05)
filtered_tost

linear_ready_m <- linear_ready %>%
  filter(library_type %in% c("enriched_m", "unenriched_m"))

models_m <- linear_ready_m %>%
  group_by(Geneid) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(log10(expression+1) ~ library_type, data = .)))

results_m <- models_m %>%
  mutate(slope = map_dbl(model, ~ summary(.x)$coefficients[2,1]),sderror= map_dbl(model,~summary(.x)$coefficients[2,2]))
view(results_m)
TOST_results_m <- lapply(1:nrow(results_m), function(i) {
  tryCatch({
    tsum_TOST(
      m1 = df_means$enriched_m_mean[i],
      m2 = df_means$unenriched_m_mean[i],
      sd1 = df_means$enriched_m_sd[i],
      sd2 = df_means$unenriched_m_sd[i],
      n1 = 3,
      n2 = 3,
      eqb = 0.05*df_means$unenriched_m_mean[i],
      eqbound_type = "raw",
      var.equal = FALSE
    )
  }, error = function(e) NA)
})
##############################################################
 TOST_results_m <- lapply(1:nrow(results_m), function(i) {
  tsum_TOST(
    m1 = results_m$slope[i],
    sd1 = results_m$sderror[i],
    n1 = 3,
    mu = 0,                # reference value
    eqb = 0.05,
    eqbound_type = "raw"
  )
})
##############################################################
names(TOST_results_m) <- results$Geneid


TOST_s_male <- map_dfr(TOST_results_m, ~ data.frame(
  estimate = .x$effsize[1,1],
  c_int_low = .x$effsize[1,3],
  c_int_high = .x$effsize[1,4],
  t_test = .x$TOST[1,4],
  p_lower = .x$TOST[2,4],
  p_upper = .x$TOST[3,4]
), .id = "Geneid")

filtered_tost_m <- filter(TOST_s_male, p_upper <0.05 & p_lower <0.05)

# FEMALE
linear_ready_f <- linear_ready %>%
  filter(library_type %in% c("enriched_f", "unenriched_f"))

models_f <- linear_ready_f %>%
  group_by(Geneid) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(log10(expression+1) ~ library_type, data = .)))

results_f <- models_f %>%
  mutate(slope = map_dbl(model, ~ summary(.x)$coefficients[2,1]),sderror= map_dbl(model,~summary(.x)$coefficients[2,2]))

TOST_results_f <- lapply(1:nrow(results_f), function(i) {
  tsum_TOST(
    m1 = results_f$slope[i],
    m2 = 0,
    sd1 = results_f$sderror[i],
    sd2 = 0,
    n1 = 3,
    n2 = 3,
    eqb = 0.05,
    eqbound_type = "raw",
    var.equal = FALSE
  )
})

names(TOST_results_f) <- results$Geneid

TOST_s_female <- map_dfr(TOST_results_f, ~ data.frame(
  estimate = .x$effsize[1,1],
  c_int_low = .x$effsize[1,3],
  c_int_high = .x$effsize[1,4],
  t_test = .x$TOST[1,4],
  p_lower = .x$TOST[2,4],
  p_upper = .x$TOST[3,4]
), .id = "Geneid")

filtered_tost_f <- filter(TOST_s_female, p_upper <0.05 & p_lower <0.05)
filtered_tost_f
