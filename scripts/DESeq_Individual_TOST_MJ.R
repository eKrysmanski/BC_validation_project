library(TOSTER)
library(readr)

df_means <- read.csv("data/clean/DESeq2_normalized_liver_combined.csv")
head(df_means)
view(df_means)
TOST_results_c <- lapply(1:nrow(df_means), function(i) {
  tsum_TOST(
    m1 = df_means$enriched_mean[i],
    m2 = df_means$unenriched_mean[i],
    sd1 = df_means$enriched_sd[i],
    sd2 = df_means$unenriched_sd[i],
    n1 = 6,
    n2 = 6,
    eqb = 0.05*df_means$unenriched_mean[i],
    eqbound_type = "raw",
    var.equal = FALSE
  )
})

names(TOST_results_c) <- df_means$Geneid

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

# Male

bad_male <- df_means %>%
  dplyr::filter(
    is.na(enriched_m_mean) |
      is.na(unenriched_m_mean) |
      is.na(enriched_m_sd) |
      is.na(unenriched_m_sd) |
      enriched_m_sd == 0 |
      unenriched_m_sd == 0
  )
bad_male
TOST_results_m <- lapply(1:nrow(df_means), function(i) {
  tsum_TOST(
    m1 = df_means$enriched_m_mean[i],
    m2 = df_means$unenriched_m_mean[i],
    sd1 = df_means$enriched_m_sd[i],
    sd2 = df_means$unenriched_m_sd[i],
    n1 = 3,
    n2 = 3,
    eqb = 0.05,
    eqbound_type = "raw",
    var.equal = FALSE
  )
})
TOST_results_m <- lapply(1:nrow(df_means), function(i) {
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
names(TOST_results_m) <- df_means$Geneid
TOST_results_m <- TOST_results_m[!sapply(TOST_results_m, function(x) all(is.na(x)))]

TOST_male <- map_dfr(TOST_results_m, ~ data.frame(
  estimate = .x$effsize[1,1],
  c_int_low = .x$effsize[1,3],
  c_int_high = .x$effsize[1,4],
  t_test = .x$TOST[1,4],
  p_lower = .x$TOST[2,4],
  p_upper = .x$TOST[3,4]
), .id = "Geneid")


# Female

TOST_results_f <- lapply(1:nrow(df_means), function(i) {
  tsum_TOST(
    m1 = df_means$enriched_f_mean[i],
    m2 = df_means$unenriched_f_mean[i],
    sd1 = df_means$enriched_f_sd[i],
    sd2 = df_means$unenriched_f_sd[i],
    n1 = 3,
    n2 = 3,
    eqb = 0.05*df_means$unenriched_f_mean[i],
    eqbound_type = "raw",
    var.equal = FALSE
  )
})
names(TOST_results_f) <- df_means$Geneid

TOST_female <- map_dfr(TOST_results_f, ~ data.frame(
  estimate = .x$effsize[1,1],
  c_int_low = .x$effsize[1,3],
  c_int_high = .x$effsize[1,4],
  t_test = .x$TOST[1,4],
  p_lower = .x$TOST[2,4],
  p_upper = .x$TOST[3,4]
), .id = "Geneid")


