# =============================================================================
# HETEROGENEITY ANALYSIS — Cochran's Q and I²
# =============================================================================

library(readxl)
library(dplyr)
library(purrr)
library(metafor)
library(tidyr)

# =============================================================================
# 1. LOAD AND CLEAN DATA
# =============================================================================

raw <- read_excel(
  "C:/Users/josh1/Downloads/Calibration curves (9).xlsx",
  sheet     = "gamma-H2AX",
  col_names = FALSE
)

colnames(raw) <- c(
  "flag", "added_by", "source", "table_fig",
  "time_h", "alpha", "beta", "se_alpha", "se_beta",
  "cov_ab", "response_dist", "dispersion", "error_var",
  "blood_type", "radiation_type", "radiation_source",
  "scoring", "lab_number", "scoring_subtype", "extra1", "extra2"
)

# Remove header row, keep only rows with numeric time and both SEs
dat_se <- raw |>
  slice(-1) |>
  filter(!is.na(time_h)) |>
  mutate(
    time_h   = as.numeric(time_h),
    alpha    = as.numeric(alpha),
    beta     = as.numeric(beta),
    se_alpha = as.numeric(se_alpha),
    se_beta  = as.numeric(se_beta)
  ) |>
  filter(
    !is.na(se_alpha),
    !is.na(se_beta),
    !is.na(alpha),
    !is.na(beta)
  )

cat("Rows after cleaning:", nrow(dat_se), "\n")
cat("Time points present:", sort(unique(dat_se$time_h)), "\n")

# =============================================================================
# 2. IDENTIFY TESTABLE TIME POINTS (k >= 2)
# =============================================================================

testable_times <- dat_se |>
  group_by(time_h) |>
  summarise(k = n(), .groups = "drop") |>
  filter(k >= 2) |>
  pull(time_h)

cat("Time points with k >= 2:", testable_times, "\n")

# =============================================================================
# 3. RUN Q-TEST AND I² AT EACH TIME POINT
# =============================================================================

het_results <- map_dfr(testable_times, function(t) {
  
  d <- dat_se |> filter(time_h == t)
  
  fit_a <- rma(
    yi     = alpha,
    sei    = se_alpha,
    method = "FE",
    data   = d
  )
  
  fit_b <- rma(
    yi     = beta,
    sei    = se_beta,
    method = "FE",
    data   = d
  )
  
  tibble(
    time_h = t,
    k      = nrow(d),
    Q_a    = round(fit_a$QE,  3),
    df     = fit_a$k - 1,
    p_a    = round(fit_a$QEp, 4),
    I2_a   = round(fit_a$I2,  1),
    Q_b    = round(fit_b$QE,  3),
    p_b    = round(fit_b$QEp, 4),
    I2_b   = round(fit_b$I2,  1)
  )
})

# =============================================================================
# 4. PRINT RESULTS
# =============================================================================


cat(sprintf("%-8s %-4s %-10s %-10s %-12s %-10s %-10s %-12s\n",
            "Time(h)", "k",
            "Q(alpha)", "p(alpha)", "I2(alpha)%",
            "Q(beta)",  "p(beta)",  "I2(beta)%"))
cat(rep("-", 78), "\n", sep = "")

for (i in seq_len(nrow(het_results))) {
  r <- het_results[i, ]
  cat(sprintf("%-8.1f %-4d %-10.3f %-10.4f %-12.1f %-10.3f %-10.4f %-12.1f\n",
              r$time_h, r$k,
              r$Q_a, r$p_a, r$I2_a,
              r$Q_b, r$p_b, r$I2_b))
}


