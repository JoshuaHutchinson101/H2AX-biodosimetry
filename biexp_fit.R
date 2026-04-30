# ==========================================================
# 05_biexp_model.R
#
# Fits and compares candidate decay models for alpha(t) and
# beta(t) using AIC, BIC, and likelihood ratio tests.
#
# Alpha candidates:
#   1. Constant
#   2. Linear
#   3. Single-exponential with floor [chosen]
#
# Beta candidates:
#   1. Constant
#   2. Single-exponential, no floor
#   3. Single-exponential, with floor
#   4. Bi-exponential, with floor (non-identifiable)
#   5. Bi-exponential, no floor [chosen]
#
# Also reports:
#   - Parameter correlation matrices (identifiability check)
#   - Pointwise curve equivalence (floor vs no floor)
#   - Repair kinetics (half-lives, proportions)
#   - Sensitivity to t = 8h observation
#
# Standard errors are derived from the observed Fisher
# information matrix: Var(theta) = sigma^2 (J'WJ)^{-1},
# where J is the Jacobian of fitted values w.r.t. parameters
# and W is the diagonal matrix of inverse-variance weights.
# Requires: pooled_mv from MV_meta_analysis.R
# Produces:
#   fit_a_pt, fit_b_fixed  -- chosen models for downstream use
#   biexp_fit_1gy.png      -- bi-exponential fit figure
# ==========================================================

library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)
library(MASS)

# ----------------------------------------------------------
# 0. Functional forms
# ----------------------------------------------------------
alpha_model <- function(t, p0, p1, k1) {
  p0 + p1 * exp(-k1 * t)
}

beta_model <- function(t, p1, p2, k1, k2) {
  p1 * exp(-k1 * t) + p2 * exp(-k2 * t)
}

# ----------------------------------------------------------
# 1. Alpha model selection
# ----------------------------------------------------------
cat(strrep("=", 60), "\n")
cat("ALPHA MODEL SELECTION\n")
cat(strrep("=", 60), "\n\n")

fit_alpha_const <- nls(
  a_mv ~ p0,
  data    = pooled_mv,
  start   = list(p0 = mean(pooled_mv$a_mv)),
  weights = 1 / pooled_mv$se_a_mv^2
)

fit_alpha_lin <- nls(
  a_mv ~ p0 + p1 * time_h,
  data    = pooled_mv,
  start   = list(p0 = max(pooled_mv$a_mv), p1 = -0.05),
  weights = 1 / pooled_mv$se_a_mv^2
)

fit_alpha_exp <- nls(
  a_mv ~ p0 + p1 * exp(-k1 * time_h),
  data      = pooled_mv,
  start     = list(p0 = 0.1, p1 = 1.5, k1 = 0.5),
  weights   = 1 / pooled_mv$se_a_mv^2,
  algorithm = "port",
  lower     = c(0, 0, 1e-6),
  control   = nls.control(maxiter = 1000)
)

alpha_model_table <- tibble(
  Model = c("Constant", "Linear", "Single-exp floor [chosen]"),
  k     = c(1L, 2L, 3L),
  AIC   = round(c(AIC(fit_alpha_const),
                  AIC(fit_alpha_lin),
                  AIC(fit_alpha_exp)), 3),
  BIC   = round(c(BIC(fit_alpha_const),
                  BIC(fit_alpha_lin),
                  BIC(fit_alpha_exp)), 3)
) |>
  arrange(AIC) |>
  mutate(
    delta_AIC = round(AIC - min(AIC), 3),
    delta_BIC = round(BIC - min(BIC), 3)
  )

cat("Alpha AIC/BIC table:\n")
print(alpha_model_table, width = Inf)

# LRT: linear vs single-exp (nested, df = 2)
lrt_alpha   <- 2 * as.numeric(logLik(fit_alpha_exp) -
                                logLik(fit_alpha_lin))
lrt_alpha_p <- pchisq(lrt_alpha, df = 2, lower.tail = FALSE)

cat(sprintf(
  "\nLRT linear vs single-exp alpha: stat = %.3f, df = 2, p = %.4f\n",
  lrt_alpha, lrt_alpha_p
))

cat("\nParameter correlation matrix (single-exp alpha):\n")
print(round(cov2cor(vcov(fit_alpha_exp)), 3))

# Export chosen alpha model
fit_a_pt <- fit_alpha_exp

# ----------------------------------------------------------
# 2. Beta model selection
# ----------------------------------------------------------
cat("\n", strrep("=", 60), "\n")
cat("BETA MODEL SELECTION\n")
cat(strrep("=", 60), "\n\n")

fit_beta_const <- nls(
  b_mv ~ p0,
  data    = pooled_mv,
  start   = list(p0 = mean(pooled_mv$b_mv)),
  weights = 1 / pooled_mv$se_b_mv^2
)

fit_beta_single <- nls(
  b_mv ~ p1 * exp(-k1 * time_h),
  data      = pooled_mv,
  start     = list(p1 = 12, k1 = 0.5),
  weights   = 1 / pooled_mv$se_b_mv^2,
  algorithm = "port",
  lower     = c(0, 1e-6),
  control   = nls.control(maxiter = 1000)
)

fit_beta_single_floor <- tryCatch(
  nls(
    b_mv ~ p0 + p1 * exp(-k1 * time_h),
    data      = pooled_mv,
    start     = list(p0 = 0.1, p1 = 12, k1 = 0.5),
    weights   = 1 / pooled_mv$se_b_mv^2,
    algorithm = "port",
    lower     = c(0, 0, 1e-6),
    control   = nls.control(maxiter = 1000)
  ),
  error = function(e) {
    message("Single-exp with floor failed: ", e$message); NULL
  }
)

# Chosen model: bi-exponential, no floor
fit_b_fixed <- nls(
  b_mv ~ p1 * exp(-k1 * time_h) + p2 * exp(-k2 * time_h),
  data      = pooled_mv,
  start     = list(p1 = 11.57, p2 = 2.45,
                   k1 = 0.97,  k2 = 0.038),
  weights   = 1 / pooled_mv$se_b_mv^2,
  algorithm = "port",
  lower     = c(0, 0, 1e-6, 1e-6),
  control   = nls.control(maxiter = 1000)
)

# Five-parameter model (non-identifiable, included for comparison only)
fit_beta_bi_floor <- nls(
  b_mv ~ p0 + p1 * exp(-k1 * time_h) + p2 * exp(-k2 * time_h),
  data      = pooled_mv,
  start     = list(p0 = 0, p1 = 11.57, p2 = 2.45,
                   k1 = 0.97, k2 = 0.038),
  weights   = 1 / pooled_mv$se_b_mv^2,
  algorithm = "port",
  lower     = c(0, 0, 0, 1e-6, 1e-6),
  control   = nls.control(maxiter = 1000)
)

beta_fits <- list(
  "Constant"                  = fit_beta_const,
  "Single-exp, no floor"      = fit_beta_single,
  "Single-exp, with floor"    = fit_beta_single_floor,
  "Bi-exp, no floor [chosen]" = fit_b_fixed,
  "Bi-exp, with floor"        = fit_beta_bi_floor
)
beta_fits <- Filter(Negate(is.null), beta_fits)

beta_model_table <- map_dfr(names(beta_fits), function(nm) {
  m <- beta_fits[[nm]]
  tibble(
    Model = nm,
    k     = length(coef(m)),
    RSS   = round(sum(residuals(m)^2), 4),
    AIC   = round(AIC(m), 3),
    BIC   = round(BIC(m), 3)
  )
}) |>
  arrange(AIC) |>
  mutate(
    delta_AIC = round(AIC - min(AIC), 3),
    delta_BIC = round(BIC - min(BIC), 3)
  )

cat("Beta AIC/BIC table:\n")
print(beta_model_table, width = Inf)

# LRT: single-exp no floor vs bi-exp no floor (nested, df = 2)
lrt_beta   <- 2 * as.numeric(logLik(fit_b_fixed) -
                               logLik(fit_beta_single))
lrt_beta_p <- pchisq(lrt_beta, df = 2, lower.tail = FALSE)

cat(sprintf(
  "\nLRT single-exp vs bi-exp (no floor): stat = %.3f, df = 2, p = %.4f\n",
  lrt_beta, lrt_beta_p
))
 
# ----------------------------------------------------------
# 3. Identifiability analysis
# ----------------------------------------------------------
cat("\n", strrep("=", 60), "\n")
cat("IDENTIFIABILITY ANALYSIS\n")
cat(strrep("=", 60), "\n\n")

cat("Correlation matrix (5-param WITH floor -- non-identifiable):\n")
print(round(cov2cor(vcov(fit_beta_bi_floor)), 3))

cat("\nCorrelation matrix (4-param NO floor -- chosen):\n")
print(round(cov2cor(vcov(fit_b_fixed)), 3))

# Pointwise curve equivalence
time_fine <- seq(0.5, 24, by = 0.1)
cc5       <- as.list(coef(fit_beta_bi_floor))
cc4       <- as.list(coef(fit_b_fixed))

beta_5par <- cc5$p0 + cc5$p1 * exp(-cc5$k1 * time_fine) +
  cc5$p2 * exp(-cc5$k2 * time_fine)
beta_4par <- cc4$p1 * exp(-cc4$k1 * time_fine) +
  cc4$p2 * exp(-cc4$k2 * time_fine)

cat(sprintf(
  "\nMax pointwise difference (5-param vs 4-param): %.2e foci/cell/Gy\n",
  max(abs(beta_5par - beta_4par))
))
cat(sprintf("AIC 5-param (with floor): %.3f\n", AIC(fit_beta_bi_floor)))
cat(sprintf("AIC 4-param (no floor):   %.3f\n", AIC(fit_b_fixed)))

# ----------------------------------------------------------
# 4. LRT summary
# ----------------------------------------------------------
cat("\n", strrep("=", 60), "\n")
cat("LRT SUMMARY\n")
cat(strrep("=", 60), "\n\n")
cat(sprintf(
  "Alpha -- linear vs single-exp:   stat = %.3f, df = 2, p = %.4f\n",
  lrt_alpha, lrt_alpha_p
))
cat(sprintf(
  "Beta  -- single-exp vs bi-exp:   stat = %.3f, df = 2, p = %.4f\n",
  lrt_beta, lrt_beta_p
))

# ----------------------------------------------------------
# 5. Parameter estimates and repair kinetics
# ----------------------------------------------------------
cat("\n", strrep("=", 60), "\n")
cat("PARAMETER ESTIMATES AND KINETICS\n")
cat(strrep("=", 60), "\n\n")

cat("Alpha (single-exp with floor):\n")
print(round(summary(fit_a_pt)$coefficients, 4))

cat("\nBeta (bi-exp, no floor):\n")
print(round(summary(fit_b_fixed)$coefficients, 4))

cc        <- as.list(coef(fit_b_fixed))
ca        <- as.list(coef(fit_a_pt))
total_amp <- cc$p1 + cc$p2

cat(sprintf("\nFast half-life:   %.1f min  (%.1f%% of signal)\n",
            log(2) / cc$k1 * 60, cc$p1 / total_amp * 100))
cat(sprintf("Slow half-life:   %.1f h    (%.1f%% of signal)\n",
            log(2) / cc$k2, cc$p2 / total_amp * 100))
cat(sprintf("Alpha half-life:  %.1f min\n",
            log(2) / ca$k1 * 60))
cat(sprintf("Background floor: %.3f foci/cell\n", ca$p0))

# ----------------------------------------------------------
# 6. Sensitivity to t = 8h observation
# ----------------------------------------------------------
cat("\n", strrep("=", 60), "\n")
cat("SENSITIVITY: OMIT t = 8h\n")
cat(strrep("=", 60), "\n\n")

pooled_mv_no8 <- pooled_mv |> filter(time_h != 8)

fit_b_no8 <- tryCatch(
  nls(
    b_mv ~ p1 * exp(-k1 * time_h) + p2 * exp(-k2 * time_h),
    data      = pooled_mv_no8,
    start     = list(p1 = 11.57, p2 = 2.45,
                     k1 = 0.97,  k2 = 0.038),
    weights   = 1 / pooled_mv_no8$se_b_mv^2,
    algorithm = "port",
    lower     = c(0, 0, 1e-6, 1e-6),
    control   = nls.control(maxiter = 1000)
  ),
  error = function(e) {
    message("No-8h refit failed: ", e$message); NULL
  }
)

if (!is.null(fit_b_no8)) {
  cc_no8  <- as.list(coef(fit_b_no8))
  amp_no8 <- cc_no8$p1 + cc_no8$p2
  
  tibble(
    Quantity    = c("Fast half-life (min)",
                    "Slow half-life (h)",
                    "Fast proportion (%)"),
    `Full data` = round(c(
      log(2) / cc$k1 * 60,
      log(2) / cc$k2,
      cc$p1 / total_amp * 100
    ), 1),
    `Omit t=8h` = round(c(
      log(2) / cc_no8$k1 * 60,
      log(2) / cc_no8$k2,
      cc_no8$p1 / amp_no8 * 100
    ), 1)
  ) |>
    mutate(Change = round(`Omit t=8h` - `Full data`, 1)) |>
    print(width = Inf)
}

# ----------------------------------------------------------
# 7. Bi-exponential fit figure
# ----------------------------------------------------------
cat("\n", strrep("=", 60), "\n")
cat("BI-EXPONENTIAL FIT FIGURE\n")
cat(strrep("=", 60), "\n\n")

# Point estimate curve: alpha(t) + beta(t) * 1 Gy
response_fit <- (ca$p0 + ca$p1 * exp(-ca$k1 * time_fine)) +
  (cc$p1 * exp(-cc$k1 * time_fine) +
     cc$p2 * exp(-cc$k2 * time_fine))

biexp_df <- tibble(
  time_h   = time_fine,
  response = response_fit
)

# Snapshot points with delta-method error bars at 1 Gy
snapshots_1gy <- pooled_mv |>
  mutate(
    y_1gy  = a_mv + b_mv,
    se_1gy = sqrt(se_a_mv^2 + se_b_mv^2 + 2 * cov_ab_mv)
  )

p_biexp <- ggplot() +
  geom_line(
    data = biexp_df,
    aes(x = time_h, y = response),
    colour = "darkblue", linewidth = 1.3
  ) +
  geom_errorbar(
    data = snapshots_1gy,
    aes(x    = time_h,
        ymin = y_1gy - 1.96 * se_1gy,
        ymax = y_1gy + 1.96 * se_1gy),
    width = 0.4, colour = "black"
  ) +
  geom_point(
    data = snapshots_1gy,
    aes(x = time_h, y = y_1gy),
    colour = "darkblue", size = 3
  ) +
  labs(
    title = expression(paste(
      gamma, "-H2AX Decay at 1 Gy: Bi-Exponential Fit")),
    x     = "Time after exposure (h)",
    y     = "Predicted response (foci/cell)"
  ) +
  theme_bw(base_size = 12)

print(p_biexp)
ggsave("biexp_fit_1gy.png", plot = p_biexp,
       width = 7, height = 4.5, dpi = 300)
cat("Saved biexp_fit_1gy.png\n")