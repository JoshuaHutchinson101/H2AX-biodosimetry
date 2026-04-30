# ==========================================================
# Dose estimation setup: functional forms, covariance matrix,
# HDI function, and core Monte Carlo estimation function.
#
# Run after MV_meta_analysis.R and biexp_fit.R
# ==========================================================

library(dplyr)
library(tibble)
library(ggplot2)
library(MASS)
library(Matrix)
library(purrr)

set.seed(53)

# ----------------------------------------------------------
# Constants
# rho_within: within-time alpha-beta correlation, matches
#             the -0.5 imputed in the meta-analysis
# rho_cross:  cross-time correlation between pooled beta
#             estimates at different time points. Estimated
#             empirically as 0.847 from the Pearson correlation
#             of beta_4h vs beta_24h across labs contributing
#             at both time points.
#             Rounded to 0.8.
# phi = 1:    Poisson special case for measurement noise.
#             Could be set higher to reflect the overdispersion
#             documented by Errington et al. (2022) if
#             individual-level count distributions were available.
# ----------------------------------------------------------
PHI_MEAS   <- 1
SIGMA_T    <- 0.5
RHO_WITHIN <- -0.5
RHO_CROSS  <- 0.8
Z95        <- 1.96

# ----------------------------------------------------------
# Functional forms
# ----------------------------------------------------------
alpha_model <- function(t, p0, p1, k1) {
  p0 + p1 * exp(-k1 * t)
}

beta_model <- function(t, p1, p2, k1, k2) {
  p1 * exp(-k1 * t) + p2 * exp(-k2 * t)
}

# ----------------------------------------------------------
# NLS wrappers 
# ----------------------------------------------------------
fit_alpha_fn <- function(time, y, se) {
  nls(
    formula   = y ~ p0 + p1 * exp(-k1 * time),
    start     = as.list(coef(fit_a_pt)),
    weights   = 1 / se^2,
    algorithm = "port",
    lower     = c(0, 0, 1e-6),
    control   = nls.control(maxiter = 1000, warnOnly = TRUE)
  )
}

fit_beta_fn <- function(time, y, se) {
  nls(
    formula   = y ~ p1 * exp(-k1 * time) + p2 * exp(-k2 * time),
    start     = as.list(coef(fit_b_fixed)),
    weights   = 1 / se^2,
    algorithm = "port",
    lower     = c(0, 0, 1e-6, 1e-6),
    control   = nls.control(maxiter = 1000, warnOnly = TRUE)
  )
}

# ----------------------------------------------------------
# Build 12x12 joint covariance matrix over
# [a_t1,...,a_t6, b_t1,...,b_t6]
#
# Diagonal (j==k): variances from multivariate RE pooling
# Off-diagonal same param (j!=k): rho_cross * SE_j * SE_k
# Off-diagonal alpha-beta within same time: rho_within * SE_a * SE_b
# Cross-param cross-time: 0
# ----------------------------------------------------------
build_joint_sigma <- function(pmv,
                              rho_within = RHO_WITHIN,
                              rho_cross  = RHO_CROSS) {
  n <- nrow(pmv)
  S <- matrix(0, 2 * n, 2 * n)
  
  for (j in seq_len(n)) {
    for (k in seq_len(n)) {
      S[j, k] <- if (j == k) pmv$se_a_mv[j]^2 else
        rho_cross * pmv$se_a_mv[j] * pmv$se_a_mv[k]
    }
  }
  
  for (j in seq_len(n)) {
    for (k in seq_len(n)) {
      jj <- j + n; kk <- k + n
      S[jj, kk] <- if (j == k) pmv$se_b_mv[j]^2 else
        rho_cross * pmv$se_b_mv[j] * pmv$se_b_mv[k]
    }
  }
  
  for (j in seq_len(n)) {
    jj <- j + n
    S[j, jj] <- rho_within * pmv$se_a_mv[j] * pmv$se_b_mv[j]
    S[jj, j] <- S[j, jj]
  }
  
  as.matrix(nearPD(S)$mat)
}

# ----------------------------------------------------------
# HDI (Chen-Shao): shortest interval containing prob*100%
# of draws.
# ----------------------------------------------------------
hdi <- function(draws, prob = 0.95) {
  sorted <- sort(draws)
  n      <- length(sorted)
  gap    <- max(1, floor(prob * n))
  widths <- sorted[(gap + 1):n] - sorted[1:(n - gap)]
  lo_idx <- which.min(widths)
  c(lower = sorted[lo_idx], upper = sorted[lo_idx + gap])
}

# ----------------------------------------------------------
# g-statistic: Fieller's check for whether beta(t) is
# significantly different from zero. g >= 1 means the
# upper CI is unbounded.
# ----------------------------------------------------------
g_stat <- function(beta_hat, se_beta, z = Z95) {
  (z^2 * se_beta^2) / beta_hat^2
}

se_beta_at_t <- function(t, pmv) {
  pmv$se_b_mv[which.min(abs(pmv$time_h - t))]
}

beta_hat_at_t <- function(t) {
  cb <- as.list(coef(fit_b_fixed))
  beta_model(t, cb$p1, cb$p2, cb$k1, cb$k2)
}

# ----------------------------------------------------------
# Pre-build joint covariance matrix once at startup —
# reused for every call to estimate_dose_hybrid_mc
# ----------------------------------------------------------
cat(sprintf(
  "Building 12x12 joint covariance matrix (rho_within=%.1f, rho_cross=%.1f)...\n",
  RHO_WITHIN, RHO_CROSS
))

SIGMA_JOINT <- build_joint_sigma(pooled_mv)
MU_JOINT    <- c(pooled_mv$a_mv, pooled_mv$b_mv)
N_TP        <- nrow(pooled_mv)

cat(sprintf("Done. Min eigenvalue: %.6f\n",
            min(eigen(SIGMA_JOINT, only.values = TRUE)$values)))

# ----------------------------------------------------------
# Core dose estimation function
# Three uncertainty sources combined per draw:
#   1. Joint MVN draw of (alpha, beta) at all time points
#      -> refit NLS models -> get perturbed decay curves
#   2. Poisson measurement noise on R_obs
#   3. Normal perturbation of t_obs
#
# Convergence tracking:
#   n_nls_ok   -- draws where both NLS refits succeeded
#   n_b_neg    -- draws with non-positive beta(t) (filtered)
#   n_neg_raw  -- draws with negative dose (filtered)
# NLS convergence rate is reported separately from the
# count of physically plausible draws used in the HDI.
# ----------------------------------------------------------
estimate_dose_hybrid_mc <- function(
    response,
    time_h,
    pmv     = pooled_mv,
    B       = 10000,
    n_cells = NA,
    phi     = PHI_MEAS,
    sigma_t = SIGMA_T
) {
  
  ca0      <- as.list(coef(fit_a_pt))
  cb0      <- as.list(coef(fit_b_fixed))
  a_pt_val <- alpha_model(time_h, ca0$p0, ca0$p1, ca0$k1)
  b_pt_val <- beta_model( time_h, cb0$p1, cb0$p2, cb0$k1, cb0$k2)
  dose_hat <- (response - a_pt_val) / b_pt_val
  
  se_b    <- se_beta_at_t(time_h, pmv)
  g       <- g_stat(b_pt_val, se_b)
  det_lim <- g >= 1
  
  # Source 2: measurement noise
  # Floor at 0.01 to handle response = 0
  if (!is.na(n_cells) && n_cells > 0) {
    se_obs  <- sqrt(phi * max(response, 0.01) / n_cells)
    r_draws <- pmax(rnorm(B, mean = response, sd = se_obs), 0)
  } else {
    se_obs  <- 0
    r_draws <- rep(response, B)
  }
  
  # Source 3: time uncertainty
  if (sigma_t > 0) {
    t_draws <- pmin(pmax(rnorm(B, mean = time_h, sd = sigma_t), 0.5), 24)
  } else {
    t_draws <- rep(time_h, B)
  }
  
  # Source 1: joint MVN draw over all time points
  joint_draws <- MASS::mvrnorm(B, mu = MU_JOINT, Sigma = SIGMA_JOINT)
  
  dose_draws <- rep(NA_real_, B)
  n_nls_ok   <- 0L
  n_b_neg    <- 0L
  n_neg_raw  <- 0L
  
  for (i in seq_len(B)) {
    fit_succeeded <- FALSE
    
    tryCatch({
      a_star <- joint_draws[i, 1:N_TP]
      b_star <- joint_draws[i, (N_TP + 1):(2 * N_TP)]
      
      fa_i <- fit_alpha_fn(pmv$time_h, a_star, pmv$se_a_mv)
      fb_i <- fit_beta_fn( pmv$time_h, b_star, pmv$se_b_mv)
      
      fit_succeeded <- TRUE
      
      ca_i <- as.list(coef(fa_i))
      cb_i <- as.list(coef(fb_i))
      
      a_i <- alpha_model(t_draws[i], ca_i$p0, ca_i$p1, ca_i$k1)
      b_i <- beta_model( t_draws[i], cb_i$p1, cb_i$p2,
                         cb_i$k1, cb_i$k2)
      
      # Filter draws with near-zero beta to avoid numerical
      # blow-up from division by very small denominators
      if (b_i > 1e-6) {
        d_i <- (r_draws[i] - a_i) / b_i
        if (d_i > 0) {
          dose_draws[i] <- d_i
        } else {
          n_neg_raw <- n_neg_raw + 1L
        }
      } else {
        n_b_neg <- n_b_neg + 1L
      }
    }, error = function(e) NULL)
    
    if (fit_succeeded) n_nls_ok <- n_nls_ok + 1L
  }
  
  nls_conv_rate <- n_nls_ok / B
  clean         <- dose_draws[!is.na(dose_draws)]
  n_pos         <- length(clean)
  
  if (nls_conv_rate < 0.90)
    warning(sprintf("Low NLS convergence: %.1f%%", 100 * nls_conv_rate))
  
  ci <- hdi(clean, prob = 0.95)
  
  list(
    summary = tibble(
      response          = response,
      time_h            = time_h,
      n_cells           = ifelse(is.na(n_cells), NA_real_, as.numeric(n_cells)),
      phi               = phi,
      se_obs            = se_obs,
      sigma_t_h         = sigma_t,
      dose_estimate     = dose_hat,
      dose_median       = median(clean),
      lower_ci          = ci[["lower"]],
      upper_ci          = ci[["upper"]],
      ci_width          = ci[["upper"]] - ci[["lower"]],
      g_statistic       = round(g, 3),
      det_limit_flag    = det_lim,
      n_total           = B,
      n_nls_converged   = n_nls_ok,
      nls_conv_pct      = round(100 * nls_conv_rate, 1),
      n_b_negative      = n_b_neg,
      n_dose_negative   = n_neg_raw,
      n_positive        = n_pos,
      p_low             = mean(clean < 1),
      p_mod             = mean(clean >= 1 & clean < 2),
      p_high            = mean(clean >= 2)
    ),
    draws = clean
  )
}

cat("Setup complete. Run dose_2_diagnostics.R next.\n")