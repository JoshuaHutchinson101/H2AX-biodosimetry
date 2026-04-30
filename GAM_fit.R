# ==========================================================
# Generalised Additive Model
#
# Fits separate GAMs to the multivariate RE pooled estimates
# of alpha(t) and beta(t), using inverse-variance weighting
# and penalised thin-plate regression splines (k = 6).
#
# Uncertainty propagation uses the delta method:
#   Var[Y(t)] = Var[a(t)] + D^2 Var[b(t)] + 2D Cov[a(t), b(t)]
#
# Requires: pooled_mv from MV_meta_analysis.R
# Produces:
#   gam_a, gam_b       -- fitted GAM objects
#   decay_gam          -- predicted decay curve at dose_plot
#   GAM fit figure     -- saved as gam_fit_1gy.png
# ==========================================================

library(dplyr)
library(tibble)
library(mgcv)
library(ggplot2)


# ----------------------------------------------------------
# 0. Settings
# ----------------------------------------------------------
DOSE_PLOT  <- 1.0          # dose (Gy) at which to plot decay curve
TIME_GRID  <- seq(0, 24, length.out = 300)
RHO_AB     <- -0.5         # assumed a-b correlation (matches meta-analysis)

# ----------------------------------------------------------
# 1. Fit GAMs to pooled alpha and beta
#    k = 6 basis functions (one per time point)
#    Weights = inverse-variance from multivariate RE
# ----------------------------------------------------------
gam_a <- gam(
  a_mv ~ s(time_h, k = 6),
  data    = pooled_mv,
  weights = 1 / se_a_mv^2,
  method  = "REML"
)

gam_b <- gam(
  b_mv ~ s(time_h, k = 6),
  data    = pooled_mv,
  weights = 1 / se_b_mv^2,
  method  = "REML"
)

# ----------------------------------------------------------
# 2. Diagnostics -- EDF and basis dimension check
# ----------------------------------------------------------
cat("=== GAM diagnostics: beta ===\n")
cat(sprintf("EDF (beta smooth): %.3f\n", sum(gam_b$edf)))
gam.check(gam_b)

cat("\n=== GAM diagnostics: alpha ===\n")
cat(sprintf("EDF (alpha smooth): %.3f\n", sum(gam_a$edf)))
gam.check(gam_a)

# ----------------------------------------------------------
# 3. Save GAM diagnostic plots
# ----------------------------------------------------------
png("gam_diagnostics_beta.png", width = 800, height = 700, res = 120)
par(mfrow = c(2, 2))
gam.check(gam_b)
dev.off()

png("gam_diagnostics_alpha.png", width = 800, height = 700, res = 120)
par(mfrow = c(2, 2))
gam.check(gam_a)
dev.off()

cat("Diagnostic plots saved.\n")

# ----------------------------------------------------------
# 4. Predict on fine time grid
# ----------------------------------------------------------
new_t  <- data.frame(time_h = TIME_GRID)
pred_a <- predict(gam_a, newdata = new_t, se.fit = TRUE)
pred_b <- predict(gam_b, newdata = new_t, se.fit = TRUE)

# Delta-method variance propagation
# Var[Y(t)] = Var[a] + D^2 Var[b] + 2D Cov[a,b]
# Cov[a,b] approximated as rho * SE(a) * SE(b)
cov_ab_grid <- RHO_AB * pred_a$se.fit * pred_b$se.fit
var_y       <- pred_a$se.fit^2 +
  DOSE_PLOT^2 * pred_b$se.fit^2 +
  2 * DOSE_PLOT * cov_ab_grid

decay_gam <- tibble(
  time_h   = TIME_GRID,
  response = pred_a$fit + pred_b$fit * DOSE_PLOT,
  lower    = response - 1.96 * sqrt(var_y),
  upper    = response + 1.96 * sqrt(var_y)
)

# ----------------------------------------------------------
# 5. Discrete snapshot points for overlay
#    Delta-method error bars using multivariate RE covariance
# ----------------------------------------------------------
points_at_dose <- pooled_mv |>
  mutate(
    resp_pt = a_mv + DOSE_PLOT * b_mv,
    var_pt  = se_a_mv^2 +
      DOSE_PLOT^2 * se_b_mv^2 +
      2 * DOSE_PLOT * cov_ab_mv,
    low_pt  = resp_pt - 1.96 * sqrt(var_pt),
    upp_pt  = resp_pt + 1.96 * sqrt(var_pt)
  )

# ----------------------------------------------------------
# 6. Plot and save
# ----------------------------------------------------------
p_gam <- ggplot() +
  geom_ribbon(
    data = decay_gam,
    aes(x = time_h, ymin = lower, ymax = upper),
    fill = "grey60", alpha = 0.25
  ) +
  geom_line(
    data = decay_gam,
    aes(x = time_h, y = response),
    colour = "black", linewidth = 1.1
  ) +
  geom_errorbar(
    data = points_at_dose,
    aes(x = time_h, ymin = low_pt, ymax = upp_pt),
    width = 0.4, colour = "darkblue"
  ) +
  geom_point(
    data = points_at_dose,
    aes(x = time_h, y = resp_pt),
    colour = "darkblue", size = 2.5
  ) +
  labs(
    title    = bquote(gamma*"-H2AX Decay at "*.(DOSE_PLOT)*" Gy: GAM Fit"),
    subtitle = "Line: continuous GAM with 95% CI  |  Points: multivariate RE pooled estimates",
    x        = "Time after exposure (h)",
    y        = "Predicted response (foci/cell)"
  ) +
  theme_minimal(base_size = 12)

print(p_gam)
ggsave("gam_fit_1gy.png", plot = p_gam,
       width = 7, height = 5, dpi = 150)
cat("GAM fit plot saved as gam_fit_1gy.png\n")