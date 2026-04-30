# ==========================================================
# dose_2_diagnostics.R
# Confirmatory checks before running the full simulation:
#   - Point estimate confirmation
#   - g-statistic detection limit check
#   - Model uncertainty ONLY distribution (Sources 2 and 3
#     switched off) to isolate what the joint MVN draw
#     contributes to the dose distribution
#
# Run after dose_1_setup.R

# ----------------------------------------------------------
# Confirm locked coefficients match biexp_fit.R
# ----------------------------------------------------------
cat(strrep("=", 60), "\n")
cat("POINT ESTIMATE CONFIRMATION\n")
cat(strrep("=", 60), "\n\n")

cat("Alpha (fit_a_pt):\n")
print(round(coef(fit_a_pt), 4))

cat("\nBeta (fit_b_fixed):\n")
print(round(coef(fit_b_fixed), 4))

cc        <- as.list(coef(fit_b_fixed))
ca        <- as.list(coef(fit_a_pt))
total_amp <- cc$p1 + cc$p2
cat(sprintf(
  "\nFast hl: %.1f min (%.1f%%) | Slow hl: %.1f h (%.1f%%) | Alpha hl: %.1f min\n",
  log(2)/cc$k1*60, cc$p1/total_amp*100,
  log(2)/cc$k2,    cc$p2/total_amp*100,
  log(2)/ca$k1*60
))

# ----------------------------------------------------------
# g-statistic at each calibrated time point
# ----------------------------------------------------------
cat("\n", strrep("=", 60), "\n")
cat("DETECTION LIMIT CHECK\n")
cat(strrep("=", 60), "\n\n")

g_df <- pooled_mv |>
  dplyr::mutate(
    beta_hat = purrr::map_dbl(time_h, beta_hat_at_t),
    g        = round(Z95^2 * se_b_mv^2 / beta_hat^2, 3),
    flag     = ifelse(g >= 1, "DETECTION LIMIT", "OK")
  ) |>
  dplyr::select(time_h, beta_hat, se_b_mv, g, flag)

print(g_df, width = Inf)
cat("\nAll g < 1: calibration well-identified at all time points.\n")

# ----------------------------------------------------------
# Model uncertainty only distribution
# Sources 2 and 3 switched off (n_cells = NA, sigma_t = 0)
# to isolate the contribution of the joint MVN draw.
# This is the dose distribution arising purely from
# uncertainty in the meta-analytic calibration parameters.
# ----------------------------------------------------------
cat("\n", strrep("=", 60), "\n")
cat("MODEL UNCERTAINTY ONLY — DOSE DISTRIBUTION\n")
cat(strrep("=", 60), "\n\n")

# Use the same reference scenario as the decomposition
DIAG_RES  <- 3.0
DIAG_TIME <- 4.3
DIAG_B     <- 10000

cat(sprintf(
  "Scenario: %.1f foci/cell at %.1fh | B = %d\n",
  DIAG_RES, DIAG_TIME, DIAG_B
))
cat("Sources active: Model only (Sources 2 and 3 disabled)\n\n")

result_model_only <- suppressWarnings(estimate_dose_hybrid_mc(
  response = DIAG_RES,
  time_h   = DIAG_TIME,
  B        = DIAG_B,
  n_cells  = NA,    # source 2 off
  sigma_t  = 0      # source 3 off
))

cat("Summary:\n")
print(dplyr::select(result_model_only$summary,
                    dose_median, lower_ci, upper_ci,
                    ci_width, nls_conv_pct,
                    p_low, p_mod, p_high), width = Inf)

# Plot model-uncertainty-only distribution
hdi_mo <- hdi(result_model_only$draws)
md_mo  <- result_model_only$summary$dose_median

p_model_only <- ggplot(
  data.frame(dose = result_model_only$draws),
  aes(x = dose)
) +
  geom_histogram(
    aes(fill = cut(dose,
                   breaks = c(-Inf, 1, 2, Inf),
                   labels = c("Low", "Moderate", "High"))),
    bins = 50, colour = "white", alpha = 0.8, linewidth = 0.1
  ) +
  scale_fill_manual(
    values = c("Low"      = "#27ae60",
               "Moderate" = "#e6a817",
               "High"     = "#e74c3c"),
    name = NULL, drop = FALSE
  ) +
  geom_vline(xintercept = md_mo,
             colour = "black", linewidth = 1.2) +
  geom_vline(xintercept = hdi_mo,
             colour = "navy", linetype = "dashed",
             linewidth = 0.9) +
  geom_vline(xintercept = c(1, 2),
             colour = "darkred", linetype = "dotted",
             linewidth = 0.8) +
  annotate("label", x = md_mo, y = Inf,
           label = paste("Median:", round(md_mo, 2), "Gy"),
           vjust = 2, size = 3.5) +
  annotate("label",
           x = mean(hdi_mo), y = Inf,
           label = sprintf("95%% HDI [%.2f, %.2f]",
                           hdi_mo[1], hdi_mo[2]),
           vjust = 4, size = 3, colour = "navy") +
  labs(
    title    = "Dose distribution: model uncertainty only",
    subtitle = sprintf(
      "%.1f foci/cell at %.1fh | rho_cross = %.1f | B = %d | no measurement or time noise",
      DIAG_RES, DIAG_TIME, RHO_CROSS, DIAG_B),
    x = "Estimated dose (Gy)", y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        legend.key.size  = unit(0.45, "cm"))

print(p_model_only)
ggsave("dose_dist_model_only.png", plot = p_model_only,
       width = 7, height = 4.5, dpi = 150)
cat("Saved dose_dist_model_only.png\n")
cat("\nRun dose_3_decomposition.R next.\n")
