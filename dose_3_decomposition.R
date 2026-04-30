# ==========================================================
# dose_3_decomposition.R
# Uncertainty decomposition: compares CI width when each
# source is added in turn, and produces the full
# three-source dose distribution figure.
#
# Run after dose_2_diagnostics.R
# ==========================================================

if (!exists("estimate_dose_hybrid_mc"))
  stop("Run dose_1_setup.R first.")

SAMPLE_RES   <- 3.0
SAMPLE_TIME  <- 4.3
SAMPLE_CELLS <- 50
DECOMP_B     <- 10000

cat(strrep("=", 60), "\n")
cat("UNCERTAINTY DECOMPOSITION\n")
cat(strrep("=", 60), "\n\n")
cat(sprintf(
  "%.1f foci/cell at %.1fh | n_cells = %d | \u03c3_t = %.1fh | B = %d\n\n",
  SAMPLE_RES, SAMPLE_TIME, SAMPLE_CELLS, SIGMA_T, DECOMP_B
))

result_none <- suppressWarnings(estimate_dose_hybrid_mc(
  SAMPLE_RES, SAMPLE_TIME, B = DECOMP_B,
  n_cells = NA, phi = PHI_MEAS, sigma_t = 0
))
result_meas <- suppressWarnings(estimate_dose_hybrid_mc(
  SAMPLE_RES, SAMPLE_TIME, B = DECOMP_B,
  n_cells = SAMPLE_CELLS, phi = PHI_MEAS, sigma_t = 0
))
result_time <- suppressWarnings(estimate_dose_hybrid_mc(
  SAMPLE_RES, SAMPLE_TIME, B = DECOMP_B,
  n_cells = NA, phi = PHI_MEAS, sigma_t = SIGMA_T
))
result_full <- suppressWarnings(estimate_dose_hybrid_mc(
  SAMPLE_RES, SAMPLE_TIME, B = DECOMP_B,
  n_cells = SAMPLE_CELLS, phi = PHI_MEAS, sigma_t = SIGMA_T
))

decomp_df <- tibble::tibble(
  Scenario = c("Model only",
               "Model + Measurement",
               "Model + Time",
               "Model + Measurement + Time"),
  Median = c(result_none$summary$dose_median,
             result_meas$summary$dose_median,
             result_time$summary$dose_median,
             result_full$summary$dose_median),
  Lower  = c(result_none$summary$lower_ci,
             result_meas$summary$lower_ci,
             result_time$summary$lower_ci,
             result_full$summary$lower_ci),
  Upper  = c(result_none$summary$upper_ci,
             result_meas$summary$upper_ci,
             result_time$summary$upper_ci,
             result_full$summary$upper_ci),
  Width  = c(result_none$summary$ci_width,
             result_meas$summary$ci_width,
             result_time$summary$ci_width,
             result_full$summary$ci_width)
)

print(decomp_df, width = Inf)

p_decomp <- ggplot(
  decomp_df |>
    dplyr::mutate(Scenario = factor(Scenario, levels = rev(Scenario))),
  aes(y = Scenario)
) +
  geom_segment(aes(x = Lower, xend = Upper, yend = Scenario),
               linewidth = 2, colour = "steelblue") +
  geom_point(aes(x = Median), size = 3, colour = "darkblue") +
  labs(
    title    = "Uncertainty decomposition by source (95% HDI)",
    subtitle = sprintf(
      "%.1f foci/cell at %.1fh | n_cells = %d | \u03c3_t = %.1fh",
      SAMPLE_RES, SAMPLE_TIME, SAMPLE_CELLS, SIGMA_T),
    x = "Estimated dose (Gy)", y = NULL
  ) +
  theme_minimal(base_size = 12)

print(p_decomp)
ggsave("uncertainty_decomposition.png", plot = p_decomp,
       width = 7, height = 4, dpi = 150)

# ----------------------------------------------------------
# Full three-source dose distribution with CI convergence
# ----------------------------------------------------------
cat("\n", strrep("=", 60), "\n")
cat("FULL DISTRIBUTION AND CONVERGENCE\n")
cat(strrep("=", 60), "\n\n")

cat(sprintf("Convergence: %.1f%% (%d / %d draws)\n",
            result_full$summary$conv_rate_pct,
            result_full$summary$n_converged,
            result_full$summary$n_total))

# CI convergence plot
stab_seq <- seq(50, length(result_full$draws), by = 50)
stab_df  <- purrr::map_dfr(stab_seq, function(n) {
  h <- hdi(result_full$draws[1:n])
  tibble::tibble(n = n, lower = h[1], upper = h[2],
                 width = h[2] - h[1])
})

p_conv <- ggplot(stab_df, aes(x = n)) +
  geom_line(aes(y = width), colour = "darkblue", linewidth = 1) +
  labs(title = "CI convergence (HDI) — all uncertainty sources",
       x = "Number of successful draws",
       y = "95% HDI width (Gy)") +
  theme_minimal(base_size = 12)

print(p_conv)
ggsave("ci_convergence.png", plot = p_conv,
       width = 7, height = 4, dpi = 150)

# Full dose distribution
hdi_full <- hdi(result_full$draws)
md_full  <- result_full$summary$dose_median

p_dist <- ggplot(
  data.frame(dose = result_full$draws),
  aes(x = dose)
) +
  geom_histogram(
    aes(fill = cut(dose,
                   breaks = c(-Inf, 1, 2, Inf),
                   labels = c("Low", "Moderate", "High"))),
    bins = 60, colour = "white", alpha = 0.8, linewidth = 0.1
  ) +
  scale_fill_manual(
    values = c("Low"      = "#27ae60",
               "Moderate" = "#e6a817",
               "High"     = "#e74c3c"),
    name = NULL, drop = FALSE
  ) +
  geom_vline(xintercept = md_full,
             colour = "black", linewidth = 1.2) +
  geom_vline(xintercept = hdi_full,
             colour = "navy", linetype = "dashed",
             linewidth = 0.9) +
  geom_vline(xintercept = c(1, 2),
             colour = "darkred", linetype = "dotted",
             linewidth = 0.8) +
  annotate("label", x = md_full, y = Inf,
           label = paste("Median:", round(md_full, 2), "Gy"),
           vjust = 2, size = 3.5) +
  labs(
    title    = "Monte Carlo dose distribution (all uncertainty sources)",
    subtitle = sprintf(
      "%.1f foci/cell at %.1fh | n_cells = %d | \u03c3_t = %.1fh | n = %d draws",
      SAMPLE_RES, SAMPLE_TIME, SAMPLE_CELLS, SIGMA_T,
      result_full$summary$n_positive),
    x = "Dose (Gy)", y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        legend.key.size  = unit(0.45, "cm"))

print(p_dist)
ggsave("dose_distribution.png", plot = p_dist,
       width = 7, height = 4.5, dpi = 150)
cat("Saved dose_distribution.png and ci_convergence.png\n")
cat("\nRun dose_4_triage.R next.\n")