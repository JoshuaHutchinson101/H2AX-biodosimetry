# ==========================================================
# dose_4_triage.R
# Triage scenario Table 5.1 in report
# Foci: 3, 5, 10 foci/cell
# Times: 1.5, 2.5, 4, 8 hours
# Run after dose_1_setup.R
# takes a few minutes to run
# ==========================================================

cat(strrep("=", 60), "\n")
cat("TRIAGE SCENARIOS\n")
cat(strrep("=", 60), "\n\n")
cat(sprintf(
  "B = 10000 | n_cells = NA | phi = %d | sigma_t = 0h\n",
  PHI_MEAS
))
cat(sprintf(
  "rho_within = %.1f | rho_cross = %.1f\n\n",
  RHO_WITHIN, RHO_CROSS
))

scenarios <- tibble::tribble(
  ~response, ~time_h,
  3,   4.0,
  3,   8.0,
  5,   2.5,
  5,   4.0,
  5,   8.0,
  10,  1.5,
  10,  2.5,
)

triage_results <- purrr::map_dfr(
  seq_len(nrow(scenarios)),
  function(i) {
    suppressWarnings(estimate_dose_hybrid_mc(
      response = scenarios$response[i],
      time_h   = scenarios$time_h[i],
      B        = 10000,
      n_cells  = NA,
      phi      = PHI_MEAS,
      sigma_t  = 0
    ))$summary
  }
)

cat("Triage results:\n")
triage_results |>
  dplyr::select(response, time_h,
                dose_median, lower_ci, upper_ci,
                p_low, p_mod, p_high,
                nls_conv_pct) |>
  dplyr::mutate(
    across(c(dose_median, lower_ci, upper_ci),
           \(x) round(x, 3)),
    across(c(p_low, p_mod, p_high),
           \(x) round(x, 3))
  ) |>
  print(width = Inf)

cat("\nConvergence:\n")
triage_results |>
  dplyr::select(response, time_h, nls_conv_pct,
                n_nls_converged, n_positive, n_dose_negative) |>
  print(width = Inf)

# ----------------------------------------------------------
# Triage results table
# ----------------------------------------------------------
cat("\n", strrep("=", 60), "\n")
cat("TRIAGE RESULTS TABLE\n")
cat(strrep("=", 60), "\n\n")

triage_table <- triage_results |>
  dplyr::mutate(
    across(c(dose_median, lower_ci, upper_ci),
           \(x) round(x, 2)),
    across(c(p_low, p_mod, p_high),
           \(x) round(x, 3)),
    ci_95 = paste0("(", lower_ci, ", ", upper_ci, ")")
  ) |>
  dplyr::select(
    response, time_h, dose_median, ci_95, p_low, p_mod, p_high
  )

print(triage_table, n = Inf, width = Inf)

cat("\nDone.\n")