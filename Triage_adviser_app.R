library(shiny)
library(ggplot2)
library(dplyr)
library(MASS)
library(tibble)
library(Matrix)
library(purrr)

# ==========================================================
# SHINY APP: BIOLOGICAL DOSIMETRY TRIAGE Adviser
#
# Sources dose_1_setup.R for all model functions and
# constants so the app backend is identical to the
# analysis pipeline.
#
# Requires in environment before running:
#   pooled_mv, fit_a_pt, fit_b_fixed
#   (from 01_meta_analysis.R and 05_biexp_model.R)
# ==========================================================

# Source setup script — loads all functions, constants,
# and builds SIGMA_JOINT / MU_JOINT
source("dose_1_setup.R")

# ----------------------------------------------------------
# Pre-compute bootstrap parameter draws at startup.
# Uses the same 12x12 joint MVN as the analysis pipeline
# so the app is consistent with dose_2/3/4 outputs.
# ----------------------------------------------------------
B_PRECOMPUTE <- 10000

message("Pre-computing ", B_PRECOMPUTE, " bootstrap parameter draws...")

alpha_params <- matrix(NA_real_, nrow = B_PRECOMPUTE, ncol = 3,
                       dimnames = list(NULL, c("p0", "p1", "k1")))
beta_params  <- matrix(NA_real_, nrow = B_PRECOMPUTE, ncol = 4,
                       dimnames = list(NULL, c("p1", "p2", "k1", "k2")))

# Draw all B rows from the joint 12x12 MVN at once
joint_draws_pre <- MASS::mvrnorm(B_PRECOMPUTE,
                                 mu    = MU_JOINT,
                                 Sigma = SIGMA_JOINT)

n_conv <- 0L
for (i in seq_len(B_PRECOMPUTE)) {
  tryCatch({
    a_star <- joint_draws_pre[i, 1:N_TP]
    b_star <- joint_draws_pre[i, (N_TP + 1):(2 * N_TP)]
    
    fa <- fit_alpha_fn(pooled_mv$time_h, a_star, pooled_mv$se_a_mv)
    fb <- fit_beta_fn( pooled_mv$time_h, b_star, pooled_mv$se_b_mv)
    
    alpha_params[i, ] <- coef(fa)
    beta_params[i, ]  <- coef(fb)
    n_conv            <- n_conv + 1L
  }, error = function(e) NULL)
}

valid_idx      <- which(!is.na(alpha_params[, 1]) &
                          !is.na(beta_params[, 1]))
alpha_params_v <- alpha_params[valid_idx, , drop = FALSE]
beta_params_v  <- beta_params[valid_idx,  , drop = FALSE]
N_VALID        <- length(valid_idx)
conv_pct       <- round(100 * N_VALID / B_PRECOMPUTE, 1)

message(sprintf("Done: %d / %d draws converged (%.1f%%)",
                N_VALID, B_PRECOMPUTE, conv_pct))

# ----------------------------------------------------------
# Query function
# Pre-computed parameter matrices handle Source 1.
# Sources 2 and 3 applied fresh per query.
# ----------------------------------------------------------
query_mc <- function(response_obs, time_obs,
                     n_cells = NA,
                     phi     = PHI_MEAS,
                     sigma_t = SIGMA_T) {
  
  N <- N_VALID
  
  # Source 2: Poisson measurement noise
  if (!is.na(n_cells) && n_cells > 0) {
    se_obs  <- sqrt(phi * max(response_obs, 0.01) / n_cells)
    r_draws <- pmax(rnorm(N, mean = response_obs, sd = se_obs), 0)
  } else {
    se_obs  <- 0
    r_draws <- rep(response_obs, N)
  }
  
  # Source 3: time uncertainty
  if (sigma_t > 0) {
    t_draws <- pmin(pmax(rnorm(N, mean = time_obs, sd = sigma_t),
                         0.5), 24)
  } else {
    t_draws <- rep(time_obs, N)
  }
  
  # Evaluate decay curves at perturbed times
  a_vals <- alpha_params_v[, "p0"] +
    alpha_params_v[, "p1"] *
    exp(-alpha_params_v[, "k1"] * t_draws)
  
  b_vals <- beta_params_v[, "p1"] *
    exp(-beta_params_v[, "k1"] * t_draws) +
    beta_params_v[, "p2"] *
    exp(-beta_params_v[, "k2"] * t_draws)
  
  # Invert
  pos_b  <- b_vals > 0
  d_vals <- (r_draws[pos_b] - a_vals[pos_b]) / b_vals[pos_b]
  d_vals <- d_vals[d_vals > 0]
  
  if (length(d_vals) < 10) {
    return(list(
      dose_pt = NA, dose_median = NA, lower = NA, upper = NA,
      p_low = NA, p_mod = NA, p_high = NA,
      se_obs = se_obs, sigma_t = sigma_t,
      n_used = length(d_vals), draws = d_vals
    ))
  }
  
  ci <- hdi(d_vals, prob = 0.95)
  
  # Point estimate from locked coefficients
  ca   <- as.list(coef(fit_a_pt))
  cb   <- as.list(coef(fit_b_fixed))
  a_pt <- alpha_model(time_obs, ca$p0, ca$p1, ca$k1)
  b_pt <- beta_model( time_obs, cb$p1, cb$p2, cb$k1, cb$k2)
  d_pt <- (response_obs - a_pt) / b_pt
  
  list(
    dose_pt     = d_pt,
    dose_median = median(d_vals),
    lower       = ci[["lower"]],
    upper       = ci[["upper"]],
    p_low       = mean(d_vals < 1),
    p_mod       = mean(d_vals >= 1 & d_vals < 2),
    p_high      = mean(d_vals >= 2),
    se_obs      = se_obs,
    sigma_t     = sigma_t,
    n_used      = length(d_vals),
    draws       = d_vals
  )
}

# ----------------------------------------------------------
# Decay curve for plot
# ----------------------------------------------------------
compute_decay_curve <- function(dose_gy,
                                times = seq(0.5, 24,
                                            length.out = 300)) {
  ca <- as.list(coef(fit_a_pt))
  cb <- as.list(coef(fit_b_fixed))
  tibble(
    time  = times,
    y_val = alpha_model(times, ca$p0, ca$p1, ca$k1) +
      dose_gy * beta_model(times, cb$p1, cb$p2,
                           cb$k1, cb$k2)
  )
}

# ==========================================================
# CSS
# ==========================================================
custom_css <- "
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:wght@300;400;500;600&family=DM+Mono:wght@400;500&display=swap');

:root {
  --navy:      #0d1b2a;
  --navy-mid:  #162032;
  --navy-card: #1c2a3a;
  --teal:      #1abc9c;
  --amber:     #e6a817;
  --red:       #e74c3c;
  --green:     #27ae60;
  --text:      #e8edf2;
  --text-muted:#8fa3b8;
  --border:    #263548;
}
html, body {
  background-color: var(--navy) !important;
  color: var(--text) !important;
  font-family: 'DM Sans', sans-serif !important;
  font-size: 14px;
}
.top-header {
  background: linear-gradient(135deg, var(--navy-mid) 0%, #0f2237 100%);
  border-bottom: 2px solid var(--teal);
  padding: 18px 28px 14px;
}
.top-header h1 {
  font-size: 1.55rem; font-weight: 600; color: var(--text);
  margin: 0 0 3px; letter-spacing: -0.3px;
}
.top-header .subtitle {
  font-size: 0.8rem; color: var(--text-muted);
  font-family: 'DM Mono', monospace; letter-spacing: 0.5px;
}
.teal-dot { color: var(--teal); }
.main-wrap { display: flex; gap: 0; min-height: calc(100vh - 80px); }
.left-panel {
  width: 290px; min-width: 290px;
  background: var(--navy-mid);
  border-right: 1px solid var(--border);
  padding: 20px 18px; overflow-y: auto;
}
.right-panel { flex: 1; padding: 20px 24px; overflow-y: auto; }
.section-label {
  font-size: 0.68rem; font-weight: 600; letter-spacing: 1.2px;
  text-transform: uppercase; color: var(--teal);
  margin: 18px 0 8px; padding-bottom: 5px;
  border-bottom: 1px solid var(--border);
}
.section-label:first-child { margin-top: 0; }
.form-control, .form-control:focus {
  background: var(--navy-card) !important;
  border: 1px solid var(--border) !important;
  color: var(--text) !important;
  border-radius: 6px !important;
  font-family: 'DM Mono', monospace !important;
  font-size: 0.92rem !important;
  box-shadow: none !important;
}
.form-control:focus { border-color: var(--teal) !important; }
label {
  color: var(--text-muted) !important; font-size: 0.82rem !important;
  font-weight: 400 !important; margin-bottom: 4px !important;
}
.checkbox label { color: var(--text-muted) !important; font-size: 0.82rem !important; }
input[type='checkbox'] { accent-color: var(--teal); }
.irs--shiny .irs-bar, .irs--shiny .irs-bar-edge {
  background: var(--teal) !important; border-color: var(--teal) !important;
}
.irs--shiny .irs-handle { border-color: var(--teal) !important; }
.irs--shiny .irs-from, .irs--shiny .irs-to, .irs--shiny .irs-single {
  background: var(--teal) !important;
}
.irs--shiny .irs-line { background: var(--border) !important; }
.irs--shiny .irs-grid-text, .irs--shiny .irs-min,
.irs--shiny .irs-max { color: var(--text-muted) !important; }
.panel-divider { border: none; border-top: 1px solid var(--border); margin: 16px 0; }
.dose-row { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 14px; margin-bottom: 20px; }
.dose-card {
  background: var(--navy-card); border: 1px solid var(--border);
  border-radius: 10px; padding: 16px 18px;
  position: relative; overflow: hidden;
}
.dose-card::before {
  content: ''; position: absolute; top: 0; left: 0; right: 0; height: 3px;
}
.dose-card.pt::before  { background: var(--teal); }
.dose-card.med::before { background: var(--amber); }
.dose-card.ci::before  { background: var(--text-muted); }
.dose-card .card-label {
  font-size: 0.7rem; text-transform: uppercase;
  letter-spacing: 1px; color: var(--text-muted); margin-bottom: 6px;
}
.dose-card .card-value {
  font-family: 'DM Mono', monospace; font-size: 1.9rem;
  font-weight: 500; line-height: 1; margin-bottom: 4px;
}
.dose-card.pt  .card-value { color: var(--teal); }
.dose-card.med .card-value { color: var(--amber); }
.dose-card.ci  .card-value { color: var(--text); font-size: 1.1rem; margin-top: 6px; }
.dose-card .card-sub { font-size: 0.72rem; color: var(--text-muted); }
.triage-row { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 10px; margin-bottom: 20px; }
.triage-badge {
  border-radius: 8px; padding: 12px 10px;
  text-align: center; border: 1px solid transparent;
}
.triage-badge .badge-pct {
  font-family: 'DM Mono', monospace; font-size: 1.5rem;
  font-weight: 500; line-height: 1; margin-bottom: 3px;
}
.triage-badge .badge-label {
  font-size: 0.68rem; text-transform: uppercase; letter-spacing: 0.8px;
}
.triage-badge.low      { background: rgba(39,174,96,0.12);  border-color: rgba(39,174,96,0.3);  }
.triage-badge.low      .badge-pct,
.triage-badge.low      .badge-label { color: #27ae60; }
.triage-badge.moderate { background: rgba(230,168,23,0.12); border-color: rgba(230,168,23,0.3); }
.triage-badge.moderate .badge-pct,
.triage-badge.moderate .badge-label { color: var(--amber); }
.triage-badge.high     { background: rgba(231,76,60,0.12);  border-color: rgba(231,76,60,0.3);  }
.triage-badge.high     .badge-pct,
.triage-badge.high     .badge-label { color: var(--red); }
.chart-row { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; margin-bottom: 16px; }
.chart-card {
  background: var(--navy-card); border: 1px solid var(--border);
  border-radius: 10px; padding: 14px 16px 10px;
}
.chart-card .chart-title {
  font-size: 0.72rem; text-transform: uppercase;
  letter-spacing: 1px; color: var(--text-muted); margin-bottom: 10px;
}
.info-bar {
  background: var(--navy-card); border: 1px solid var(--border);
  border-radius: 8px; padding: 10px 16px;
  font-size: 0.75rem; color: var(--text-muted);
  font-family: 'DM Mono', monospace; margin-top: 4px;
}
.se-pill {
  display: inline-block; background: rgba(26,188,156,0.1);
  border: 1px solid rgba(26,188,156,0.3); border-radius: 4px;
  padding: 3px 8px; font-family: 'DM Mono', monospace;
  font-size: 0.78rem; color: var(--teal); margin-top: 4px;
}
.shiny-plot-output { background: transparent !important; }
"

# ==========================================================
# UI
# ==========================================================
ui <- fluidPage(
  tags$head(
    tags$style(HTML(custom_css)),
    tags$meta(name  = "viewport",
              content = "width=device-width, initial-scale=1")
  ),
  
  div(class = "top-header",
      h1(HTML(paste0(
        '<span class="teal-dot">\u03b3</span>',
        '-H2AX Biological Dosimetry Triage Adviser'
      ))),
      div(class = "subtitle",
          HTML(paste0(
            "Multivariate RE meta-analysis \u00b7 ",
            "Bi-exponential kinetic model \u00b7 ",
            "B\u2009=\u2009", N_VALID, " bootstrap draws \u00b7 ",
            "Convergence: ", conv_pct, "% \u00b7 ",
            "\u03c1_c\u2009=\u20090.8 \u00b7 95% HDI"
          ))
      )
  ),
  
  div(class = "main-wrap",
      
      # Left panel — inputs
      div(class = "left-panel",
          
          div(class = "section-label", "Patient Measurement"),
          numericInput("resp", "Measured Foci / Cell",
                       value = 4.5, min = 0, step = 0.1),
          numericInput("time", "Time Post-Exposure (h)",
                       value = 6.0, min = 0.5, max = 24, step = 0.1),
          
          hr(class = "panel-divider"),
          div(class = "section-label", "Measurement Uncertainty"),
          numericInput("n_cells", "Cells Scored",
                       value = 50, min = 1, step = 10),
          uiOutput("se_obs_display"),
          checkboxInput(
            "use_meas_unc",
            HTML(paste0(
              "Apply measurement uncertainty<br>",
              "<span style='color:#637d94;font-size:0.78rem'>",
              "R* ~ N(foci, SE\u00b2) \u00b7 Poisson SE</span>"
            )),
            value = TRUE
          ),
          
          hr(class = "panel-divider"),
          div(class = "section-label", "Time Uncertainty"),
          sliderInput("sigma_t", "Time SD (hours)",
                      min = 0, max = 2, value = 0.5, step = 0.1),
          checkboxInput(
            "use_time_unc",
            HTML(paste0(
              "Apply time uncertainty<br>",
              "<span style='color:#637d94;font-size:0.78rem'>",
              "t* ~ N(t, \u03c3\u1d57\u00b2) \u00b7 default \u00b130\u2009min</span>"
            )),
            value = TRUE
          ),
          
          hr(class = "panel-divider"),
          div(class = "section-label", "Active Uncertainty Sources"),
          uiOutput("unc_sources_summary")
      ),
      
      # Right panel — outputs
      div(class = "right-panel",
          
          div(class = "dose-row",
              div(class = "dose-card pt",
                  div(class = "card-label", "Point Estimate"),
                  div(class = "card-value", textOutput("dose_pt_out")),
                  div(class = "card-sub",
                      "Gy \u00b7 exact inputs, no uncertainty")
              ),
              div(class = "dose-card med",
                  div(class = "card-label", "MC Median"),
                  div(class = "card-value", textOutput("dose_med_out")),
                  div(class = "card-sub",
                      "Gy \u00b7 all uncertainty sources")
              ),
              div(class = "dose-card ci",
                  div(class = "card-label", "95% HDI"),
                  div(class = "card-value", uiOutput("ci_value")),
                  div(class = "card-sub",   uiOutput("n_draws_note"))
              )
          ),
          
          div(class = "section-label",
              style = "margin-bottom:12px",
              "Triage Classification"),
          div(class = "triage-row",
              div(class = "triage-badge low",
                  div(class = "badge-pct", textOutput("p_low")),
                  div(class = "badge-label", "Low \u00b7 < 1\u2009Gy")
              ),
              div(class = "triage-badge moderate",
                  div(class = "badge-pct", textOutput("p_mod")),
                  div(class = "badge-label",
                      "Moderate \u00b7 1\u20132\u2009Gy")
              ),
              div(class = "triage-badge high",
                  div(class = "badge-pct", textOutput("p_high")),
                  div(class = "badge-label", "High \u00b7 > 2\u2009Gy")
              )
          ),
          
          div(class = "chart-row",
              div(class = "chart-card",
                  div(class = "chart-title",
                      "\u03b3-H2AX Decay Pathway"),
                  plotOutput("decayPlot", height = "280px")
              ),
              div(class = "chart-card",
                  div(class = "chart-title",
                      "Monte Carlo Dose Distribution"),
                  plotOutput("distPlot", height = "280px")
              )
          ),
          
          div(class = "info-bar",
              HTML(paste0(
                "Model: 4-param bi-exp \u03b2(t) + single-exp \u03b1(t) \u00b7 ",
                "Meta-analysis: multivariate RE, lab-clustered \u00b7 ",
                "Source\u20091: 12\u00d712 joint MVN (\u03c1_within=\u22120.5, ",
                "\u03c1_cross=0.8) \u00b7 ",
                "Source\u20092: Poisson SE \u00b7 ",
                "Source\u20093: Normal time perturbation \u00b7 ",
                "CI: 95% HDI (Chen\u2013Shao)"
              ))
          )
      )
  )
)

# ==========================================================
# SERVER
# ==========================================================
theme_dosim <- function() {
  theme_minimal(base_family = "sans", base_size = 12) +
    theme(
      plot.background   = element_rect(fill = "#1c2a3a", colour = NA),
      panel.background  = element_rect(fill = "#1c2a3a", colour = NA),
      panel.grid.major  = element_line(colour = "#263548", linewidth = 0.4),
      panel.grid.minor  = element_blank(),
      axis.text         = element_text(colour = "#8fa3b8", size = 9),
      axis.title        = element_text(colour = "#8fa3b8", size = 10),
      plot.title        = element_blank(),
      legend.background = element_rect(fill = "#1c2a3a", colour = NA),
      legend.text       = element_text(colour = "#8fa3b8", size = 9),
      legend.title      = element_text(colour = "#8fa3b8", size = 9),
      plot.margin       = margin(4, 8, 4, 4)
    )
}

server <- function(input, output, session) {
  
  eff_n_cells <- reactive({
    if (isTRUE(input$use_meas_unc)) input$n_cells else NA_real_
  })
  
  eff_sigma_t <- reactive({
    if (isTRUE(input$use_time_unc)) input$sigma_t else 0
  })
  
  output$se_obs_display <- renderUI({
    nc <- input$n_cells; r <- input$resp
    if (is.na(nc) || nc <= 0 || is.na(r) || r <= 0) return(NULL)
    se <- sqrt(PHI_MEAS * max(r, 0.01) / nc)
    div(class = "se-pill",
        paste0("SE\u2092\u2093\u2099 = ", round(se, 3),
               " foci/cell"))
  })
  
  output$unc_sources_summary <- renderUI({
    sources <- c(
      "\u2022 Model (12\u00d712 joint MVN)",
      if (isTRUE(input$use_meas_unc))
        paste0("\u2022 Measurement (Poisson, n=",
               input$n_cells, ")"),
      if (isTRUE(input$use_time_unc) && input$sigma_t > 0)
        paste0("\u2022 Time (\u03c3=", input$sigma_t, "h)")
    )
    tags$p(HTML(paste(sources, collapse = "<br>")),
           style = paste0("font-size:0.8rem; color:#8fa3b8;",
                          " line-height:1.6; margin:0;"))
  })
  
  result <- reactive({
    req(input$resp, input$time)
    validate(
      need(input$resp    >= 0,   "Foci/Cell must be \u2265 0"),
      need(input$time    >= 0.5, "Time must be \u2265 0.5h"),
      need(input$time    <= 24,  "Time must be \u2264 24h"),
      need(input$n_cells >= 1,   "Cells scored must be \u2265 1")
    )
    suppressWarnings(query_mc(
      response_obs = input$resp,
      time_obs     = input$time,
      n_cells      = eff_n_cells(),
      phi          = PHI_MEAS,
      sigma_t      = eff_sigma_t()
    ))
  })
  
  output$dose_pt_out <- renderText({
    d <- result()$dose_pt
    if (is.na(d) || d < 0) "< 0" else sprintf("%.2f", d)
  })
  
  output$dose_med_out <- renderText({
    d <- result()$dose_median
    if (is.na(d)) "N/A" else sprintf("%.2f", d)
  })
  
  output$ci_value <- renderUI({
    r <- result()
    if (is.na(r$lower))
      return(div(style = "color:#8fa3b8; font-size:0.9rem",
                 "Insufficient draws"))
    HTML(sprintf("[%.2f \u2013 %.2f]\u2009Gy",
                 max(0, r$lower), r$upper))
  })
  
  output$n_draws_note <- renderUI({
    HTML(paste0(result()$n_used, " / ", N_VALID, " draws"))
  })
  
  fmt_pct <- function(x) {
    if (is.na(x)) "\u2014" else
      paste0(round(x * 100, 1), "%")
  }
  output$p_low  <- renderText({ fmt_pct(result()$p_low)  })
  output$p_mod  <- renderText({ fmt_pct(result()$p_mod)  })
  output$p_high <- renderText({ fmt_pct(result()$p_high) })
  
  output$decayPlot <- renderPlot({
    r     <- result()
    d_est <- if (is.na(r$dose_pt) || r$dose_pt < 0) 0
    else r$dose_pt
    curve <- compute_decay_curve(d_est)
    
    ci_str <- if (!is.na(r$lower))
      paste0("95% HDI [",
             round(max(0, r$lower), 2), "\u2013",
             round(r$upper, 2), "] Gy")
    else ""
    
    p <- ggplot(curve, aes(x = time, y = y_val)) +
      geom_line(colour = "#1abc9c", linewidth = 1.4)
    
    if (eff_sigma_t() > 0) {
      p <- p + annotate("rect",
                        xmin = max(0.5, input$time - eff_sigma_t()),
                        xmax = min(24,  input$time + eff_sigma_t()),
                        ymin = -Inf, ymax = Inf,
                        fill = "#e74c3c", alpha = 0.07)
    }
    
    p +
      geom_vline(xintercept = input$time,
                 linetype = "dashed", colour = "#e74c3c",
                 linewidth = 0.5) +
      geom_point(aes(x = input$time, y = input$resp),
                 colour = "#e74c3c", size = 5, shape = 16) +
      annotate("label",
               x     = pmin(input$time + 1.5, 21),
               y     = max(curve$y_val, na.rm = TRUE) * 0.88,
               label = paste0(round(d_est, 2), " Gy\n", ci_str),
               fill  = "#162032", colour = "#e8edf2",
               size  = 3, label.size = 0.3, lineheight = 1.3) +
      labs(x = "Hours post-exposure",
           y = "Predicted Foci / Cell") +
      theme_dosim()
  }, bg = "#1c2a3a")
  
  output$distPlot <- renderPlot({
    r <- result()
    req(!is.null(r$draws), length(r$draws) >= 10)
    df <- data.frame(dose = r$draws)
    
    ggplot(df, aes(x = dose)) +
      geom_histogram(
        aes(fill = cut(dose,
                       breaks = c(-Inf, 1, 2, Inf),
                       labels = c("Low", "Moderate", "High"))),
        bins = 55, colour = "#1c2a3a", linewidth = 0.15
      ) +
      scale_fill_manual(
        values = c("Low"      = "#27ae60",
                   "Moderate" = "#e6a817",
                   "High"     = "#e74c3c"),
        name = NULL, drop = FALSE
      ) +
      geom_vline(xintercept = r$dose_median,
                 colour = "#e8edf2", linewidth = 1) +
      geom_vline(xintercept = c(r$lower, r$upper),
                 colour = "#8fa3b8", linewidth = 0.6,
                 linetype = "dashed") +
      annotate("text",
               x = r$dose_median, y = Inf,
               label = paste0("med=", round(r$dose_median, 2)),
               vjust = 2, hjust = -0.1,
               colour = "#e8edf2", size = 3) +
      labs(x = "Estimated Dose (Gy)", y = "Count") +
      theme_dosim() +
      theme(legend.position = "top",
            legend.key.size  = unit(0.45, "cm"))
  }, bg = "#1c2a3a")
}

# ==========================================================
shinyApp(ui, server)
