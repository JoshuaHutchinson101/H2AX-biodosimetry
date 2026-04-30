library(readxl)
library(dplyr)
library(janitor)
library(tidyr)
library(Matrix)
library(metafor)
library(ggplot2)

# ----------------------------------------------------------
# 1. Load and clean
# ----------------------------------------------------------
gamma <- read_excel(
  "C:/Users/josh1/Downloads/Calibration curves (9).xlsx",
  sheet = "gamma-H2AX"
) |>
  clean_names() |>
  dplyr::mutate(
    time_h   = as.numeric(time_h),
    lab      = factor(lab_number),
    curve_id = dplyr::row_number()
  )

# ----------------------------------------------------------
# 2. Within-curve covariance for (a_i, b_i)
# ----------------------------------------------------------
rho_ab <- -0.5  # assumed Corr(a,b) 
gamma <- gamma |>
  dplyr::mutate(cov_ab_mv = rho_ab * se_a * se_b)

# ----------------------------------------------------------
# 3. Multivariate RE meta-analysis per time_h (clustered by lab)
# ----------------------------------------------------------
fit_mv_by_time <- function(dat_time) {
  
  dat_time <- dat_time |> dplyr::arrange(curve_id)
  k <- nrow(dat_time)
  
  if (k == 0) {
    return(tibble(
      time_h    = dat_time$time_h[1],
      a_mv      = NA_real_,
      b_mv      = NA_real_,
      se_a_mv   = NA_real_,
      se_b_mv   = NA_real_,
      cov_ab_mv = NA_real_
    ))
  }
  
  long <- dat_time |>
    dplyr::select(curve_id, lab, a, se_a, b, se_b, cov_ab_mv) |>
    tidyr::pivot_longer(
      cols      = c(a, b),
      names_to  = "param",
      values_to = "yi"
    ) |>
    dplyr::mutate(param = factor(param, levels = c("a", "b")))
  
  blocks <- lapply(seq_len(k), function(i) {
    sa2 <- dat_time$se_a[i]^2
    sb2 <- dat_time$se_b[i]^2
    cab <- dat_time$cov_ab_mv[i]
    matrix(c(sa2, cab,
             cab, sb2), nrow = 2, ncol = 2)
  })
  V <- Matrix::bdiag(blocks)
  
  res_mv <- metafor::rma.mv(
    yi     = long$yi,
    V      = V,
    mods   = ~ param - 1,    # pooled a and pooled b
    random = ~ param | lab,  # lab-clustered RE; UN between-lab cov
    struct = "UN",
    data   = long
  )
  
  beta <- coef(res_mv)
  vc   <- vcov(res_mv)[1:2, 1:2]
  
  tibble(
    time_h    = dat_time$time_h[1],
    a_mv      = unname(beta["parama"]),
    b_mv      = unname(beta["paramb"]),
    se_a_mv   = sqrt(vc[1, 1]),
    se_b_mv   = sqrt(vc[2, 2]),
    cov_ab_mv = vc[1, 2]
  )
}

# ----------------------------------------------------------
# 4. Pool by time_h
# ----------------------------------------------------------
pooled_mv <- gamma |>
  dplyr::filter(!is.na(time_h)) |>
  dplyr::group_by(time_h) |>
  dplyr::group_modify(~ fit_mv_by_time(.x)) |>
  dplyr::ungroup() |>
  dplyr::arrange(time_h)

cat("Pooled multivariate RE by time (clustered by lab):\n")
print(pooled_mv)

# ----------------------------------------------------------
# 5. Plot calibration curves from pooled (a,b) per time
# ----------------------------------------------------------
dose_grid <- seq(0, 5, length.out = 100)

cal_data_mv <- pooled_mv |>
  dplyr::select(time_h, a_mv, b_mv)  |>   
  tidyr::expand_grid(dose = dose_grid) |>
  dplyr::mutate(response = a_mv + b_mv * dose)

ggplot(cal_data_mv, aes(x = dose, y = response, colour = factor(time_h))) +
  geom_line(linewidth = 1.1) +
  labs(
    title  = "Gamma-H2AX calibration curves\nMV random-effects pooled by time (lab-clustered)",
    x      = "Dose (Gy)",
    y      = "Predicted response",
    colour = "Time (h)"
  ) +
  theme_minimal()
