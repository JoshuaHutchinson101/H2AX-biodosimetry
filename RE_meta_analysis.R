library(readxl)
library(dplyr)
library(janitor)
library(ggplot2)
library(tidyr)

# ----------------------------------------------------------
# 1. Load and clean
# ----------------------------------------------------------
gamma <- read_excel(
  "C:/Users/josh1/Downloads/Calibration curves (9).xlsx",
  sheet = "gamma-H2AX"
) |>
  clean_names()

gamma$time_h <- as.numeric(gamma$time_h)  

# -----------------------------------------------------------
# 2. DL random-effects function
# ----------------------------------------------------------
dl_re <- function(y, se) {
  ok <- is.finite(y) & is.finite(se)
  y  <- y[ok]
  se <- se[ok]
  
  k <- length(y)
  if (k == 0) {
    return(list(mu = NA_real_, se_mu = NA_real_, tau2 = NA_real_, k = k))
  }
  if (k == 1) {
    return(list(mu = y[1], se_mu = se[1], tau2 = 0, k = k))
  }
  
  w <- 1 / (se^2)
  mu_fe <- sum(w * y) / sum(w)
  Q  <- sum(w * (y - mu_fe)^2)
  df <- k - 1
  c_val <- sum(w) - sum(w^2) / sum(w)
  tau2 <- ifelse(c_val > 0, max(0, (Q - df) / c_val), 0)
  
  w_re  <- 1 / (se^2 + tau2)
  mu_re <- sum(w_re * y) / sum(w_re)
  se_mu <- sqrt(1 / sum(w_re))
  
  list(mu = mu_re, se_mu = se_mu, tau2 = tau2, k = k)
}

# ----------------------------------------------------------
# 3. Pool by time_h using DL RE
# ----------------------------------------------------------
pooled_re <- gamma |>
  group_by(time_h) |>
  summarise(
    {
      res_a <- dl_re(a, se_a)
      res_b <- dl_re(b, se_b)
      
      tibble(
        a_re    = res_a$mu,
        se_a_re = res_a$se_mu,
        tau2_a  = res_a$tau2,
        n_a     = res_a$k,
        
        b_re    = res_b$mu,
        se_b_re = res_b$se_mu,
        tau2_b  = res_b$tau2,
        n_b     = res_b$k
      )
    },
    .groups = "drop"
  ) |>
  arrange(time_h)

cat("Pooled RE (with NA row):\n")
print(pooled_re)

# drop the NA-time row
pooled_re <- pooled_re |> filter(!is.na(time_h))

cat("\nPooled RE after dropping NA time:\n")
print(pooled_re)

# ----------------------------------------------------------
# 4. Build calibration curves (RE pooled)
# ----------------------------------------------------------
dose_grid <- seq(0, 5, length.out = 100)

cal_data_re <- pooled_re |>
  dplyr::select(time_h, a_re, b_re) |>       
  tidyr::expand_grid(dose = dose_grid) |>
  mutate(response = a_re + b_re * dose)

# ----------------------------------------------------------
# 5. Plot calibration curves
# ----------------------------------------------------------
ggplot(cal_data_re, aes(x = dose, y = response, colour = factor(time_h))) +
  geom_line(linewidth = 1.1) +
  labs(
    title = "Gamma-H2AX calibration curves\nRandom-effects (DL) pooled by time",
    x = "Dose (Gy)",
    y = "Predicted response",
    colour = "Time (h)"
  ) +
  theme_minimal()
