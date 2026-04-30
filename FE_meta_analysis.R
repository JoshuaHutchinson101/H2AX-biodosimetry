## ----------------------------------------------------------
## 0. Load packages
## ----------------------------------------------------------
library(readxl)
library(dplyr)
library(janitor)
library(ggplot2)
library(tidyr)

## ----------------------------------------------------------
## 1. Read and clean
## ----------------------------------------------------------
gamma <- read_excel(
  "C:/Users/josh1/Downloads/Calibration curves (9).xlsx", #change to users specific file
  sheet = "gamma-H2AX"
) |>
  clean_names()

gamma$time_h <- as.numeric(gamma$time_h)

## ----------------------------------------------------------
## 2. Keep required columns
## ----------------------------------------------------------
gamma_sub <- gamma |>  
  dplyr::select(time_h, a, b, se_a, se_b)  |>   
  mutate(time_h = as.numeric(time_h)) |>
  filter(!is.na(se_a), !is.na(se_b)) |> #discard if missing SE
  filter(!is.na(a), !is.na(b)) #discard missing values
## ----------------------------------------------------------
## 3. Fixed-effect pooling
## ----------------------------------------------------------
pooled_fe <- gamma_sub |>
  group_by(time_h) |>
  summarise(
    a_fe    = sum(a / se_a^2) / sum(1 / se_a^2),
    se_a_fe = 1 / sqrt(sum(1 / se_a^2)),
    b_fe    = sum(b / se_b^2) / sum(1 / se_b^2),
    se_b_fe = 1 / sqrt(sum(1 / se_b^2)),
    n       = n(),
    .groups = "drop"
  ) |>
  arrange(time_h)

print(pooled_fe)

## ----------------------------------------------------------
## 4. Build calibration curves
## ----------------------------------------------------------
dose_grid <- seq(0, 5, length.out = 100)

cal_data <- pooled_fe |>
  dplyr::select(time_h, a_fe, b_fe) |>        
  tidyr::expand_grid(dose = dose_grid) |>
  mutate(response = a_fe + b_fe * dose)

## ----------------------------------------------------------
## 5. Plot
## ----------------------------------------------------------
ggplot(cal_data, aes(x = dose, y = response, colour = factor(time_h))) +
  geom_line(linewidth = 1) +
  labs(
    title = "Gamma-H2AX calibration curves\n(Fixed-effect pooling by time)",
    x = "Dose (Gy)",
    y = "Predicted response",
    colour = "Time (h)"
  ) +
  theme_minimal()


