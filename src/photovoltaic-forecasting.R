# =============================================================================
# TFG: Photovoltaic Installation Analysis — Navantia Module Workshop
# Multivariate and Functional Data Analysis
#
# Author : Hugo Balado Mosquera
# Degree : B.Sc. Data Science & Engineering — Universidade da Coruña
# Grade  : Honours (Matrícula de Honor)
# =============================================================================

library(fda.usc)
library(tidyverse)
library(lubridate)
library(data.table)
library(energy)
library(dbscan)
library(psych)
library(cluster)
library(ranger)
library(caret)
library(corrplot)
library(MASS)
library(car)
library(leaps)
library(mctest)
library(ggm)
library(Hmisc)
library(lmtest)
library(tseries)
library(nortest)

select <- dplyr::select  # avoid MASS::select masking

# =============================================================================
# PATHS  —  edit these to match your local data directory
# =============================================================================

DATA_DIR      <- "data/"
PATH_INV1     <- file.path(DATA_DIR, "FV_Rep1.csv")
PATH_INV2     <- file.path(DATA_DIR, "FV_Rep2.csv")
PATH_INV3     <- file.path(DATA_DIR, "FV_Rep3.csv")
PATH_METEO    <- file.path(DATA_DIR, "METEO_MAESTRA_TFG_2024_2026.csv")
PATH_OUTPUTS  <- "outputs/"
dir.create(PATH_OUTPUTS, showWarnings = FALSE)

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

COL_CLUSTER <- c("1" = "#CD5C5C", "2" = "#FDAE61", "3" = "#4682B4")
COL_LINEA   <- "#34495E"
COL_PUNTO   <- "#2C3E50"
COL_MES     <- setNames(
  colorRampPalette(c("#4682B4","#92C5DE","#F4A582","#CD5C5C","#B2182B"))(12),
  c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
)

FESTIVOS_GALICIA <- as.Date(c(
  "2024-11-01","2024-12-06","2024-12-25",
  "2025-01-01","2025-01-06","2025-04-17","2025-04-18","2025-05-01",
  "2025-05-17","2025-07-25","2025-08-15","2025-11-01",
  "2025-12-06","2025-12-08","2025-12-25",
  "2026-01-01","2026-01-06","2026-02-17"
))

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Parse a Navantia SCADA sparse CSV.
#' Header: 3 lines. Separator: ";". Date format: "YYYY.MM.DD HH:MM:SS.sss".
parsear_scada_sparse <- function(path, nf_target = 2, col_name = "value") {
  lines  <- readLines(path, encoding = "UTF-8")
  lines  <- lines[-(1:3)]
  fields <- strsplit(lines, ";")
  ok     <- fields[lengths(fields) == nf_target]

  if (length(ok) == 0) {
    warning("No rows with NF=", nf_target, " in ", basename(path))
    return(tibble(Fecha = as.POSIXct(character()), !!col_name := numeric()))
  }

  tibble(
    Fecha     = with_tz(
      as.POSIXct(vapply(ok, `[[`, character(1), 1),
                 format = "%Y.%m.%d %H:%M:%OS", tz = "Etc/GMT-2"),
      tzone = "UTC"
    ),
    !!col_name := as.numeric(gsub(",", ".", vapply(ok, `[[`, character(1), nf_target)))
  ) |>
    filter(!is.na(Fecha), !is.na(.data[[col_name]])) |>
    arrange(Fecha)
}

#' Regularise async data to a 10-minute grid (mean or last value).
regularizar_10min <- function(df, method = c("mean", "last")) {
  method <- match.arg(method)
  df |>
    mutate(Fecha = floor_date(Fecha, "10 minutes")) |>
    group_by(Fecha) |>
    summarise(
      across(everything(),
             if (method == "mean") ~ mean(., na.rm = TRUE) else ~ last(.)),
      .groups = "drop"
    ) |>
    arrange(Fecha)
}

#' Reorder cluster labels so Cluster 1 always has the highest mean of
#' `col_order`. Ensures consistent colour coding across analyses.
reordenar_clusters <- function(df, col_cluster, col_order = "Pot_Media") {
  medias <- df |>
    filter(!is.na(.data[[col_order]])) |>
    group_by(across(all_of(col_cluster))) |>
    summarise(mean_ord = mean(.data[[col_order]], na.rm = TRUE), .groups = "drop") |>
    arrange(desc(mean_ord)) |>
    mutate(new_label = as.factor(row_number()))

  mapping <- setNames(as.character(medias$new_label),
                      as.character(medias[[col_cluster]]))
  df[[col_cluster]] <- factor(mapping[as.character(df[[col_cluster]])])
  df
}

# =============================================================================
# STEP 1 — IMPORT AND CLEANING
# =============================================================================

df_inv1 <- parsear_scada_sparse(PATH_INV1, 2, "Potencia_AC_Inv1")
df_inv2 <- parsear_scada_sparse(PATH_INV2, 2, "Potencia_AC_Inv2")
df_inv3 <- parsear_scada_sparse(PATH_INV3, 2, "Potencia_AC_Inv3")

cat("Inverter 1:", nrow(df_inv1), "records |",
    format(min(df_inv1$Fecha), "%Y-%m-%d"), "to",
    format(max(df_inv1$Fecha), "%Y-%m-%d"), "\n")
cat("Inverter 2:", nrow(df_inv2), "records |",
    format(min(df_inv2$Fecha), "%Y-%m-%d"), "to",
    format(max(df_inv2$Fecha), "%Y-%m-%d"), "\n")
cat("Inverter 3:", nrow(df_inv3), "records |",
    format(min(df_inv3$Fecha), "%Y-%m-%d"), "to",
    format(max(df_inv3$Fecha), "%Y-%m-%d"), "\n\n")

# MeteoGalicia encodes missing values as -9999.
# Note: METEO_MAESTRA has a double-conversion UTC offset bug that is corrected
# by reading as UTC, re-interpreting as Europe/Madrid, then forcing back to UTC.
df_meteo <- fread(PATH_METEO, sep = ",", header = TRUE) |>
  as_tibble() |>
  mutate(
    Fecha = ymd_hms(Fecha, tz = "UTC"),
    Fecha = with_tz(Fecha, tzone = "Europe/Madrid"),
    Fecha = force_tz(Fecha, tzone = "UTC")
  ) |>
  mutate(across(where(is.numeric), ~ if_else(. <= -9999, NA_real_, .))) |>
  filter(!is.na(Fecha)) |>
  arrange(Fecha)

cat("Meteo:", nrow(df_meteo), "records |",
    format(min(df_meteo$Fecha), "%Y-%m-%d"), "to",
    format(max(df_meteo$Fecha), "%Y-%m-%d"), "\n\n")

# =============================================================================
# STEP 2 — REGULARISE TO 10-MINUTE GRID AND BUILD FINAL TABLE
# =============================================================================

df_potencia_10min <- df_inv1 |>
  full_join(df_inv2, by = "Fecha") |>
  full_join(df_inv3, by = "Fecha") |>
  arrange(Fecha) |>
  regularizar_10min("mean") |>
  mutate(Potencia_AC_Total = rowSums(
    select(., starts_with("Potencia_AC_Inv")), na.rm = TRUE
  ))

cat("Power 10-min grid:", nrow(df_potencia_10min), "records\n\n")

df_final <- df_meteo |>
  inner_join(df_potencia_10min, by = "Fecha") |>
  mutate(
    Mes             = lubridate::month(with_tz(Fecha, "Europe/Madrid"), label = TRUE, abbr = TRUE),
    Dia_Semana      = lubridate::wday(with_tz(Fecha, "Europe/Madrid"), label = TRUE, abbr = TRUE, week_start = 1),
    Hora            = lubridate::hour(Fecha),
    Es_Finde        = lubridate::wday(with_tz(Fecha, "Europe/Madrid"), week_start = 1) >= 6,
    Es_Festivo      = as.Date(with_tz(Fecha, "Europe/Madrid")) %in% FESTIVOS_GALICIA,
    Es_No_Laborable = Es_Finde | Es_Festivo
  )

cat("Final records:", nrow(df_final), "\n")
cat("Period:", format(min(df_final$Fecha), "%Y-%m-%d"), "to",
    format(max(df_final$Fecha), "%Y-%m-%d"), "\n")
cat("Columns:", ncol(df_final), "\n\n")

fwrite(df_final, file.path(PATH_OUTPUTS, "df_final_10min.csv"), dateTimeAs = "write.csv")

# =============================================================================
# STEP 3 — DAILY AGGREGATION AND ANOMALY DETECTION
# =============================================================================
# Meteorological variables are aggregated from df_meteo (full 24h) to avoid
# the diurnal bias introduced by the inner_join in df_final.
# Power variables come from df_final (daylight hours only).

fecha_inicio_prod <- as.Date(min(df_final$Fecha),  tz = "Europe/Madrid")
fecha_fin_prod    <- as.Date(max(df_final$Fecha),  tz = "Europe/Madrid")

meteo_diario <- df_meteo |>
  mutate(Fecha_Dia = as.Date(with_tz(Fecha, "Europe/Madrid"))) |>
  filter(between(Fecha_Dia, fecha_inicio_prod, fecha_fin_prod)) |>
  group_by(Fecha_Dia) |>
  summarise(
    Rad_Media        = mean(Radiacion_Solar,  na.rm = TRUE),
    Temp_Media       = mean(Temperatura_1.5m, na.rm = TRUE),
    Humedad_Relativa = mean(Humedad_Relativa, na.rm = TRUE),
    Viento_Velocidad = mean(Viento_Velocidad, na.rm = TRUE),
    Viento_Direccion = mean(Viento_Direccion, na.rm = TRUE),
    Presion          = mean(Presion,          na.rm = TRUE),
    Radiacion_UV     = mean(Radiacion_UV,     na.rm = TRUE),
    Precipitacion    = sum(Precipitacion,     na.rm = TRUE),
    N_registros_meteo = n(),
    .groups = "drop"
  )

potencia_diaria <- df_final |>
  mutate(Fecha_Dia = as.Date(with_tz(Fecha, "Europe/Madrid"))) |>
  group_by(Fecha_Dia) |>
  summarise(
    Pot_Media        = mean(Potencia_AC_Total, na.rm = TRUE),
    Potencia_AC_Inv1 = mean(Potencia_AC_Inv1, na.rm = TRUE),
    Potencia_AC_Inv2 = mean(Potencia_AC_Inv2, na.rm = TRUE),
    Potencia_AC_Inv3 = mean(Potencia_AC_Inv3, na.rm = TRUE),
    N_registros_prod = n(),
    Mes              = first(Mes),
    Dia_Semana       = first(Dia_Semana),
    Es_No_Laborable  = first(Es_No_Laborable),
    Es_Festivo       = first(Es_Festivo),
    .groups = "drop"
  )

df_daily <- meteo_diario |>
  left_join(potencia_diaria, by = "Fecha_Dia") |>
  filter(!is.na(Rad_Media)) |>
  mutate(
    Mes             = lubridate::month(Fecha_Dia, label = TRUE, abbr = TRUE),
    Dia_Semana      = lubridate::wday(Fecha_Dia, label = TRUE, abbr = TRUE, week_start = 1),
    Es_Finde        = lubridate::wday(Fecha_Dia, week_start = 1) >= 6,
    Es_Festivo      = Fecha_Dia %in% FESTIVOS_GALICIA,
    Es_No_Laborable = Es_Finde | Es_Festivo
  )

# Anomaly detection: solar days with near-zero power (technical shutdowns).
df_daily <- df_daily |>
  mutate(Es_Anomalo = !is.na(Pot_Media) & (Pot_Media < 1 & Rad_Media > 50))

dias_anomalos <- filter(df_daily, Es_Anomalo)
cat("Anomalous days detected:", nrow(dias_anomalos), "\n")

df_daily_clean <- df_daily |>
  filter(!Es_Anomalo, !is.na(Pot_Media)) |>
  select(-Radiacion_UV, -Viento_Direccion, -N_registros_meteo, -N_registros_prod)

cat("Valid days for analysis:", nrow(df_daily_clean), "\n\n")

fwrite(df_daily,       file.path(PATH_OUTPUTS, "df_final_daily.csv"),   dateTimeAs = "write.csv")
fwrite(df_daily_clean, file.path(PATH_OUTPUTS, "df_daily_clean.csv"),   dateTimeAs = "write.csv")

# =============================================================================
# STEP 4 — UNIVARIATE EDA
# =============================================================================

# --- 4A. % NAs in meteorological data ---
pct_na_meteo <- df_meteo |>
  select(where(is.numeric)) |>
  summarise(across(everything(), ~ round(100 * mean(is.na(.)), 2))) |>
  pivot_longer(everything(), names_to = "Variable", values_to = "Pct_NA") |>
  arrange(desc(Pct_NA))

cat("% NAs in meteorological data:\n"); print(pct_na_meteo, n = 20)

# --- 4B. Numerical summary ---
vars_summary <- c("Pot_Media","Rad_Media","Temp_Media","Humedad_Relativa",
                  "Precipitacion","Viento_Velocidad","Presion")

tabla_resumen <- df_daily_clean |>
  select(all_of(vars_summary)) |>
  pivot_longer(everything(), names_to = "Variable", values_to = "Valor") |>
  group_by(Variable) |>
  summarise(
    n       = sum(!is.na(Valor)),
    Media   = round(mean(Valor,          na.rm = TRUE), 2),
    Mediana = round(median(Valor,        na.rm = TRUE), 2),
    DE      = round(sd(Valor,            na.rm = TRUE), 2),
    Min     = round(min(Valor,           na.rm = TRUE), 2),
    Q1      = round(quantile(Valor, .25, na.rm = TRUE), 2),
    Q3      = round(quantile(Valor, .75, na.rm = TRUE), 2),
    Max     = round(max(Valor,           na.rm = TRUE), 2),
    .groups = "drop"
  ) |>
  mutate(Variable = factor(Variable, levels = vars_summary)) |>
  arrange(Variable)

cat("\nNumerical summary (n =", nrow(df_daily_clean), "):\n")
print(tabla_resumen, width = Inf)

# --- 4C. Histograms with density ---
plot_hist <- function(data, var, label, unit, bins = 30) {
  ggplot(data, aes(x = .data[[var]])) +
    geom_histogram(aes(y = after_stat(density)), bins = bins,
                   fill = "#4682B4", color = "white", alpha = 0.7) +
    geom_density(color = "#CD5C5C", linewidth = 0.8) +
    labs(x = paste0(label, " (", unit, ")"), y = "Density") +
    theme_minimal(base_size = 10)
}

print(plot_hist(df_daily_clean, "Pot_Media",        "Mean daily AC power",       "kW"))
print(plot_hist(df_daily_clean, "Rad_Media",        "Mean daily solar radiation", "W/m²"))
print(plot_hist(df_daily_clean, "Temp_Media",       "Mean daily temperature",    "°C"))
print(plot_hist(df_daily_clean, "Humedad_Relativa", "Mean daily humidity",       "%"))
print(plot_hist(df_daily_clean, "Precipitacion",    "Daily precipitation",       "L/m²"))
print(plot_hist(df_daily_clean, "Viento_Velocidad", "Mean daily wind speed",     "km/h"))
print(plot_hist(df_daily_clean, "Presion",          "Mean daily pressure",       "hPa"))

# --- 4D. Monthly boxplots ---
MESES_ORDEN <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

df_box <- df_daily_clean |>
  mutate(Mes = factor(Mes, levels = MESES_ORDEN)) |>
  left_join(
    df_daily_clean |>
      group_by(Mes) |>
      summarise(Pot_Mediana = median(Pot_Media, na.rm = TRUE), .groups = "drop"),
    by = "Mes"
  )

print(
  ggplot(df_box, aes(x = Mes, y = Pot_Media, fill = Pot_Mediana)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1.5) +
    scale_fill_gradient(low = "#4682B4", high = "#CD5C5C", name = "Median\nPow (kW)") +
    labs(x = "Month", y = "Mean daily AC power (kW)") +
    theme_minimal(base_size = 10)
)

# =============================================================================
# STEP 5 — BIVARIATE EDA
# =============================================================================

plot_biv <- function(df, x_var, x_label, x_unit) {
  ggplot(df, aes(x = .data[[x_var]], y = Pot_Media)) +
    geom_point(alpha = 0.5, color = "#4682B4", size = 1.8) +
    geom_smooth(method = "lm", se = TRUE, color = "#CD5C5C",
                linewidth = 0.8, linetype = "dashed") +
    labs(x = paste0(x_label, " (", x_unit, ")"),
         y = "Mean daily AC power (kW)") +
    theme_minimal(base_size = 10)
}

print(plot_biv(df_daily_clean, "Rad_Media",        "Mean daily solar radiation", "W/m²"))
print(plot_biv(df_daily_clean, "Temp_Media",       "Mean daily temperature",    "°C"))
print(plot_biv(df_daily_clean, "Humedad_Relativa", "Mean daily humidity",       "%"))
print(plot_biv(df_daily_clean, "Viento_Velocidad", "Mean daily wind speed",     "km/h"))
print(plot_biv(df_daily_clean, "Precipitacion",    "Daily precipitation",       "L/m²"))
print(plot_biv(df_daily_clean, "Presion",          "Mean daily pressure",       "hPa"))

# =============================================================================
# STEP 6 — PEARSON CORRELATION MATRIX
# =============================================================================

df_cor <- df_daily_clean |>
  select(Rad_Media, Temp_Media, Humedad_Relativa, Precipitacion,
         Viento_Velocidad, Presion, Potencia_AC_Inv1, Potencia_AC_Inv2,
         Potencia_AC_Inv3, Pot_Media) |>
  rename(Hum_Rel = Humedad_Relativa, Precip = Precipitacion,
         Viento_Vel = Viento_Velocidad,
         Pot_Inv1 = Potencia_AC_Inv1, Pot_Inv2 = Potencia_AC_Inv2,
         Pot_Inv3 = Potencia_AC_Inv3)

cor_matrix <- cor(df_cor, use = "pairwise.complete.obs") |> round(2)

cat("\nPearson correlations with Pot_Media:\n")
cor_matrix["Pot_Media", ] |>
  enframe(name = "Variable", value = "r_Pearson") |>
  filter(Variable != "Pot_Media") |>
  arrange(desc(abs(r_Pearson))) |>
  print(n = 20)

# Heatmap
cor_long <- cor_matrix |>
  as.data.frame() |>
  rownames_to_column("Var1") |>
  pivot_longer(-Var1, names_to = "Var2", values_to = "Cor") |>
  mutate(Var1 = factor(Var1, levels = colnames(cor_matrix)),
         Var2 = factor(Var2, levels = rev(colnames(cor_matrix))))

print(
  ggplot(cor_long, aes(x = Var1, y = Var2, fill = Cor)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = Cor), size = 2.5) +
    scale_fill_gradient2(low = "#4682B4", mid = "white", high = "#CD5C5C",
                         midpoint = 0, limits = c(-1, 1), name = "r Pearson") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          panel.grid  = element_blank())
)

# --- Székely distance correlation (non-linear dependence) ---
vars_dcor <- c("Rad_Media","Temp_Media","Hum_Rel","Precip","Viento_Vel","Presion","Pot_Media")

datos_dcor <- df_daily_clean |>
  select(Rad_Media, Temp_Media, Humedad_Relativa, Precipitacion,
         Viento_Velocidad, Presion, Pot_Media) |>
  rename(Hum_Rel = Humedad_Relativa, Precip = Precipitacion,
         Viento_Vel = Viento_Velocidad) |>
  as.data.frame()

dcor_mat <- matrix(1, length(vars_dcor), length(vars_dcor),
                   dimnames = list(vars_dcor, vars_dcor))
for (i in 1:(length(vars_dcor) - 1)) {
  for (j in (i + 1):length(vars_dcor)) {
    v <- dcor(datos_dcor[, i], datos_dcor[, j])
    dcor_mat[i, j] <- dcor_mat[j, i] <- v
  }
}

cat("\nSzékely distance correlation matrix:\n"); print(round(dcor_mat, 2))

# =============================================================================
# STEP 7 — UNIVARIATE OUTLIER DETECTION (TUKEY IQR)
# =============================================================================

vars_outliers <- c("Pot_Media","Rad_Media","Temp_Media","Humedad_Relativa",
                   "Precipitacion","Viento_Velocidad","Presion")

outliers_univ <- map_dfr(vars_outliers, function(v) {
  x   <- df_daily_clean[[v]]
  q1  <- quantile(x, .25, na.rm = TRUE)
  q3  <- quantile(x, .75, na.rm = TRUE)
  iqr <- q3 - q1

  lim_mild_lo <- q1 - 1.5 * iqr;  lim_mild_hi <- q3 + 1.5 * iqr
  lim_ext_lo  <- q1 - 3.0 * iqr;  lim_ext_hi  <- q3 + 3.0 * iqr

  idx <- which(x < lim_mild_lo | x > lim_mild_hi)
  if (length(idx) == 0) return(tibble())

  tibble(
    Variable        = v,
    Fecha_Dia       = df_daily_clean$Fecha_Dia[idx],
    Valor           = x[idx],
    Lim_Mild_Lo     = round(lim_mild_lo, 2), Lim_Mild_Hi  = round(lim_mild_hi, 2),
    Lim_Ext_Lo      = round(lim_ext_lo,  2), Lim_Ext_Hi   = round(lim_ext_hi,  2),
    Direction       = if_else(x[idx] < lim_mild_lo, "Lower", "Upper"),
    Severity        = if_else(x[idx] < lim_ext_lo | x[idx] > lim_ext_hi, "Extreme", "Mild")
  )
})

dias_atipicos_univ <- outliers_univ |>
  group_by(Fecha_Dia) |>
  summarise(N_vars     = n_distinct(Variable),
            Max_Sev    = if_else(any(Severity == "Extreme"), "Extreme", "Mild"),
            Variables  = paste(unique(Variable), collapse = ", "),
            .groups    = "drop") |>
  arrange(desc(Max_Sev), desc(N_vars))

cat("\nUnique outlier days:", nrow(dias_atipicos_univ), "of", nrow(df_daily_clean), "\n")
cat("  Extreme:", sum(dias_atipicos_univ$Max_Sev == "Extreme"), "\n")
cat("  Mild only:", sum(dias_atipicos_univ$Max_Sev == "Mild"), "\n\n")

df_daily_clean <- df_daily_clean |>
  mutate(
    Es_Atipico_Univ = Fecha_Dia %in% dias_atipicos_univ$Fecha_Dia,
    Es_Extremo_Univ = Fecha_Dia %in% filter(dias_atipicos_univ, Max_Sev == "Extreme")$Fecha_Dia
  )

# =============================================================================
# STEP 8 — MULTIVARIATE OUTLIER DETECTION (LOF)
# =============================================================================

vars_lof <- c("Pot_Media","Rad_Media","Temp_Media","Humedad_Relativa",
              "Precipitacion","Viento_Velocidad","Presion")

mat_lof <- df_daily_clean |> select(all_of(vars_lof)) |> scale() |> as.data.frame()

K_LOF <- 10
lof_scores <- lof(mat_lof, minPts = K_LOF)
df_daily_clean$LOF_Score <- lof_scores

umbral_lof <- 2
df_daily_clean <- df_daily_clean |>
  mutate(Es_Atipico_LOF    = LOF_Score > umbral_lof,
         Es_Atipico_Alguno = Es_Atipico_Univ | Es_Atipico_LOF)

cat("LOF outliers (score >", umbral_lof, "):", sum(df_daily_clean$Es_Atipico_LOF), "\n")

# Dataset for MDS / clustering: exclude LOF outliers
df_analysis    <- filter(df_daily_clean, !Es_Atipico_LOF)
fechas_lof     <- filter(df_daily_clean, Es_Atipico_LOF)$Fecha_Dia
meteo_analysis <- filter(meteo_diario, !(Fecha_Dia %in% fechas_lof))

cat("Days for MDS/clustering:", nrow(df_analysis),
    "(removed", sum(df_daily_clean$Es_Atipico_LOF), "LOF outliers)\n\n")

fwrite(df_daily_clean, file.path(PATH_OUTPUTS, "df_daily_annotated.csv"), dateTimeAs = "write.csv")

# =============================================================================
# STEP 9 — PRINCIPAL COMPONENT ANALYSIS (PCA)
# =============================================================================

# Helper: variance table
tabla_varianza <- function(pca) {
  ev   <- pca$sdev^2
  prop <- ev / sum(ev)
  data.frame(PC         = paste0("PC", seq_along(ev)),
             Eigenvalue = round(ev,        4),
             Var_Pct    = round(100 * prop, 2),
             Var_Cum    = round(100 * cumsum(prop), 2))
}

# Helper: scree plot
plot_scree <- function(tv) {
  n <- nrow(tv)
  ggplot(tv, aes(seq_len(n), Eigenvalue)) +
    geom_line(color = COL_LINEA, linewidth = 1) +
    geom_point(color = COL_PUNTO, size = 3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "#CD5C5C", linewidth = 0.8) +
    scale_x_continuous(breaks = seq_len(n), labels = tv$PC) +
    labs(x = "Principal component", y = "Eigenvalue") +
    theme_minimal(base_size = 10)
}

# Helper: loading heatmap
plot_loadings <- function(pca, n_pc = NULL) {
  rot <- pca$rotation
  if (!is.null(n_pc)) rot <- rot[, 1:n_pc, drop = FALSE]
  as.data.frame(rot) |>
    rownames_to_column("Variable") |>
    pivot_longer(-Variable, names_to = "PC", values_to = "Loading") |>
    mutate(PC = factor(PC, levels = colnames(rot))) |>
    ggplot(aes(x = PC, y = Variable, fill = Loading)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = round(Loading, 2)), size = 2.8) +
    scale_fill_gradient2(low = "#4682B4", mid = "white", high = "#CD5C5C",
                         midpoint = 0, limits = c(-1, 1)) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank())
}

# --- 9A. PCA on inverters ---
mat_inv <- df_analysis |> select(Potencia_AC_Inv1, Potencia_AC_Inv2, Potencia_AC_Inv3) |> filter(complete.cases(.))
R_inv   <- cor(mat_inv, use = "pairwise.complete.obs")

cat("Inter-inverter correlations:\n"); print(round(R_inv, 4))
cat("Bartlett test:\n"); print(cortest.bartlett(R_inv, n = nrow(mat_inv)))
print(KMO(R_inv))

acp_inv <- prcomp(mat_inv, center = TRUE, scale. = TRUE)
tv_a    <- tabla_varianza(acp_inv)
cat("\nVariance by component:\n"); print(tv_a, row.names = FALSE)
cat("Kaiser criterion: retain", sum(acp_inv$sdev^2 > 1), "component(s)\n\n")
cat("Loadings:\n"); print(round(acp_inv$rotation, 4))
print(plot_scree(tv_a))

# --- 9B. PCA on all variables (meteo + total power) ---
vars_b  <- c("Pot_Media","Rad_Media","Temp_Media","Humedad_Relativa",
             "Precipitacion","Viento_Velocidad","Presion")
mat_b   <- df_analysis |> select(all_of(vars_b))

cat("Bartlett test (9B):\n")
R_b <- cor(mat_b, use = "pairwise.complete.obs")
print(cortest.bartlett(R_b, n = nrow(mat_b)))
print(KMO(R_b))

acp_b <- prcomp(mat_b, center = TRUE, scale. = TRUE)
tv_b  <- tabla_varianza(acp_b)
cat("\nVariance by component:\n"); print(tv_b, row.names = FALSE)
cat("Kaiser criterion: retain", sum(acp_b$sdev^2 > 1), "component(s)\n")
cat("Loadings:\n"); print(round(acp_b$rotation, 4))
print(plot_scree(tv_b))

# --- 9C. PCA on meteorological variables only ---
vars_c <- c("Rad_Media","Temp_Media","Humedad_Relativa","Precipitacion",
            "Viento_Velocidad","Presion")
mat_c  <- meteo_analysis |> select(all_of(vars_c))

R_c <- cor(mat_c, use = "pairwise.complete.obs")
cat("Bartlett test (9C):\n"); print(cortest.bartlett(R_c, n = nrow(mat_c)))
print(KMO(R_c))

acp_c <- prcomp(mat_c, center = TRUE, scale. = TRUE)
tv_c  <- tabla_varianza(acp_c)
cat("\nVariance by component:\n"); print(tv_c, row.names = FALSE)
cat("Kaiser criterion: retain", sum(acp_c$sdev^2 > 1), "component(s)\n")
cat("Loadings:\n"); print(round(acp_c$rotation, 4))
print(plot_scree(tv_c))

# =============================================================================
# STEP 10 — CLASSICAL MDS ON VARIABLES
# =============================================================================
# Distance derived from Pearson r: d(Xi,Xj) = sqrt(2*(1 - r_ij))

vars_mds     <- c("Pot_Media","Rad_Media","Temp_Media",
                  "Hum_Rel","Precip","Viento_Vel","Presion")
vars_mds_src <- c("Pot_Media","Rad_Media","Temp_Media",
                  "Humedad_Relativa","Precipitacion","Viento_Velocidad","Presion")
n_v <- length(vars_mds)

mat_cor_mds <- df_analysis |>
  select(all_of(vars_mds_src)) |>
  rename(Hum_Rel = Humedad_Relativa, Precip = Precipitacion,
         Viento_Vel = Viento_Velocidad) |>
  cor(use = "pairwise.complete.obs")

D_vars <- sqrt(2 * (1 - mat_cor_mds))

# Classical MDS: B = HAH
A <- (-1/2) * D_vars^2
H <- diag(n_v) - matrix(1/n_v, n_v, n_v)
B <- H %*% A %*% H

evd        <- eigen(B)
lambda_mds <- evd$values
autovec    <- evd$vectors

lambda_pos <- lambda_mds[lambda_mds > 0]
prop_acum  <- cumsum(lambda_pos) / sum(lambda_pos)
cat("Cumulative proportion (positive eigenvalues):\n")
for (i in seq_along(prop_acum)) cat("  k =", i, "->", round(prop_acum[i], 3), "\n")

k_mds    <- 2
result   <- autovec[, 1:k_mds] %*% diag(sqrt(lambda_mds[1:k_mds]))
rownames(result) <- vars_mds

gof_mds   <- sum(lambda_mds[1:k_mds]) / sum(lambda_pos)
mds_check <- cmdscale(D_vars, k = k_mds, eig = TRUE)
cat("GOF (k=2):", round(gof_mds, 3),
    "| cmdscale GOF[2]:", round(mds_check$GOF[2], 3), "\n\n")

par(mar = c(4, 4, 2, 1))
plot(result[, 1], result[, 2], pch = 19, col = "#4682B4",
     xlab = "x", ylab = "y",
     xlim = range(result[, 1]) + c(-.2, .2),
     ylim = range(result[, 2]) + c(-.1, .15))
text(result[, 1], result[, 2] + 0.05, labels = vars_mds, cex = 0.9, font = 2)
abline(h = 0, v = 0, lty = 2, col = "grey70")

# =============================================================================
# STEP 11 — CLUSTERING
# =============================================================================

K_MAX <- ceiling(sqrt(nrow(df_analysis)))
cat("Days available:", nrow(df_analysis), "| K_max:", K_MAX, "\n\n")

run_kmeans_sweep <- function(mat, k_max, label) {
  wss <- bet <- rep(NA, k_max)
  for (i in seq_len(k_max)) {
    km     <- kmeans(mat, i, nstart = 25, algorithm = "MacQueen")
    wss[i] <- km$tot.withinss
    bet[i] <- km$betweenss
  }
  p_wss <- ggplot(data.frame(k = seq_len(k_max), WSS = wss), aes(k, WSS)) +
    geom_line(color = COL_LINEA, linewidth = 1) + geom_point(color = COL_PUNTO, size = 2.5) +
    scale_x_continuous(breaks = seq(1, k_max, 2)) +
    labs(x = "k", y = "Total Within SS", title = label) + theme_minimal(base_size = 10)
  print(p_wss)
  list(wss = wss, bet = bet)
}

# --- 11A. All variables ---
vars_a <- c("Pot_Media","Rad_Media","Temp_Media","Humedad_Relativa",
            "Precipitacion","Viento_Velocidad","Presion")
mat_a  <- df_analysis |> select(all_of(vars_a)) |> scale()
run_kmeans_sweep(mat_a, K_MAX, "All variables")

set.seed(123)
k_a  <- 3
km_a <- kmeans(mat_a, k_a, nstart = 25, algorithm = "MacQueen")
df_analysis$Cluster_A <- as.factor(km_a$cluster)
df_analysis <- reordenar_clusters(df_analysis, "Cluster_A")

cat("Cluster A distribution (k=", k_a, "):\n"); print(table(df_analysis$Cluster_A))
cat("BSS/TSS:", round(100 * km_a$betweenss / km_a$totss, 1), "%\n\n")

# --- 11B. Power only ---
mat_b <- df_analysis |> select(Pot_Media) |> scale()
run_kmeans_sweep(mat_b, K_MAX, "Power only")

set.seed(123)
k_b  <- 3
km_b <- kmeans(mat_b, k_b, nstart = 25, algorithm = "MacQueen")
df_analysis$Cluster_B <- as.factor(km_b$cluster)
df_analysis <- reordenar_clusters(df_analysis, "Cluster_B")

cat("Cluster B distribution (k=", k_b, "):\n"); print(table(df_analysis$Cluster_B))
cat("BSS/TSS:", round(100 * km_b$betweenss / km_b$totss, 1), "%\n\n")
cat("Contingency table A vs B:\n"); print(table(df_analysis$Cluster_A, df_analysis$Cluster_B))

fwrite(df_analysis, file.path(PATH_OUTPUTS, "df_analysis_clusters.csv"), dateTimeAs = "write.csv")

# =============================================================================
# STEP 12 — MULTIPLE LINEAR REGRESSION (MLR)
# =============================================================================

vars_rlm <- c("Pot_Media","Rad_Media","Temp_Media","Humedad_Relativa",
              "Precipitacion","Viento_Velocidad","Presion","Mes","Es_No_Laborable")

df_rlm <- df_analysis |>
  select(Fecha_Dia, all_of(vars_rlm)) |>
  mutate(
    mes_num = as.integer(factor(Mes, levels = MESES_ORDEN)),
    Mes_sin = sin(2 * pi * mes_num / 12),
    Mes_cos = cos(2 * pi * mes_num / 12),
    Laborable = factor(Es_No_Laborable, levels = c(FALSE, TRUE), labels = c("Si","No"))
  ) |>
  select(-Mes, -mes_num, -Es_No_Laborable) |>
  filter(complete.cases(.))

set.seed(123)
idx_train <- sample(nrow(df_rlm), round(.80 * nrow(df_rlm)))
df_train  <- df_rlm[idx_train, ]
df_test   <- df_rlm[-idx_train, ]
D_tr      <- select(df_train, -Fecha_Dia)

cat("Train:", nrow(df_train), "| Test:", nrow(df_test), "\n\n")

# --- 12.3 Full model ---
Mod_FULL <- lm(Pot_Media ~ ., data = D_tr)
print(summary(Mod_FULL))
print(Anova(Mod_FULL, type = "II"))

# --- 12.5 Model selection ---
Mod_NULL    <- lm(Pot_Media ~ 1, data = D_tr)
Mod_SEL     <- stepAIC(Mod_NULL, direction = "both", trace = FALSE,
                       scope = list(lower = Mod_NULL, upper = Mod_FULL))
cat("\nSelected model formula:", deparse(formula(Mod_SEL)), "\n")
print(summary(Mod_SEL))
cat("VIF:\n"); print(round(car::vif(Mod_SEL), 3))

# --- 12.6 Diagnostics ---
par(mfrow = c(2, 2)); plot(Mod_SEL); par(mfrow = c(1, 1))
car::qqPlot(residuals(Mod_SEL), col = "#4575B4", col.lines = "#CD5C5C",
            pch = 16, cex = 0.7, id = list(n = 3))
print(nortest::lillie.test(residuals(Mod_SEL)))
print(lmtest::bptest(Mod_SEL))
print(tseries::runs.test(factor(sign(residuals(Mod_SEL)))))

# Cook's distance check
cook_d   <- cooks.distance(Mod_SEL)
idx_cook <- which(cook_d > 4 / nrow(D_tr))
cat("Cook influential observations (>4/n):", length(idx_cook), "\n")

# --- 12.7 Cross-validation ---
set.seed(123)
ctrl_cv  <- trainControl(method = "repeatedcv", number = 5, repeats = 10)
cv_SEL   <- train(formula(Mod_SEL), data = D_tr, method = "lm", trControl = ctrl_cv)
print(cv_SEL$results)

pred_test <- predict(Mod_SEL, newdata = df_test)
obs_test  <- df_test$Pot_Media
RMSE_test <- sqrt(mean((obs_test - pred_test)^2))
R2_test   <- cor(obs_test, pred_test)^2
MAE_test  <- mean(abs(obs_test - pred_test))
cat(sprintf("MLR test — RMSE: %.3f | R2: %.3f | MAE: %.3f\n\n",
            RMSE_test, R2_test, MAE_test))

# =============================================================================
# STEP 13 — RANDOM FOREST
# =============================================================================

grid_rf <- expand.grid(mtry = 2:7, splitrule = "variance", min.node.size = 5)
set.seed(123)
rf_tuning <- train(formula(Mod_SEL), data = D_tr, method = "ranger",
                   trControl = ctrl_cv, tuneGrid = grid_rf,
                   num.trees = 500, importance = "impurity")

cat("Optimal mtry:", rf_tuning$bestTune$mtry, "\n")

set.seed(123)
rf_final <- ranger(formula   = formula(Mod_SEL), data = D_tr,
                   num.trees = 500, mtry = rf_tuning$bestTune$mtry,
                   importance = "impurity")

cat("OOB R²:", round(rf_final$r.squared, 3),
    "| OOB RMSE:", round(sqrt(rf_final$prediction.error), 3), "\n")

# Variable importance
importancia <- data.frame(Variable   = names(rf_final$variable.importance),
                          Importance = rf_final$variable.importance) |>
  arrange(desc(Importance)) |>
  mutate(Variable = factor(Variable, levels = rev(Variable)))

print(
  ggplot(importancia, aes(x = Variable, y = Importance)) +
    geom_col(fill = "#4682B4", alpha = 0.8) + coord_flip() +
    labs(x = NULL, y = "Importance (node impurity)") + theme_minimal(base_size = 10)
)

pred_rf_test  <- predict(rf_final, data = df_test)$predictions
RMSE_rf       <- sqrt(mean((obs_test - pred_rf_test)^2))
R2_rf         <- cor(obs_test, pred_rf_test)^2
MAE_rf        <- mean(abs(obs_test - pred_rf_test))
cat(sprintf("RF test — RMSE: %.3f | R2: %.3f | MAE: %.3f\n\n",
            RMSE_rf, R2_rf, MAE_rf))

saveRDS(rf_final, file.path(PATH_OUTPUTS, "rf_final.rds"))

# =============================================================================
# STEP 14 — FUNCTIONAL DATA ANALYSIS (FDA)
# =============================================================================
# The FDA section is structurally identical to Sections 14.1–14.4 of the
# original analysis. Key parameters are defined below for easy adjustment.

col_mes_pal <- colorRampPalette(c("#4682B4","#92C5DE","#CD5C5C","#B2182B"))(12)
select      <- dplyr::select

# Grid definitions
GRID_PROD <- seq(5, 20.5, by = 1/6); GRID_PROD <- GRID_PROD[GRID_PROD < 20.5]
GRID_24H  <- seq(0, 24,   by = 1/6); GRID_24H  <- GRID_24H[GRID_24H  < 24]
dias_ref  <- sort(as.Date(df_daily_clean$Fecha_Dia))

# --- 14.1 Build functional objects ---
construir_prod <- function(df10, dias, col, grid) {
  gr <- round(grid, 4)
  df_obs <- df10 |>
    mutate(Fecha_Dia = as.Date(Fecha, tz = "UTC"),
           Hora_Dec  = round(hour(Fecha) + minute(Fecha) / 60, 4)) |>
    filter(Fecha_Dia %in% dias, between(Hora_Dec, min(gr), max(gr))) |>
    select(Fecha_Dia, Hora_Dec, valor = all_of(col))

  expand_grid(Fecha_Dia = dias, Hora_Dec = gr) |>
    left_join(df_obs, by = c("Fecha_Dia","Hora_Dec")) |>
    mutate(valor = replace_na(valor, 0)) |>
    pivot_wider(names_from = Hora_Dec, values_from = valor) |>
    arrange(Fecha_Dia) -> wide

  list(matrix = as.matrix(wide[, -1]), fechas = wide$Fecha_Dia,
       valid  = rep(TRUE, nrow(wide)))
}

construir_meteo <- function(df_met, dias, col, grid, na_thr = 0.05) {
  gr <- round(grid, 4)
  df_obs <- df_met |>
    mutate(Fecha_Dia = as.Date(Fecha, tz = "UTC"),
           Hora_Dec  = round(hour(Fecha) + minute(Fecha) / 60, 4)) |>
    filter(Fecha_Dia %in% dias, Hora_Dec %in% gr) |>
    select(Fecha_Dia, Hora_Dec, valor = all_of(col))

  expand_grid(Fecha_Dia = dias, Hora_Dec = gr) |>
    left_join(df_obs, by = c("Fecha_Dia","Hora_Dec")) |>
    pivot_wider(names_from = Hora_Dec, values_from = valor) |>
    arrange(Fecha_Dia) -> wide

  M   <- as.matrix(wide[, -1])
  pct <- rowSums(is.na(M)) / ncol(M)

  for (i in seq_len(nrow(M))) {
    if (pct[i] > 0 && pct[i] <= na_thr) {
      x <- M[i, ]
      if (sum(!is.na(x)) >= 2)
        M[i, ] <- approx(which(!is.na(x)), x[!is.na(x)],
                         xout = seq_along(x), rule = 2)$y
    }
  }
  list(matrix = M, fechas = wide$Fecha_Dia, valid = pct <= na_thr)
}

cat("Building functional matrices...\n")
r_pot    <- construir_prod(df_final, dias_ref, "Potencia_AC_Total", GRID_PROD)
r_rad    <- construir_prod(df_final, dias_ref, "Radiacion_Solar",   GRID_PROD)
r_temp   <- construir_meteo(df_meteo, dias_ref, "Temperatura_1.5m", GRID_24H)
r_hum    <- construir_meteo(df_meteo, dias_ref, "Humedad_Relativa",  GRID_24H)
r_viento <- construir_meteo(df_meteo, dias_ref, "Viento_Velocidad",  GRID_24H)
r_pres   <- construir_meteo(df_meteo, dias_ref, "Presion",           GRID_24H)

validos_all <- r_pot$valid & r_rad$valid & r_temp$valid &
               r_hum$valid & r_viento$valid & r_pres$valid
cat("Valid days across all variables:", sum(validos_all), "of", length(validos_all), "\n")

df_meta <- df_daily_clean |>
  mutate(Fecha_Dia = as.Date(Fecha_Dia)) |>
  arrange(Fecha_Dia) |>
  select(Fecha_Dia, Mes, Dia_Semana, Es_No_Laborable, Es_Festivo,
         Pot_Media, Rad_Media, Temp_Media, Humedad_Relativa,
         Precipitacion, Viento_Velocidad, Presion) |>
  mutate(
    Mes             = factor(Mes, levels = MESES_ORDEN),
    Es_No_Laborable = as.logical(Es_No_Laborable),
    Es_Festivo      = as.logical(Es_Festivo),
    Llueve          = factor(Precipitacion > 0.5, labels = c("Dry","Rain")),
    Cielo_Despejado = factor(Rad_Media > 250,     labels = c("Overcast","Clear"))
  )

idx    <- which(validos_all)
fechas <- dias_ref[validos_all]
df_f   <- as.data.frame(df_meta[validos_all, ])
rownames(df_f) <- format(fechas)

hacer_fdata <- function(res, idx, grid, rango, label, unit) {
  M <- res$matrix[idx, , drop = FALSE]; rownames(M) <- format(fechas)
  fdata(M, round(grid, 4), rango,
        names = list(main = paste0("Daily curves — ", label),
                     xlab = "Hour (UTC)", ylab = paste0(label, " (", unit, ")")))
}

ldata_fv <- ldata(
  "df"       = df_f,
  "potencia" = hacer_fdata(r_pot,    idx, GRID_PROD, c(5,20.5), "AC Power",     "kW"),
  "rad"      = hacer_fdata(r_rad,    idx, GRID_PROD, c(5,20.5), "Radiation",    "W/m²"),
  "temp"     = hacer_fdata(r_temp,   idx, GRID_24H,  c(0,24),   "Temperature",  "°C"),
  "hum"      = hacer_fdata(r_hum,    idx, GRID_24H,  c(0,24),   "Humidity",     "%"),
  "viento"   = hacer_fdata(r_viento, idx, GRID_24H,  c(0,24),   "Wind speed",   "km/h"),
  "presion"  = hacer_fdata(r_pres,   idx, GRID_24H,  c(0,24),   "Pressure",     "hPa")
)

saveRDS(ldata_fv, file.path(PATH_OUTPUTS, "ldata_fv.rds"))

# --- 14.2 Smoothing ---
# Power & radiation: smooth.pos (guarantees positivity via exp(W(t))).
# Meteorological variables: B-spline with 1-SE rule (Breiman criterion).

suavizar_pos <- function(fdobj, nbasis = 15, lambda = 0.1) {
  t   <- fdobj$argvals; rng <- fdobj$rangeval; n <- nrow(fdobj$data)
  bsp <- create.bspline.basis(rng, nbasis = nbasis)
  wp  <- fdPar(fd(matrix(0, nbasis, 1), bsp), Lfdobj = 2, lambda = lambda)
  M   <- matrix(0, n, length(t))
  for (i in seq_len(n)) {
    yi    <- pmax(fdobj$data[i, ], 0.01)
    res_i <- tryCatch(smooth.pos(t, yi, wp), error = function(e) NULL)
    M[i, ] <- if (!is.null(res_i)) as.vector(exp(eval.fd(t, res_i$Wfdobj)))
              else pmax(fdobj$data[i, ], 0)
  }
  rownames(M) <- rownames(fdobj$data)
  fdata(M, t, rng, names = fdobj$names)
}

suavizar_1se <- function(fd, k_range, seed = 42) {
  set.seed(seed)
  res_full  <- optim.basis(fd, numbasis = k_range, lambda = 0)
  gcv_vec   <- res_full$gcv
  K_min     <- k_range[which.min(gcv_vec)]
  gcv_min   <- min(gcv_vec)
  resid_mat <- fd$data - res_full$fdata.est$data
  se        <- sd(rowMeans(resid_mat^2) / (1 - K_min/ncol(fd$data))^2) / sqrt(nrow(fd$data))
  K_1se     <- k_range[min(which(gcv_vec <= gcv_min + se))]
  res_1se   <- optim.basis(fd, numbasis = K_1se, lambda = 0)
  cat(sprintf("  K_min=%2d (GCV=%.4f) | K_1se=%2d\n", K_min, gcv_min, K_1se))
  list(fdata_1se = res_1se$fdata.est, K_min = K_min, K_1se = K_1se,
       gcv = gcv_vec, gcv_min = gcv_min, se = se,
       umbral = gcv_min + se, rango = k_range)
}

cat("\nSmoothing power and radiation (smooth.pos)...\n")
pot_pos <- suavizar_pos(ldata_fv$potencia)
rad_pos <- suavizar_pos(ldata_fv$rad)

cat("\nB-spline + 1-SE rule for meteorological variables:\n")
cat("· temperature: ");  res_s_temp   <- suavizar_1se(ldata_fv$temp,    11:100)
cat("· humidity:    ");  res_s_hum    <- suavizar_1se(ldata_fv$hum,     11:100)
cat("· wind:        ");  res_s_viento <- suavizar_1se(ldata_fv$viento,  11:50)
cat("· pressure:    ");  res_s_pres   <- suavizar_1se(ldata_fv$presion, 11:100)

ldata_fv$potencia <- pot_pos
ldata_fv$rad      <- rad_pos
ldata_fv$temp     <- res_s_temp$fdata_1se
ldata_fv$hum      <- res_s_hum$fdata_1se
ldata_fv$viento   <- res_s_viento$fdata_1se
ldata_fv$presion  <- res_s_pres$fdata_1se

saveRDS(ldata_fv, file.path(PATH_OUTPUTS, "ldata_fv_smooth.rds"))

# --- 14.3 Functional depth and outlier detection ---
set.seed(123)
ldata_raw  <- readRDS(file.path(PATH_OUTPUTS, "ldata_fv.rds"))
fdat_raw   <- ldata_raw$potencia
col_mes    <- col_mes_pal[as.integer(format(ldata_raw$df$Fecha_Dia, "%m"))]

depth_fm   <- depth.FM(fdat_raw, trim = 0.25)
cat("FM depth summary:\n"); print(summary(depth_fm$dep))

# Seasonal stratification
mes_num   <- as.integer(format(ldata_raw$df$Fecha_Dia, "%m"))
estacion  <- case_when(
  mes_num %in% c(12,1,2)  ~ "Winter",
  mes_num %in% c(3,4,5)   ~ "Spring",
  mes_num %in% c(6,7,8)   ~ "Summer",
  mes_num %in% c(9,10,11) ~ "Autumn"
)

for (s in c("Winter","Spring","Summer","Autumn")) {
  idx_s <- which(estacion == s); fd_s <- fdat_raw[idx_s]
  pond  <- outliers.depth.pond(fd_s, nb=200, smo=0.05, trim=0.01, dfunc=depth.FM)
  trim  <- outliers.depth.trim(fd_s, nb=200, smo=0.05, trim=0.01, dfunc=depth.FM)
  n_com <- length(intersect(pond$outliers, trim$outliers))
  cat(sprintf("%s (n=%d): pond=%d | trim=%d | confirmed=%d\n",
              s, length(idx_s), length(pond$outliers), length(trim$outliers), n_com))
}

# Confirmed outliers (manually verified after seasonal analysis)
FECHAS_OUTLIERS_FDA <- as.Date(c("2026-02-23","2026-02-26","2025-07-27"))

ldata_smo <- readRDS(file.path(PATH_OUTPUTS, "ldata_fv_smooth.rds"))
idx_out   <- match(as.character(FECHAS_OUTLIERS_FDA),
                   as.character(ldata_smo$df$Fecha_Dia))
idx_ok    <- setdiff(seq_len(nrow(ldata_smo$df)), idx_out)

ldata_clean <- list(
  potencia = ldata_smo$potencia[idx_ok],
  rad      = ldata_smo$rad[idx_ok],
  temp     = ldata_smo$temp[idx_ok],
  hum      = ldata_smo$hum[idx_ok],
  viento   = ldata_smo$viento[idx_ok],
  presion  = ldata_smo$presion[idx_ok],
  df       = ldata_smo$df[idx_ok, ]
)

cat("Curves for FPCA:", nrow(ldata_clean$df),
    "(removed", length(idx_out), "outliers)\n\n")

# --- 14.4 Functional PCA ---
K_POT <- 3; K_TEMP <- 2; K_RAD <- 3
K_HUM <- 5; K_VIENTO <- 5; K_PRES <- 1

pc_pot    <- create.pc.basis(ldata_clean$potencia, 1:8)
pc_temp   <- create.pc.basis(ldata_clean$temp,     1:6)
pc_rad    <- create.pc.basis(ldata_clean$rad,      1:8)
pc_hum    <- create.pc.basis(ldata_clean$hum,      1:6)
pc_viento <- create.pc.basis(ldata_clean$viento,   1:6)
pc_presion<- create.pc.basis(ldata_clean$presion,  1:6)

cat("FPCA variance summaries:\n")
cat("Power:      "); summary(pc_pot)
cat("Radiation:  "); summary(pc_rad)
cat("Temperature:"); summary(pc_temp)

# Sign convention: flip PC1 of temperature so positive scores = warm days
pc_temp$basis$data[1, ] <- -pc_temp$basis$data[1, ]

# --- 14.5 Functional regression ---
fechas_comunes  <- intersect(df_rlm$Fecha_Dia, ldata_clean$df$Fecha_Dia)
idx_ldata       <- which(ldata_clean$df$Fecha_Dia %in% fechas_comunes)
idx_rlm_comun   <- which(df_rlm$Fecha_Dia %in% fechas_comunes)

set.seed(123)
n_comun         <- length(fechas_comunes)
idx_tr          <- sample(n_comun, round(.80 * n_comun))
idx_te          <- setdiff(seq_len(n_comun), idx_tr)

subset_ldata <- function(ld, idx) {
  lapply(names(ld), function(nm) {
    if (nm == "df") ld$df[idx, ] else ld[[nm]][idx, ]
  }) |> setNames(names(ld))
}
ldata_train <- subset_ldata(ldata_clean, idx_ldata[idx_tr])
ldata_test  <- subset_ldata(ldata_clean, idx_ldata[idx_te])
y_train     <- ldata_train$df$Pot_Media
y_test      <- ldata_test$df$Pot_Media

LAM_GRID <- 10^seq(-2, 4, length.out = 10)

# Recalculate PC bases on training set
pc_rad_tr     <- create.pc.basis(ldata_train$rad,     1:K_RAD)
pc_temp_tr    <- create.pc.basis(ldata_train$temp,    1:K_TEMP)
pc_hum_tr     <- create.pc.basis(ldata_train$hum,     1:K_HUM)
pc_viento_tr  <- create.pc.basis(ldata_train$viento,  1:K_VIENTO)
pc_presion_tr <- create.pc.basis(ldata_train$presion, 1:K_PRES)

eval_model <- function(mod, y_obs, fd_test, tag) {
  pred  <- predict(mod, fd_test)
  rmse  <- sqrt(mean((y_obs - pred)^2))
  r2    <- 1 - sum((y_obs - pred)^2) / sum((y_obs - mean(y_obs))^2)
  cat(sprintf("%-6s test — R2: %.4f | RMSE: %.4f kW\n", tag, r2, rmse))
  invisible(list(pred=pred, rmse=rmse, r2=r2))
}

# M1a: radiation
m1a <- fregre.pc(ldata_train$rad, y_train, l = 1:K_RAD)
eval_model(m1a, y_test, ldata_test$rad, "M1a")

# M1b: temperature
m1b <- fregre.pc(ldata_train$temp, y_train, l = 1:K_TEMP)
eval_model(m1b, y_test, ldata_test$temp, "M1b")

# M2a: radiation + CV
m2a_cv <- fregre.pc.cv(ldata_train$rad, y_train, kmax=6, lambda=LAM_GRID, P=c(0,0,1))
m2a    <- m2a_cv$fregre.pc
eval_model(m2a, y_test, ldata_test$rad, "M2a")

# M2b: temperature + CV
m2b_cv <- fregre.pc.cv(ldata_train$temp, y_train, kmax=6, lambda=LAM_GRID, P=c(0,0,1))
m2b    <- m2b_cv$fregre.pc
eval_model(m2b, y_test, ldata_test$temp, "M2b")

# M3: rad + temp + hum
m3 <- fregre.lm(Pot_Media ~ rad + temp + hum, data = ldata_train,
                basis.x = list(rad=pc_rad_tr, temp=pc_temp_tr, hum=pc_hum_tr))
eval_model(m3, y_test, ldata_test, "M3")

# M4: rad + temp + hum + viento + presion + Rad_Media
m4 <- fregre.lm(
  Pot_Media ~ rad + temp + hum + viento + presion + Rad_Media,
  data    = ldata_train,
  basis.x = list(rad=pc_rad_tr, temp=pc_temp_tr, hum=pc_hum_tr,
                 viento=pc_viento_tr, presion=pc_presion_tr)
)
eval_model(m4, y_test, ldata_test, "M4")

# Scalar reference (Mod_SEL) on the same test subset
df_test_comun  <- df_rlm[idx_rlm_comun[idx_te], ]
pred_m0        <- predict(Mod_SEL, newdata = df_test_comun)
obs_comun      <- df_test_comun$Pot_Media
cat(sprintf("M0 (MLR) test  — R2: %.4f | RMSE: %.4f kW\n",
            cor(obs_comun, pred_m0)^2,
            sqrt(mean((obs_comun - pred_m0)^2))))

# RF on same subset
pred_rf_comun <- predict(rf_final, data = df_test_comun)$predictions
cat(sprintf("RF test (same) — R2: %.4f | RMSE: %.4f kW\n",
            cor(obs_comun, pred_rf_comun)^2,
            sqrt(mean((obs_comun - pred_rf_comun)^2))))

cat("\n=== Analysis complete ===\n")
cat("Outputs saved to:", PATH_OUTPUTS, "\n")
