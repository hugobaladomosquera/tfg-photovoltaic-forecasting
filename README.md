# Photovoltaic Generation Forecasting at Navantia Ferrol
### Multivariate and Functional Data Analysis on Industrial SCADA Data

**Bachelor's Thesis (TFG) — Honours distinction (Matrícula de Honor)**  
B.Sc. Data Science & Engineering · Universidade da Coruña · 2026  
Author: Hugo Balado Mosquera

---

## Overview

This project develops a full predictive pipeline for photovoltaic energy generation at Navantia Ferrol's module workshop, using 432 days of asynchronous SCADA records from three inverters combined with meteorological data from the MeteoGalicia CIS Ferrol station.

The analysis is split into two complementary approaches:

- **Multivariate analysis** — PCA, MDS, k-means clustering, multiple linear regression and Random Forest on daily aggregated data.
- **Functional data analysis (FDA)** — daily power curves treated as functional objects, smoothed via B-splines and `smooth.pos`, with functional PCA and functional regression.

**Results on test set (85 days):**

| Model | R² | RMSE (kW) |
|---|---|---|
| M0 — MLR scalar reference | 0.812 | 15.30 |
| M1a — radiation curve (functional) | 0.881 | 12.17 |
| M3 — radiation + temperature + humidity | 0.902 | 11.08 |
| M4 — M3 + wind + pressure | 0.906 | 10.84 |
| RF — Random Forest | **0.911** | **10.54** |

M3 is selected as the final functional model (best parsimony/performance trade-off among functional models).

---

## Repository structure

```
tfg-photovoltaic-forecasting/
├── README.md
├── src/
│   └── analisis_fotovoltaico.R   ← main analysis script
├── docs/
│   └── TFG_GCED_Hugo_Balado.pdf  ← full thesis (Spanish)
└── data/
    └── README.md                 ← data description (files not included)
```

---

## Data

The raw SCADA and meteorological files are **not included** in this repository as they contain proprietary Navantia operational data.

The script expects the following files under `data/`:

| File | Description |
|---|---|
| `FV_Rep1.csv` | Inverter 1 — AC power time series (SCADA sparse format) |
| `FV_Rep2.csv` | Inverter 2 — AC power time series |
| `FV_Rep3.csv` | Inverter 3 — AC power time series |
| `METEO_MAESTRA_TFG_2024_2026.csv` | Meteorological data — MeteoGalicia CIS Ferrol station |

To run the script with your own data, update the path constants at the top of `src/analisis_fotovoltaico.R`:

```r
DATA_DIR <- "data/"
```

---

## Methods

### Multivariate pipeline
- **Preprocessing** — SCADA sparse format parsing, 10-minute regularisation, UTC timezone correction for MeteoGalicia data
- **Anomaly detection** — technical shutdown detection (solar days with near-zero power), univariate outlier detection (Tukey IQR), multivariate outlier detection (LOF)
- **Exploratory analysis** — Pearson and Székely distance correlations, univariate/bivariate EDA, monthly boxplots
- **Dimensionality reduction** — PCA (3 variants: inverters, all variables, meteorological only), classical MDS on variables
- **Clustering** — k-means (k=3) on all variables, meteorological variables only, and power only
- **Regression** — Multiple linear regression with stepwise AIC selection, diagnostic tests (Lilliefors, Breusch-Pagan, runs test), Cook's distance analysis
- **Random Forest** — `ranger` with 5-fold repeated CV tuning of `mtry`

### Functional Data Analysis
- **Smoothing** — `smooth.pos` for power and radiation (positivity guaranteed via exp(W(t))); B-spline with 1-SE rule (Breiman criterion) for meteorological variables
- **Functional depth** — Fraiman-Muniz depth with seasonal stratification for outlier detection
- **FPCA** — functional principal component analysis on all six variables; components selected by ≥90% cumulative variance threshold
- **Functional regression** — scalar-on-function models (M1–M4) via `fregre.pc`, `fregre.pc.cv` and `fregre.lm`

---

## Tech stack

- **R** — `fda.usc`, `tidyverse`, `data.table`, `lubridate`, `ranger`, `caret`, `car`, `MASS`, `dbscan`, `energy`, `corrplot`, `lmtest`, `nortest`, `tseries`
- **Python** — used in parallel SCADA pipeline work at Navantia (not included here)

---

## How to run

1. Clone the repository:
```bash
git clone https://github.com/hugobaladomosquera/tfg-photovoltaic-forecasting.git
cd tfg-photovoltaic-forecasting
```

2. Install R dependencies:
```r
install.packages(c(
  "fda.usc", "tidyverse", "lubridate", "data.table", "energy",
  "dbscan", "psych", "cluster", "ranger", "caret", "corrplot",
  "MASS", "car", "leaps", "mctest", "ggm", "Hmisc",
  "lmtest", "tseries", "nortest"
))
```

3. Place your data files in `data/` and run:
```r
source("src/analisis_fotovoltaico.R")
```

Outputs (plots, CSVs, model `.rds` files) are saved to `outputs/` automatically.

---

## License

Academic project. Code is open for reference and reuse. Raw data is proprietary to Navantia S.A. and not redistributable.

---

*For questions or collaboration: [hugobalado2004@gmail.com](mailto:hugobalado2004@gmail.com) · [LinkedIn](https://linkedin.com/in/hugo-balado-mosquera)*
