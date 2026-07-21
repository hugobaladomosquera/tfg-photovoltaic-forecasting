# Data

The raw data files are not included in this repository as they contain
proprietary operational data from Navantia S.A. (Ferrol, Spain).

## Expected files

| File | Description |
|---|---|
| `FV_Rep1.csv` | Inverter 1 — AC power (SCADA sparse format) |
| `FV_Rep2.csv` | Inverter 2 — AC power |
| `FV_Rep3.csv` | Inverter 3 — AC power |
| `METEO_MAESTRA_TFG_2024_2026.csv` | Meteorological data — MeteoGalicia CIS Ferrol |

## Format

SCADA files use a 3-line header, `;` separator and
`YYYY.MM.DD HH:MM:SS.sss` timestamp format.
Meteorological missing values are encoded as `-9999`.
