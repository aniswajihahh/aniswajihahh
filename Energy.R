# ===============================
# Load required libraries
# ===============================
library(tsibble)
library(dplyr)
library(fable)
library(ggplot2)
library(arrow)
library(tseries)
library(vars)
library(urca)
library(lmtest)

# ===============================
# 1. Load electricity consumption dataset
# ===============================
EC <- read_parquet("https://storage.data.gov.my/energy/electricity_consumption.parquet")

# ===============================
# 2. Convert dataset into a tsibble object
# - date as time index
# - sector as key
# ===============================
EC_t <- EC %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d")) %>% 
  as_tsibble(key = sector, index = date)

# ======================================================
# 3. Time Series Plots of Electricity Consumption by Sector
# ======================================================
# -------------------------------
# a) Commercial Sector
# -------------------------------
  # Filter data for the local commercial sector
  local_commercial <- EC_t %>%
    filter(sector == "local_commercial")
  
  # Plot monthly electricity consumption
  ggplot(local_commercial, aes(x = date, y = consumption)) +
    geom_line(color = "steelblue", linewidth = 1) +
    # Dashed line indicates the sample mean for visual reference
    geom_hline(
      yintercept = mean(local_commercial$consumption, na.rm = TRUE),
      color = "red",
      linetype = "dashed"
    ) +
    labs(
      title = "Monthly Electricity Consumption: Commercial Sector",
      x = "Month",
      y = "Electricity Consumption (kWh)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

# -------------------------------
# b) Domestic Sector
# -------------------------------
# Filter data for the local domestic sector
local_domestic <- EC_t %>%
  filter(sector == "local_domestic")

# Plot monthly electricity consumption
ggplot(local_domestic, aes(x = date, y = consumption)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  # Dashed line indicates the sample mean for visual reference
  geom_hline(
    yintercept = mean(local_domestic$consumption, na.rm = TRUE),
    color = "red",
    linetype = "dashed"
  ) +
  labs(
    title = "Monthly Electricity Consumption: Domestic Sector",
    x = "Month",
    y = "Electricity Consumption (kWh)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# -------------------------------
# c) Export Sector
# -------------------------------
# Filter data for the export sector
exports <- EC_t %>%
  filter(sector == "exports")

# Plot monthly electricity consumption
ggplot(exports, aes(x = date, y = consumption)) +
  geom_line(color = "purple", linewidth = 1) +
  # Dashed line indicates the sample mean for visual reference
  geom_hline(
    yintercept = mean(exports$consumption, na.rm = TRUE),
    color = "red",
    linetype = "dashed"
  ) +
  labs(
    title = "Monthly Electricity Consumption: Export Sector",
    x = "Month",
    y = "Electricity Consumption (kWh)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# -------------------------------
# d) Energy Losses
# -------------------------------
# Filter data for electricity losses
losses <- EC_t %>%
  filter(sector == "losses")

# Plot monthly electricity losses
ggplot(losses, aes(x = date, y = consumption)) +
  geom_line(color = "orange", linewidth = 1) +
  # Dashed line indicates the sample mean for visual reference
  geom_hline(
    yintercept = mean(losses$consumption, na.rm = TRUE),
    color = "red",
    linetype = "dashed"
  ) +
  labs(
    title = "Monthly Electricity Losses",
    x = "Month",
    y = "Electricity Losses (kWh)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# ===============================
# 4. Extract consumption variables by sector
# ===============================
Commercial <- local_commercial$consumption
Domestic   <- local_domestic$consumption 
Export     <- exports$consumption       
Loss       <- losses$consumption

# ===============================
# 5. Augmented Dickey-Fuller (ADF) Test
# Purpose: To test stationarity of each time series
# ===============================

# Commercial sector
adf.test(Commercial)
dCommercial <- diff(Commercial)
adf.test(dCommercial)

# Domestic sector
adf.test(Domestic)

# Export sector
adf.test(Export)
dExport <- diff(Export)
adf.test(dExport)

# Energy losses
adf.test(Loss)
dLoss <- diff(Loss)
adf.test(dLoss)

# ===============================
# Prepare data for cointegration and VAR analysis
# ===============================

# Level data (non-stationary series)
data_level <- data.frame(
  Local_Commercial = Commercial,
  Exports = Export,
  Losses = Loss
)

# First-differenced (stationary) data
data_pegun <- data.frame(dCommercial, dExport, dLoss)
colnames(data_pegun) <- c("Local_Commercial", "Exports", "Losses")

# ===============================
# 6. Optimal Lag Length Selection for VAR
# ===============================
VARselect(data_pegun, lag.max = 12, type = "const")

# ===============================
# 7. Johansen Cointegration Test
# ===============================

# Trace test
johansen_test <- ca.jo(
  data_level,
  type = "trace",
  ecdet = "const",
  K = 1
)
summary(johansen_test)

# Maximum eigenvalue test
johansen_test2 <- ca.jo(
  data_level,
  type = "eigen",
  ecdet = "const",
  K = 1
)
summary(johansen_test2)

# ===============================
# 8. Vector Error Correction Model (VECM)
# ===============================
vecm_model1 <- cajorls(johansen_test, r = 1)
summary(vecm_model1$rlm)

# ===============================
# Residual Diagnostics for VECM
# ===============================
residuals_vecm <- resid(vecm_model1$rlm)
residuals_vecm <- as.matrix(residuals_vecm)
colnames(residuals_vecm)

par(mfrow = c(2,1))

# Commercial sector residuals
acf(residuals_vecm[,1],   main = "Residual ACF: Commercial")
pacf(residuals_vecm[,1], main = "Residual PACF: Commercial")

# Export sector residuals
acf(residuals_vecm[,2],   main = "Residual ACF: Export")
pacf(residuals_vecm[,2], main = "Residual PACF: Export")

# Energy losses residuals
acf(residuals_vecm[,3],   main = "Residual ACF: Energy Losses")
pacf(residuals_vecm[,3], main = "Residual PACF: Energy Losses")

# ===============================
# 9. Granger Causality Analysis
# Purpose: To examine causal relationships between sectors
# ===============================
grangertest(Local_Commercial ~ Exports, order = 1, data = data_level)
grangertest(Local_Commercial ~ Losses, order = 1, data = data_level)

grangertest(Exports ~ Local_Commercial, order = 1, data = data_level)
grangertest(Exports ~ Losses, order = 1, data = data_level)

grangertest(Losses ~ Local_Commercial, order = 1, data = data_level)
grangertest(Losses ~ Exports, order = 1, data = data_level)

