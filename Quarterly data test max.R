# =============================================================================
# Script: US Inflation Data Preparation and Analysis (Quarterly Version)
# Description: Only uses CPInonSA to compute quarterly inflation.
# Uses only end-of-quarter months and exogenous variables from FREDQ.csv.
# =============================================================================

# Load required libraries
library(tidyverse)
library(readxl)
library(lubridate)

# Load non-seasonally adjusted CPI data from Excel
CPI_US_labour_dataset <- read_excel("CPI US labour dataset.xlsx", 
                                    range = "A12:M124")

# Reshape from wide to long format and create date column
CPI_US_labour_long <- CPI_US_labour_dataset %>%
  pivot_longer(
    cols = -1,
    names_to = "Month",
    values_to = "CPInonSA"
  ) %>%
  rename(Year = 1) %>%
  filter(Year >= 1959) %>%
  mutate(
    Month_num = match(Month, month.abb),
    Date = as.Date(paste(Year, Month_num, "01", sep = "-"))
  ) %>%
  arrange(Date) %>%
  filter(month(Date) %in% c(3, 6, 9, 12))  # Keep only end-of-quarter months

# Load FREDQ quarterly data (seasonally adjusted exogenous vars)
fredq <- read.csv("FREDQ.csv")
fredq <- fredq %>%
  mutate(sasdate = as.Date(sasdate, format = "%m/%d/%Y")) %>%
  filter(sasdate >= as.Date("1959-03-01"), sasdate < as.Date("2025-01-01"))

# Merge with non-seasonally adjusted CPI
inflation_df_quarterly <- fredq %>%
  left_join(
    CPI_US_labour_long %>%
      select(Date, CPInonSA),
    by = c("sasdate" = "Date")
  ) %>%
  arrange(sasdate) %>%
  mutate(
    Quarter = quarter(sasdate),
    Year = year(sasdate),
    CPInonSA = as.numeric(CPInonSA)
  ) %>%
  # Calculate inflation: log difference of CPInonSA
  mutate(
    inflationNonSA = 100 * (log(CPInonSA) - log(lag(CPInonSA)))
  ) %>%
  drop_na(inflationNonSA)

# Preview or save
# save(inflation_df_quarterly, file = "inflation_df_quarterly.RData")

# =============================================================================
# Script: MAR Model Estimation and Forecast Evaluation
# Description:
#   This script estimates and compares Mixed Causal-Noncausal Autoregressive (MAR)
#   models for US inflation using non-seasonally adjusted data. It includes:
#     - Model selection via information criteria
#     - Residual diagnostics for AR models
#     - Forecasting with various model configurations
#     - Out-of-sample forecast performance evaluation using RMSE
# =============================================================================

# Load required libraries
library(MASS)          # For statistical distributions and matrix functions
source("MARX_functions.R")  # Custom functions for MAR model estimation
source("MART.R")            # MART model training and forecasting routines (includes information criteria calculations)
library(forecast)      # For ARIMA modeling and forecast tools
library(pbmcapply)     # For parallel processing with progress bar
library(stats)

#basic MAR
marx(inflation_df_quarterly$inflationNonSA, NULL, p_max = 6, sig_level = 0.1)

# Define maximum lag orders for causal (p_C) and noncausal (p_NC) components
p_C_max <- 6
p_NC_max <- 6

# Initialize matrices to store information criteria and model diagnostics
LL <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)
AIC <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)
BIC <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)
HQ <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)
N <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)

# Loop over all possible (p_C, p_NC) combinations and compute info criteria
for (i in 0:p_C_max) {
  for (j in 0:p_NC_max) {
    marx_loop <- marx.t(inflation_df_quarterly$inflationNonSA, NULL, p_C = i, p_NC = j)
    information <- information.criteria(type = "MARX", marx_loop)
    LL[(i+1),(j+1)] <- information$loglikelihood
    AIC[(i+1),(j+1)] <- information$aic
    BIC[(i+1),(j+1)] <- information$bic
    HQ[(i+1),(j+1)] <- information$hq
    N[(i+1),(j+1)] <- information$n
  }
}


# -----------------------------------------------------------------------------
# Multi-horizon out-of-sample forecast evaluation
# -----------------------------------------------------------------------------

# Forecasting parameters
h <- 4         # Forecast horizon
N <- 15000        # Posterior draws
M <- 50          # MA truncation

# Model specifications
p_C_mixed <- 2;  p_NC_mixed <- 2    # Mixed MAR(2,2)
p_C_causal <- 4; p_NC_causal <- 0   # Purely causal AR(4)

# Define forecast evaluation window
data_series <- inflation_df_quarterly$inflationNonSA
start_index <- 100
end_index <- length(data_series) - h
forecast_indices <- start_index:end_index

# -----------------------------------------------------------------------------
# Forecasting loop using all 3 model configurations
# -----------------------------------------------------------------------------

results_list <- pbmclapply(
  X = forecast_indices,
  FUN = function(t) {
    y_window <- data_series[1:t]  # Expanding window up to time t
    
    forecast_mixed <- forecast.marx(
      y = y_window,
      p_C = p_C_mixed,
      p_NC = p_NC_mixed,
      h = h,
      M = M,
      N = N
    )
    
    forecast_causal <- forecast.marx(
      y = y_window,
      p_C = p_C_causal,
      p_NC = p_NC_causal,
      h = h,
      M = M,
      N = N
    )
    
    actual <- data_series[(t + 1):(t + h)]  # Future observed values
    
    return(list(
      mixed = forecast_mixed,
      causal = forecast_causal,
      actual = actual
    ))
  },
  mc.cores = parallel::detectCores() - 1  # Use all but one core
)

# -----------------------------------------------------------------------------
# Organize forecast results into matrices
# -----------------------------------------------------------------------------

forecast_mixed <- do.call(rbind, lapply(results_list, `[[`, "mixed"))
forecast_causal <- do.call(rbind, lapply(results_list, `[[`, "causal"))
forecast_mid <- do.call(rbind, lapply(results_list, `[[`, "mid"))
actual_matrix <- do.call(rbind, lapply(results_list, `[[`, "actual"))

colnames(forecast_mixed) <- paste0("h", 1:h)
colnames(forecast_causal) <- paste0("h", 1:h)
colnames(forecast_mid) <- paste0("h", 1:h)
colnames(actual_matrix) <- paste0("h", 1:h)
