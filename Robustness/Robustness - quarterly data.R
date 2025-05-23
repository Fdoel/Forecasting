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
#marx(inflation_df_quarterly$inflationNonSA, NULL, p_max = 6, sig_level = 0.1)

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
data <- as.matrix(inflation_df_quarterly$inflationNonSA)
mar_model_quarterly <- marx.t(data, NULL, p_C = p_C_mixed, p_NC = p_NC_mixed)
p_C_causal <- 4; p_NC_causal <- 0   # Purely causal AR(4)
ar_model_quarterly <- marx.t(data, NULL, p_C = p_C_causal, p_NC = p_NC_causal)

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
actual_matrix <- do.call(rbind, lapply(results_list, `[[`, "actual"))

colnames(forecast_mixed) <- paste0("h", 1:h)
colnames(forecast_causal) <- paste0("h", 1:h)
colnames(actual_matrix) <- paste0("h", 1:h)

# -----------------------------------------------------------------------------
# Compute RMSE for each model across horizons
# -----------------------------------------------------------------------------

rmse <- function(forecast, actual) {
  sqrt(colMeans((forecast - actual)^2, na.rm = TRUE))
}

rmse_mixed <- rmse(forecast_mixed, actual_matrix)
rmse_causal <- rmse(forecast_causal, actual_matrix)

# -----------------------------------------------------------------------------
# Compute Diebold-Mariano test p-values
# -----------------------------------------------------------------------------

compute_dm_tests <- function(forecast1, forecast2, actual, h) {
  p_values <- numeric(h)
  for (i in 1:h) {
    e1 <- actual[, i] - forecast1[, i]
    e2 <- actual[, i] - forecast2[, i]
    valid <- complete.cases(e1, e2)
    e1 <- e1[valid]
    e2 <- e2[valid]
    
    if (length(e1) > 10) {
      dm <- tryCatch(
        dm.test(e1, e2, alternative = "two.sided", h = 1, power = 2),
        error = function(e) return(NA)
      )
      p_values[i] <- ifelse(is.list(dm), dm$p.value, NA)
    } else {
      p_values[i] <- NA
    }
  }
  return(p_values)
}

dm_mixed_vs_causal <- compute_dm_tests(forecast_mixed, forecast_causal, actual_matrix, h)

# -----------------------------------------------------------------------------
# Combine RMSEs and DM p-values into a tidy data frame
# -----------------------------------------------------------------------------

rmse_df <- data.frame(
  horizon = 1:h,
  RMSE_mixed = rmse_mixed,
  RMSE_causal = rmse_causal,
  DM_mixed_vs_causal = dm_mixed_vs_causal
)

# Print RMSE and DM test comparison
print(rmse_df)


## SETAR AND MART ESTIMATION AND FORECASTING
# -----------------------------------------------------------------------------
# Grid search SETAR
# -----------------------------------------------------------------------------

library(pbmcapply)
source("MART.R")

# Set model parameters
thresholds <- seq(0.4, 1.2, by = 0.1) 
ds <- seq(1, 6)
p_C_max <- 4
p_NC_max <- 4

# Set number of cores based on OS
if (.Platform$OS.type == "windows") {
  n_cores <- 1
} else {
  n_cores <- 3 
  RNGkind("L'Ecuyer-CMRG")  # Safe parallel RNG
}

# Create a parameter grid of all combinations
param_grid <- expand.grid(
  threshold = thresholds,
  d = ds,
  #g = gammas,
  i = 0:p_C_max, 
  j = 0
)

# Function to run MART and get BIC value
run_model <- function(params) {
  t <- params$threshold
  d <- params$d
  #gamma <- params$g
  i <- params$i
  j <- params$j
  
  # Run MART model
  MART_d <- MART(data, NULL, i, j, t, d)
  bic_value <- information.criteria("MART", MART_d)
  return(data.frame(threshold=t, d=d, i=i, j=j, bic=bic_value))
}

# Use pbmclapply with a progress bar
bic_results <- pbmclapply(
  1:nrow(param_grid),
  function(idx) run_model(param_grid[idx, ]),
  mc.cores = n_cores 
)

# Combine all results into one data frame
ART_Q <- do.call(rbind, bic_results)

# Voeg de gamma-kolom toe aan het uiteindelijke data frame
#bic_smart_tdig_df$gamma <- param_grid$g


# Save the results
save(ART_Q, file = "BIC_ART_Q.RData")

View(ART_Q)


# -----------------------------------------------------------------------------
# Grid search MART
# -----------------------------------------------------------------------------

library(pbmcapply)
source("MART.R")

# Set model parameters
thresholds <- seq(0.4, 1.2, by = 0.1) 
ds <- seq(1, 6)
p_C_max <- 4
p_NC_max <- 4

# Set number of cores based on OS
if (.Platform$OS.type == "windows") {
  n_cores <- 1
} else {
  n_cores <- 3 
  RNGkind("L'Ecuyer-CMRG")  # Safe parallel RNG
}

# Create a parameter grid of all combinations
param_grid <- expand.grid(
  threshold = thresholds,
  d = ds,
  #g = gammas,
  i = 0:p_C_max, 
  j = 0:p_NC_max
)

# Function to run MART and get BIC value
run_model <- function(params) {
  t <- params$threshold
  d <- params$d
  #gamma <- params$g
  i <- params$i
  j <- params$j
  
  # Run MART model
  MART_d <- MART(data, NULL, i, j, t, d)
  bic_value <- information.criteria("MART", MART_d)
  return(data.frame(threshold=t, d=d, i=i, j=j, bic = bic_value))
}

# Use pbmclapply with a progress bar
bic_results <- pbmclapply(
  1:nrow(param_grid),
  function(idx) run_model(param_grid[idx, ]),
  mc.cores = n_cores
)

# Combine all results into one data frame
MART_Q <- do.call(rbind, bic_results)

# Voeg de gamma-kolom toe aan het uiteindelijke data frame
#bic_smart_tdig_df$gamma <- param_grid$g


# Save the results
save(MART_Q, file = "BIC_MART_Q.RData")

View(MART_Q)


source("MART.R")
load("inflation_df_monthly.RData")
library(pbmcapply)
library(forecast)
# -----------------------------------------------------------------------------
# Multi-horizon out-of-sample forecast evaluation
# -----------------------------------------------------------------------------

# Forecasting parameters
h <- 4         # Forecast horizon
N <- 15000        # Posterior draws
M <- 50          # MA truncation

# Model specifications
p_C_mart <- 2;  p_NC_mart <- 2    # Mixed MAR(2,2)
p_C_art <- 4; p_NC_art<- 0   # Purely causal SETAR(4,0)
c_mart <- 1.2
c_setar <- 0.9
d_mart <- 1
d_setar <- 1

# Define forecast evaluation window
data_series <- data
start_index <- 100
end_index <- length(data_series) - h
forecast_indices <- start_index:end_index

# -----------------------------------------------------------------------------
# Forecasting loop using the 2 model configurations
# -----------------------------------------------------------------------------

results_list <- pbmclapply(
  X = forecast_indices,
  FUN = function(t) {
    tryCatch({
      y_window <- data_series[1:t]
      
      # Call forecast.MART with correct parameters
      forecast_mart <- forecast.MART(
        y = y_window,
        p_C = p_C_mart,
        p_NC = p_NC_mart,
        c = c_mart,
        d = d_mart,
        h = h,
        M = M,
        N = N
      )
      
      # Make sure we get only the forecast component
      mart_forecast <- forecast_mart$forecast
      mart_defaulted <- forecast_mart$defaulted
      
      forecast_art <- forecast.MART(
        y = y_window,
        p_C = p_C_art,
        p_NC = p_NC_art,
        c = c_setar,
        d = d_setar,
        h = h,
        M = M,
        N = N
      )
      
      # Make sure we get only the forecast component
      art_forecast <- forecast_art$forecast
      art_defaulted <- forecast_art$defaulted
      
      actual <- data_series[(t + 1):(t + h)]
      
      return(list(
        mart = mart_forecast,         # Use mart_forecast directly here
        mart_defaulted = mart_defaulted, # Use mart_defaulted directly here
        art = art_forecast,           # Use art_forecast directly here
        art_defaulted = art_defaulted,   # Use art_defaulted directly here
        actual = actual
      ))
    }, error = function(e) {
      message(sprintf("Error at t = %d: %s", t, e$message))
      return(NULL)
    })
  },
  mc.cores = parallel::detectCores() - 1
)
# -----------------------------------------------------------------------------
# Organize forecast results into matrices
# -----------------------------------------------------------------------------

# Extract default flags
mart_quarterly_default_flags <- sapply(results_list, function(x) if (!is.null(x)) x$mart_defaulted else NA)
art_quarterly_default_flags <- sapply(results_list, function(x) if (!is.null(x)) x$art_defaulted else NA)

# Compute default percentages
pct_default_mart_quarterly <- mean(mart_quarterly_default_flags, na.rm = TRUE) * 100
pct_default_art_quarterly <- mean(art_quarterly_default_flags, na.rm = TRUE) * 100

# Safely extract components and skip NULLs
forecast_mart_quarterly <- do.call(rbind, lapply(results_list, function(x) {
  if (!is.null(x) && !is.null(x$mart)) return(x$mart)
  return(matrix(NA, nrow = 1, ncol = h))  # Fallback for failed cases
}))

forecast_art_quarterly <- do.call(rbind, lapply(results_list, function(x) {
  if (!is.null(x) && !is.null(x$art)) return(x$art)
  return(matrix(NA, nrow = 1, ncol = h))
}))

actual_matrix <- do.call(rbind, lapply(results_list, function(x) {
  if (!is.null(x) && !is.null(x$actual)) return(matrix(x$actual, nrow = 1))
  return(matrix(NA, nrow = 1, ncol = h))
}))

colnames(forecast_mart_quarterly) <- paste0("h", 1:h)
colnames(forecast_art_quarterly) <- paste0("h", 1:h)
colnames(actual_matrix) <- paste0("h", 1:h)

save(forecast_mart_quarterly, forecast_art_quarterly, actual_matrix, file = "forecast_quarterly_results.RData")

# -----------------------------------------------------------------------------
# Compute RMSE for each model across horizons
# -----------------------------------------------------------------------------

rmse <- function(forecast, actual) {
  sqrt(colMeans((forecast - actual)^2, na.rm = TRUE))
}

rmse_mart <- rmse(forecast_mart_quarterly, actual_matrix)
rmse_art <- rmse(forecast_art_quarterly, actual_matrix)


# -----------------------------------------------------------------------------
# Compute Diebold-Mariano test p-values
# -----------------------------------------------------------------------------

compute_dm_tests <- function(forecast1, forecast2, actual, h) {
  p_values <- numeric(h)
  for (i in 1:h) {
    e1 <- actual[, i] - forecast1[, i]
    e2 <- actual[, i] - forecast2[, i]
    valid <- complete.cases(e1, e2)
    e1 <- e1[valid]
    e2 <- e2[valid]
    
    if (length(e1) > 10) {
      dm <- tryCatch(
        dm.test(e1, e2, alternative = "two.sided", h = 1, power = 2),
        error = function(e) return(NA)
      )
      p_values[i] <- ifelse(is.list(dm), dm$p.value, NA)
    } else {
      p_values[i] <- NA
    }
  }
  return(p_values)
}

dm_mart_vs_causal <- compute_dm_tests(forecast_mart_quarterly, forecast_art_quarterly, actual_matrix, h)


# -----------------------------------------------------------------------------
# Combine RMSEs and DM p-values into a tidy data frame
# -----------------------------------------------------------------------------

rmse_df <- data.frame(
  horizon = 1:h,
  RMSE_mart = rmse_mart,
  RMSE_art = rmse_art,
  DM_mart_vs_art = dm_mart_vs_causal
)



# Print RMSE and DM test comparison
print(rmse_df)

