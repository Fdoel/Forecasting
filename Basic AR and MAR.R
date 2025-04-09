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

# -----------------------------------------------------------------------------
# Model order selection for MAR model using information criteria
# -----------------------------------------------------------------------------

# Initial model diagnostic: scan up to lag 18 with 10% significance level
marx(inflation_df_monthly$inflationNonSA, NULL, p_max = 18, sig_level = 0.1)

# Define maximum lag orders for causal (p_C) and noncausal (p_NC) components
p_C_max <- 12
p_NC_max <- 12

# Initialize matrices to store information criteria and model diagnostics
LL <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)
AIC <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)
BIC <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)
HQ <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)
N <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)

# Loop over all possible (p_C, p_NC) combinations and compute info criteria
for (i in 0:p_C_max) {
  for (j in 0:p_NC_max) {
    marx_loop <- marx.t(inflation_df_monthly$inflationNonSA, NULL, p_C = i, p_NC = j)
    information <- information.criteria(type = "MARX", marx_loop)
    LL[(i+1),(j+1)] <- information$loglikelihood
    AIC[(i+1),(j+1)] <- information$aic
    BIC[(i+1),(j+1)] <- information$bic
    HQ[(i+1),(j+1)] <- information$hq
    N[(i+1),(j+1)] <- information$n
  }
}

# -----------------------------------------------------------------------------
# Residual diagnostics: test for independence of squared residuals (ARCH test)
# -----------------------------------------------------------------------------

# Fit a 12-lag AR model to the inflation series
model_ar12 <- Arima(inflation_df_monthly$inflationNonSA, order = c(12, 0, 0))
resids_ar12 <- model_ar12$residuals  # Extract residuals

# Step 2: Square the residuals for ARCH effect detection
resids_sq <- resids_ar12^2

# Step 3: Create lag matrix of squared residuals (lags 1 through m)
m <- 12
X <- embed(resids_sq, m + 1)
y <- X[, 1]            # Current value
X_lags <- X[, -1]      # Lagged squared residuals

# Step 4: Regress current squared residual on its lags
model_test <- lm(y ~ X_lags)

# Step 5: Test for joint significance of lag coefficients (H0: no ARCH effect)
test_statistic <- summary(model_test)$r.squared * length(y)
p_value <- pchisq(test_statistic, df = m, lower.tail = FALSE)

# Step 6: Output results of the chi-squared test
cat("Chi-squared test statistic:", test_statistic, "\n")
cat("p-value:", p_value, "\n")

# -----------------------------------------------------------------------------
# Model estimation: compare AR(12) vs mixed MAR(1,11)
# -----------------------------------------------------------------------------

# Estimate AR(12) in MAR framework (purely causal)
ar12 <- marx.t(inflation_df_monthly$inflationNonSA, NULL, p_C = 12, p_NC = 0)
ar12_crit <- information.criteria("MARX", ar12)  # Info criteria for AR(12)

# Estimate MAR(1,11): mixed causal-noncausal model
mar_final <- marx.t(inflation_df_monthly$inflationNonSA, NULL, p_C = 1, p_NC = 11)
mar111_crit <- information.criteria("MARX", mar_final)  # Info criteria for MAR(1,11)

# -----------------------------------------------------------------------------
# Multi-horizon out-of-sample forecast evaluation
# -----------------------------------------------------------------------------

# Forecasting parameters
h <- 24          # Forecast horizon
N <- 1000        # Posterior draws
M <- 50          # Simulations per draw

# Model specifications
p_C_mixed <- 1;  p_NC_mixed <- 11    # Mixed MAR(1,11)
p_C_causal <- 12; p_NC_causal <- 0   # Purely causal AR(12)
p_C_mid <- 11; p_NC_mid <- 1         # Asymmetric MAR(11,1)

# Define forecast evaluation window
data_series <- inflation_df_monthly$inflationNonSA
start_index <- 100
end_index <- length(data_series) - h
forecast_indices <- start_index:end_index

# Parallel forecast loop using all 3 model configurations
results_list <- pbmclapply(
  X = forecast_indices,
  FUN = function(t) {
    y_window <- data_series[1:t]  # Expand window up to time t
    
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
    
    forecast_mid <- forecast.marx(
      y = y_window,
      p_C = p_C_mid,
      p_NC = p_NC_mid,
      h = h,
      M = M,
      N = N
    )
    
    actual <- data_series[(t + 1):(t + h)]  # True values to compare with forecasts
    
    return(list(
      mixed = forecast_mixed,
      causal = forecast_causal,
      mid = forecast_mid,
      actual = actual
    ))
  },
  mc.cores = parallel::detectCores() - 1  # Use all cores minus one
)

# -----------------------------------------------------------------------------
# Organize forecast results into matrices
# -----------------------------------------------------------------------------

# Combine lists into forecast matrices
forecast_mixed <- do.call(rbind, lapply(results_list, `[[`, "mixed"))
forecast_causal <- do.call(rbind, lapply(results_list, `[[`, "causal"))
forecast_mid <- do.call(rbind, lapply(results_list, `[[`, "mid"))
actual_matrix <- do.call(rbind, lapply(results_list, `[[`, "actual"))

# Label forecast horizons
colnames(forecast_mixed) <- paste0("h", 1:h)
colnames(forecast_causal) <- paste0("h", 1:h)
colnames(forecast_mid) <- paste0("h", 1:h)
colnames(actual_matrix) <- paste0("h", 1:h)

# -----------------------------------------------------------------------------
# Compute RMSE for each model across horizons
# -----------------------------------------------------------------------------

# RMSE calculation function
rmse <- function(forecast, actual) {
  sqrt(colMeans((forecast - actual)^2, na.rm = TRUE))
}

# Compute RMSEs
rmse_mixed <- rmse(forecast_mixed, actual_matrix)
rmse_causal <- rmse(forecast_causal, actual_matrix)
rmse_mid <- rmse(forecast_mid, actual_matrix)

# Combine RMSEs into a tidy data frame
rmse_df <- data.frame(
  horizon = 1:h,
  RMSE_mixed = rmse_mixed,
  RMSE_causal = rmse_causal,
  RMSE_mid = rmse_mid
)

# Print RMSE comparison for each horizon
print(rmse_df)