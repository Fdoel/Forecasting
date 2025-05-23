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
load("inflation_df_monthly.RData")  # Load the inflation dataset

set.seed(20240421)

# -----------------------------------------------------------------------------
# Model order selection for MAR model using information criteria
# -----------------------------------------------------------------------------

# Initial model diagnostic: scan up to lag 18 with 10% significance level
#marx(inflation_df_monthly$inflationNonSA, NULL, p_max = 18, sig_level = 0.1)

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

## -----------------------------------------------------------------------------
# Residual diagnostics: test for independence of AR(p) residuals (Hecq et al. 2016) and test for no serial correlation (MARX package paper of HEcq et al.)
# -----------------------------------------------------------------------------

# Fit a 12-lag AR model to the inflation series
model_ar12 <- Arima(inflation_df_monthly$inflationNonSA, order = c(12, 0, 0))
resids_ar12 <- model_ar12$residuals  # Extract residuals

## - short intermezzo for residuals plotting
# Extract numeric values
resids_vals <- as.numeric(resids_ar12)

# Calculate mean and standard deviation
mu <- mean(resids_vals)
sigma <- sd(resids_vals)

# Create histogram with density
hist(resids_vals, breaks = 30, probability = TRUE,
     col = "lightblue", border = "white",
     main = "",
     xlab = "Residuals")

# Overlay normal distribution using sample mean and sd
curve(dnorm(x, mean = mu, sd = sigma),
      col = "red", lwd = 2, add = TRUE)

# Optionally, add legend
legend("topright", legend = c("Normal Density"), col = "red", lwd = 2, cex = 0.8)
## ----------------

# Step 2: Square the residuals for use as regressors
resids_sq <- resids_ar12^2

# Step 3: Create lag matrix manually
m <- 12
n <- length(resids_ar12)

# Create the response variable y (residuals from t = m+1 to n)
y <- resids_ar12[(m + 1):n]

# Create lagged squared residuals matrix
X_lags <- matrix(NA, nrow = n - m, ncol = m)
for (i in 1:m) {
  X_lags[, i] <- resids_sq[(m + 1 - i):(n - i)]
}

# Step 4: Regress current residual on lagged squared residuals
model_test <- lm(y ~ X_lags)

# Step 5: Test for joint significance of lag coefficients (H0: residuals are i.i.d.)
test_statistic <- summary(model_test)$r.squared * length(y)

# Step 6 simulate critical value using scaled t distribution
alpha <- 0.05
n_sim <- 10000
test_statistic_sim <- rep(NA, n_sim)
for(i in 1:n_sim) {
  e <- (stats::rt(length(y), df = 5)*sigma)
  e_2 <- e^2
  m <- 12
  X_lags <- matrix(NA, nrow = n - m, ncol = m)
  for (j in 1:m) {
    X_lags[, j] <- e_2[(m + 1 - j):(n - j)]
  }
  model_test_sim <- lm(e ~ X_lags)
  test_statistic_sim[i] <- summary(model_test_sim)$r.squared * length(y)
}

# Get the 95th percentile for the critical value
critical_value <- quantile(test_statistic_sim, 1 - alpha)

# Step 6: Compare test statistic with critical value
if (test_statistic > critical_value) {
  cat("Reject H0: residuals are not i.i.d.\n")
} else {
  cat("Fail to reject H0: residuals are i.i.d.\n")
}
# Step 7: Perfome Ljung-Box test
Box.test(resids_ar12, lag = 12, type = "Ljung-Box")

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
h <- 12          # Forecast horizon
N <- 15000        # Posterior draws
M <- 50          # MA truncation

# Model specifications
p_C_mixed <- 1;  p_NC_mixed <- 11    # Mixed MAR(1,11)
p_C_causal <- 12; p_NC_causal <- 0   # Purely causal AR(12)
p_C_mid <- 11; p_NC_mid <- 1         # Asymmetric MAR(11,1)

# Define forecast evaluation window
data_series <- inflation_df_monthly$inflationNonSA
start_index <- 250
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
    
    forecast_mid <- forecast.marx(
      y = y_window,
      p_C = p_C_mid,
      p_NC = p_NC_mid,
      h = h,
      M = M,
      N = N
    )
    
    actual <- data_series[(t + 1):(t + h)]  # Future observed values
    
    return(list(
      mixed = forecast_mixed,
      causal = forecast_causal,
      mid = forecast_mid,
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

# -----------------------------------------------------------------------------
# Compute RMSE for each model across horizons
# -----------------------------------------------------------------------------

rmse <- function(forecast, actual) {
  sqrt(colMeans((forecast - actual)^2, na.rm = TRUE))
}

rmse_mixed <- rmse(forecast_mixed, actual_matrix)
rmse_causal <- rmse(forecast_causal, actual_matrix)
rmse_mid <- rmse(forecast_mid, actual_matrix)

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
dm_mixed_vs_mid    <- compute_dm_tests(forecast_mixed, forecast_mid, actual_matrix, h)
dm_causal_vs_mid   <- compute_dm_tests(forecast_causal, forecast_mid, actual_matrix, h)

# -----------------------------------------------------------------------------
# Combine RMSEs and DM p-values into a tidy data frame
# -----------------------------------------------------------------------------

rmse_df <- data.frame(
  horizon = 1:h,
  RMSE_mixed = rmse_mixed,
  RMSE_causal = rmse_causal,
  RMSE_mid = rmse_mid,
  DM_mixed_vs_causal = dm_mixed_vs_causal,
  DM_mixed_vs_mid = dm_mixed_vs_mid,
  DM_causal_vs_mid = dm_causal_vs_mid
)

# Print RMSE and DM test comparison
print(rmse_df)