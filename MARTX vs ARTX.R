source("MART.R")
load("inflation_df_monthly.RData")
library(pbmcapply)
library(forecast)
library(tidyverse)
library(ggplot2)

# -----------------------------------------------------------------------------
# Multi-horizon out-of-sample forecast evaluation
# -----------------------------------------------------------------------------

# Forecasting parameters
h <- 12         # Forecast horizon
N <- 15000      # Posterior draws
M <- 50         # MA truncation

# Model specifications
p_C_mart_x_grid <- 1; p_NC_mart_x_grid <- 3       # Purely causal mart_pseudo(2,0)
p_C_mart_x_pseudo <- 1; p_NC_mart_x_pseudo <- 4      # Purely causal mart_pseudo(2,0)
c_mart_x_pseudo <- median(inflation_df_monthly$inflationNonSA)
c_mart_x_grid <- 0.6
d_mart_x_pseudo <- 1
d_mart_x_grid <- 1

# Define exogenous regressors
exo <- cbind(inflation_df_monthly$ldGS1, inflation_df_monthly$dRPI, inflation_df_monthly$dRETAIL)

# Define forecast evaluation window
data_series <- inflation_df_monthly$inflationNonSA
start_index <- 250
end_index <- length(data_series) - h
forecast_indices <- start_index:end_index

# -----------------------------------------------------------------------------
# Forecasting loop using both model configurations
# -----------------------------------------------------------------------------

results_list <- pbmclapply(
  X = forecast_indices,
  FUN = function(t) {
    tryCatch({
      y_window <- data_series[1:t]
      

      x_window <- exo[1:t,1:3]
      # Fit AR(1) models to each regressor
      GS1_model <- arima(exo[1:t, 1], order = c(1,0,0))
      dRPI_model <- arima(exo[1:t, 2], order = c(1,0,0))
      dRETAIL_model <- arima(exo[1:t, 3], order = c(1,0,0))
      GS1_for <- predict(GS1_model, n.ahead = M)$pred
      dRPI_for <- predict(dRPI_model, n.ahead = M)$pred
      dRETAIL_for <- predict(dRETAIL_model, n.ahead = M)$pred

      forecast_mart_x_pseudo <- forecast.MART(
        y = y_window,
        X = x_window,
        p_C = p_C_mart_x_pseudo,
        p_NC = p_NC_mart_x_pseudo,
        c = c_mart_x_pseudo,
        d = d_mart_x_pseudo,
        X.for = cbind(GS1_for, dRPI_for, dRETAIL_for),
        h = h,
        M = M,
        N = N
      )
      
      forecast_mart_x_grid <- forecast.MART(
        y = y_window,
        X = x_window,
        p_C = p_C_mart_x_grid,
        p_NC = p_NC_mart_x_grid,
        c = c_mart_x_grid,
        d = d_mart_x_grid,
        X.for = cbind(GS1_for, dRPI_for, dRETAIL_for),
        h = h,
        M = M,
        N = N
      )
      actual <- data_series[(t + 1):(t + h)]
      
      return(list(
        mart_x_pseudo = forecast_mart_x_pseudo$forecast,       
        mart_x_pseudo_defaulted = forecast_mart_x_pseudo$defaulted,
        mart_grid = forecast_mart_x_grid$forecast,
        mart_grid_defaulted = forecast_mart_x_grid$defaulted,
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

mart_x_pseudo_default_flags <- sapply(results_list, function(x) if (!is.null(x)) x$mart_x_pseudo_defaulted else NA)

pct_default_mart_x_pseudo <- mean(mart_x_pseudo_default_flags, na.rm = TRUE) * 100

forecast_mart_x_pseudo <- do.call(rbind, lapply(results_list, function(x) {
  if (!is.null(x) && !is.null(x$mart_x_pseudo)) return(x$mart_x_pseudo)
  return(matrix(NA, nrow = 1, ncol = h))
}))

mart_grid_default_flags <- sapply(results_list, function(x) if (!is.null(x)) x$mart_grid_defaulted else NA)

pct_default_mart_x_grid <- mean(mart_grid_default_flags, na.rm = TRUE) * 100

forecast_mart_x_grid <- do.call(rbind, lapply(results_list, function(x) {
  if (!is.null(x) && !is.null(x$mart_grid)) return(x$mart_grid)
  return(matrix(NA, nrow = 1, ncol = h))
}))

actual_matrix <- do.call(rbind, lapply(results_list, function(x) {
  if (!is.null(x) && !is.null(x$actual)) return(matrix(x$actual, nrow = 1))
  return(matrix(NA, nrow = 1, ncol = h))
}))

colnames(forecast_mart_x_pseudo) <- paste0("h", 1:h)
colnames(forecast_mart_x_grid) <- paste0("h", 1:h)
colnames(actual_matrix) <- paste0("h", 1:h)


save(forecast_mart_x_pseudo, forecast_mart_x_grid, actual_matrix, file = "forecast_x_pseudo3_results.RData")
save(pct_default_mart_x_pseudo, pct_default_mart_x_grid, file = "pct_default_MART_X.RData")

# -----------------------------------------------------------------------------
# Compute RMSE for each model across horizons
# -----------------------------------------------------------------------------

rmse <- function(forecast, actual) {
  sqrt(colMeans((forecast - actual)^2, na.rm = TRUE))
}

rmse_mart_pseudo <- rmse(forecast_mart_pseudo, actual_matrix)
rmse_mart_grid <- rmse(forecast_mart_grid, actual_matrix)

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

dm_test_results <- compute_dm_tests(forecast_mart_pseudo, forecast_mart_grid, actual_matrix, h)

# -----------------------------------------------------------------------------
# Combine RMSEs and DM p-values into a tidy data frame
# -----------------------------------------------------------------------------

rmse_df <- data.frame(
  horizon = 1:h,
  RMSE_mart_pseudo = rmse_mart_pseudo,
  RMSE_mart_grid = rmse_mart_grid,
  DM_p_value = dm_test_results
)

# Print RMSE and DM test comparison
print(rmse_df)
