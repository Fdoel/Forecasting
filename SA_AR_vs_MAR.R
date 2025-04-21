
# -----------------------------------------------------------------------------
# Multi-horizon out-of-sample forecast evaluation
# -----------------------------------------------------------------------------

# Forecasting parameters
h <- 12          # Forecast horizon
N <- 15000        # Posterior draws
M <- 50          # MA truncation

# Model specifications
p_C_SA_mar <- 1;  p_NC_SA_mar <- 9   
p_C_SA_ar <- 12; p_NC_SA_ar<- 0   


# Define forecast evaluation window
data_series <- inflation_df_monthly$inflationSA
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
    
    forecast_mar_SA <- forecast.marx(
      y = y_window,
      p_C = p_C_SA_mar,
      p_NC = p_NC_SA_mar,
      h = h,
      M = M,
      N = N
    )
    
    forecast_ar_SA <- forecast.marx(
      y = y_window,
      p_C = p_C_SA_ar,
      p_NC = p_NC_SA_ar,
      h = h,
      M = M,
      N = N
    )

    actual <- data_series[(t + 1):(t + h)]  # Future observed values
    
    return(list(
      MAR_SA = forecast_SA_MAR,
      AR_SA = forecast_SA_AR,
      actual = actual
    ))
  },
  mc.cores = parallel::detectCores() - 1  # Use all but one core
)

# -----------------------------------------------------------------------------
# Organize forecast results into matrices
# -----------------------------------------------------------------------------

forecast_SA_MAR <- do.call(rbind, lapply(results_list, `[[`, "MAR_SA"))
forecast_SA_AR <- do.call(rbind, lapply(results_list, `[[`, "AR_SA"))
actual_matrix <- do.call(rbind, lapply(results_list, `[[`, "actual"))

colnames(forecast_SA_MAR) <- paste0("h", 1:h)
colnames(forecast_SA_AR) <- paste0("h", 1:h)
colnames(actual_matrix) <- paste0("h", 1:h)

save(forecast_SA_MAR, forecast_SA_AR, actual_matrix, file = "forecast_results_MAR_SA.RData")

# -----------------------------------------------------------------------------
# Compute RMSE for each model across horizons
# -----------------------------------------------------------------------------

rmse <- function(forecast, actual) {
  sqrt(colMeans((forecast - actual)^2, na.rm = TRUE))
}

rmse_SA_MAR <- rmse(forecast_SA_MAR, actual_matrix)
rmse_SA_AR <- rmse(forecast_SA_AR, actual_matrix)


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

dm_SA_MAR_vs_AR <- compute_dm_tests(forecast_SA_MAR, forecast_SA_AR, actual_matrix, h)

# -----------------------------------------------------------------------------
# Combine RMSEs and DM p-values into a tidy data frame
# -----------------------------------------------------------------------------

rmse_df <- data.frame(
  horizon = 1:h,
  RMSE_SA_MAR = rmse_SA_MAR,
  RMSE_SA_AR = rmse_SA_AR,
  DM_SA_MAR_vs_AR = dm_SA_MAR_vs_AR

)

# Print RMSE and DM test comparison
print(rmse_df)