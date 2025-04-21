

# -----------------------------------------------------------------------------
# Multi-horizon out-of-sample forecast evaluation
# -----------------------------------------------------------------------------

# Forecasting parameters
h <- 6        # Forecast horizon
N <- 1000        # Posterior draws
M <- 50          # MA truncation
d <- 6
gamma <- 15
c <- median(inflation_df_monthly$inflationNonSA)

# Model specifications
p_C_smart <- 0;  p_NC_smart <- 2    # Mixed MAR(1,1)
p_C_sart <- 2; p_NC_sart <- 0   # Purely causal AR(12)

# Define forecast evaluation window
data_series <- inflation_df_monthly$inflationNonSA
start_index <- 100
end_index <- length(data_series) - h
forecast_indices <- start_index:end_index

# -----------------------------------------------------------------------------
# Forecasting loop using all 3 model configurations
# -----------------------------------------------------------------------------

results_list <- pbmclapply(
  X = forecast_indices,
  FUN = function(t) {
    tryCatch({
      y_window <- data_series[1:t]
      
      forecast_mart <- forecast.SMART(
        y = y_window,
        p_C = p_C_mart,
        p_NC = p_NC_mart,
        c = c,
        gamma = gamma,
        d = d,
        h = h,
        M = M,
        N = N
      )
      forecast_art <-  forecast_art <- forecast.MART(
        y = y_window,
        p_C = p_C_art,
        p_NC = p_NC_art,
        c = c,
        gamma = gamma,
        d = d,
        h = h,
        M = M,
        N = N
      )
      
      actual <- data_series[(t + 1):(t + h)]
      
      return(list(mart = forecast_smart, art = forecast_sart, actual = actual))
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

# Safely extract components and skip NULLs
forecast_smart <- do.call(rbind, lapply(results_list, function(x) {
  if (!is.null(x) && !is.null(x$smart)) return(x$smart)
  return(matrix(NA, nrow = 1, ncol = h))  # Fallback for failed cases
}))

forecast_sart <- do.call(rbind, lapply(results_list, function(x) {
  if (!is.null(x) && !is.null(x$sart)) return(x$sart)
  return(matrix(NA, nrow = 1, ncol = h))
}))

actual_matrix <- do.call(rbind, lapply(results_list, function(x) {
  if (!is.null(x) && !is.null(x$actual)) return(matrix(x$actual, nrow = 1))
  return(matrix(NA, nrow = 1, ncol = h))
}))

colnames(forecast_smart) <- paste0("h", 1:h)
colnames(forecast_sart) <- paste0("h", 1:h)
colnames(actual_matrix) <- paste0("h", 1:h)

# -----------------------------------------------------------------------------
# Compute RMSE for each model across horizons
# -----------------------------------------------------------------------------

rmse <- function(forecast, actual) {
  sqrt(colMeans((forecast - actual)^2, na.rm = TRUE))
}

rmse_smart <- rmse(forecast_smart, actual_matrix)
rmse_sart <- rmse(forecast_sart, actual_matrix)


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

dm_smart_vs_sart <- compute_dm_tests(forecast_smart, forecast_sart, actual_matrix, h)


# -----------------------------------------------------------------------------
# Combine RMSEs and DM p-values into a tidy data frame
# -----------------------------------------------------------------------------

rmse_df <- data.frame(
  horizon = 1:h,
  RMSE_smart = rmse_smart,
  RMSE_sart = rmse_sart,
  DM_smart_vs_sart = dm_smart_vs_sart
)

# Print RMSE and DM test comparison
print(rmse_df)