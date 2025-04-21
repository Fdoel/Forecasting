
# -----------------------------------------------------------------------------
# Multi-horizon out-of-sample forecast evaluation
# -----------------------------------------------------------------------------

# Forecasting parameters
h <- 12         # Forecast horizon
N <- 15000        # Posterior draws
M <- 50          # MA truncation

# Model specifications
p_C_SA_mart <- 5;  p_NC_SA_mart <- 1    # Mixed MAR(1,1)
p_C_SA_art <- 4; p_NC_SA_art<- 0   # Purely causal SETAR(2,0)
c_SA_mart <- 0.1
c_SA_art <- 0.6
d_SA_mart <- 1
d_SA_art <- 1


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
    tryCatch({
      y_window <- data_series[1:t]

      # Call forecast.MART with correct parameters
      forecast_mart <- forecast.MART(
        y = y_window,
        p_C = p_C_SA_mart,
        p_NC = p_NC_SA_mart,
        c = c_SA_mart,
        d = d_SA_mart,
        h = h,
        M = M,
        N = N
      )
      
      # Make sure we get only the forecast component
      SA_mart_forecast <- forecast_mart$forecast
      SA_mart_defaulted <- forecast_mart$defaulted
      
      forecast_art <- forecast.MART(
        y = y_window,
        p_C = p_C_SA_art,
        p_NC = p_NC_SA_art,
        c = c_SA_art,
        d = d_SA_art,
        h = h,
        M = M,
        N = N
      )
      
      # Make sure we get only the forecast component
      SA_art_forecast <- forecast_art$forecast
      SA_art_defaulted <- forecast_art$defaulted
      
      actual <- data_series[(t + 1):(t + h)]
      
      return(list(
        SA_mart = SA_mart_forecast,         # Use mart_forecast directly here
        SA_mart_defaulted = SA_mart_defaulted, # Use mart_defaulted directly here
        SA_art = SA_art_forecast,           # Use art_forecast directly here
        SA_art_defaulted = SA_art_defaulted,   # Use art_defaulted directly here
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
SA_mart_default_flags <- sapply(results_list, function(x) if (!is.null(x)) x$SA_mart_defaulted else NA)
SA_art_default_flags <- sapply(results_list, function(x) if (!is.null(x)) x$SA_art_defaulted else NA)

# Compute default percentages
pct_default_SA_mart <- mean(SA_mart_default_flags, na.rm = TRUE) * 100
pct_default_SA_art <- mean(SA_art_default_flags, na.rm = TRUE) * 100

# Safely extract components and skip NULLs
forecast_SA_mart <- do.call(rbind, lapply(results_list, function(x) {
  if (!is.null(x) && !is.null(x$SA_mart)) return(x$SA_mart)
  return(matrix(NA, nrow = 1, ncol = h))  # Fallback for failed cases
}))

forecast_SA_art <- do.call(rbind, lapply(results_list, function(x) {
  if (!is.null(x) && !is.null(x$SA_art)) return(x$SA_art)
  return(matrix(NA, nrow = 1, ncol = h))
}))

actual_matrix <- do.call(rbind, lapply(results_list, function(x) {
  if (!is.null(x) && !is.null(x$actual)) return(matrix(x$actual, nrow = 1))
  return(matrix(NA, nrow = 1, ncol = h))
}))

colnames(forecast_SA_mart) <- paste0("h", 1:h)
colnames(forecast_SA_art) <- paste0("h", 1:h)
colnames(actual_matrix) <- paste0("h", 1:h)

save(forecast_SA_mart, forecast_SA_art, actual_matrix, file = "forecast_results_SA.RData")

# -----------------------------------------------------------------------------
# Compute RMSE for each model across horizons
# -----------------------------------------------------------------------------

rmse <- function(forecast, actual) {
  sqrt(colMeans((forecast - actual)^2, na.rm = TRUE))
}

rmse_SA_mart <- rmse(forecast_SA_mart, actual_matrix)
rmse_SA_art <- rmse(forecast_SA_art, actual_matrix)


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

dm_SA_mart_vs_art <- compute_dm_tests(forecast_SA_mart, forecast_SA_art, actual_matrix, h)


# -----------------------------------------------------------------------------
# Combine RMSEs and DM p-values into a tidy data frame
# -----------------------------------------------------------------------------

rmse_df <- data.frame(
  horizon = 1:h,
  RMSE_SA_mart = rmse_SA_mart,
  RMSE_SA_art = rmse_SA_art,
  DM_SA_mart_vs_art = dm_SA_mart_vs_art
)



# Print RMSE and DM test comparison
print(rmse_df)