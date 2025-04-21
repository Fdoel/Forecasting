load("Forecasting results/forecast_MART_results.RData") # Results from MART pseudo and GS
load("Forecasting results/forecast_MARTX_results.RData") # Results from MART X pseudo and GS
load("Forecasting results/forecast_ART_results.RData") # Results from ART GS
load("Forecasting results/forecast_ARTX_results.RData") # Results from ARX GS
load("Forecasting results/forecast_MAR_results.RData") # Results from MAR and AR
load("Forecasting results/forecast_ARX_results.RData") # Results from AR and MAR
load("Forecasting results/pct_default_MART.RData") # Percentage MART models defaulted
load("Forecasting results/pct_default_MART_X.RData") # Percentage MART models defaulted
# Load default percentages for MART models


# -----------------------------------------------------------------------------
# Compute RMSE for each model across horizons
# -----------------------------------------------------------------------------

rmse <- function(forecast, actual) {
  sqrt(colMeans((forecast - actual)^2, na.rm = TRUE))
}

rmse_mar_grid_pseudo <- rmse(forecast_MAR, actual_matrix)
rmse_ar <- rmse(forecast_AR, actual_matrix)
rmse_mart_grid <- rmse(forecast_mart_grid, actual_matrix)
rmse_mart_pseudo <- rmse(forecast_mart_pseudo, actual_matrix)
rmse_mart_x_grid <- rmse(forecast_mart_x_grid, actual_matrix)
rmse_mart_x_pseudo <- rmse(forecast_mart_x_pseudo, actual_matrix)
rmse_art <- rmse(forecast_art, actual_matrix)
rmse_art_x <- rmse(forecast_art_x, actual_matrix)
rsme_arx <- rmse(forecast_ARX, actual_matrix)

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


# -----------------------------------------------------------------------------
# Combine RMSEs into a tidy data frame
# -----------------------------------------------------------------------------
h <- 12
rmse_df <- data.frame(
  horizon = 1:h,
  RMSE_mar_grid_pseudo = rmse_mar_grid_pseudo,
  RMSE_ar = rmse_ar,
  RMSE_mart_grid = rmse_mart_grid,
  RMSE_mart_pseudo = rmse_mart_pseudo,
  RMSE_mart_x_grid = rmse_mart_x_grid,
  RMSE_mart_x_pseudo = rmse_mart_x_pseudo,
  RMSE_art = rmse_art,
  RMSE_art_x = rmse_art_x,
  RMSE_arx = rsme_arx
)
# Print RMSE an
print(rmse_df)

# Compute select DM test results

# MART pseudo vs MART grid
dm_test_mart_pseudo_grid <- compute_dm_tests(forecast_mart_pseudo, forecast_mart_grid, actual_matrix, h)
# Print with horizon
dm_test_mart_pseudo_grid_df <- data.frame(
  horizon = 1:h,
  DM_p_value = dm_test_mart_pseudo_grid
)
print(dm_test_mart_pseudo_grid_df)

# MART grid vs ART
dm_test_mart_grid_art <- compute_dm_tests(forecast_mart_grid, forecast_art, actual_matrix, h)
# Print with horizon
dm_test_mart_grid_art_df <- data.frame(
  horizon = 1:h,
  DM_p_value = dm_test_mart_grid_art
)
print(dm_test_mart_grid_art_df)


