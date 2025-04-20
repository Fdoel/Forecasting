load("forecast_MART_results.RData") # Results from MART pseudo and GS
load("forecast_MARTX_results.RData") # Results from MART X pseudo and GS
load("forecast_ART_results.RData") # Results from ART GS
load("forecast_ARTX_results.RData") # Results from ARX GS

# -----------------------------------------------------------------------------
# Compute RMSE for each model across horizons
# -----------------------------------------------------------------------------

rmse <- function(forecast, actual) {
  sqrt(colMeans((forecast - actual)^2, na.rm = TRUE))
}

rmse_mart_grid <- rmse(forecast_mart_grid, actual_matrix)
rmse_mart_pseudo <- rmse(forecast_mart_pseudo, actual_matrix)
rmse_mart_x_grid <- rmse(forecast_mart_x_grid, actual_matrix)
rmse_mart_x_pseudo <- rmse(forecast_mart_x_pseudo, actual_matrix)

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
# Combine RMSEs and DM p-values into a tidy data frame
# -----------------------------------------------------------------------------
h <- 12
rmse_df <- data.frame(
  horizon = 1:h,
  RMSE_mart_grid = rmse_mart_grid,
  RMSE_mart_pseudo = rmse_mart_pseudo,
  RMSE_mart_x_grid = rmse_mart_x_grid,
  RMSE_mart_x_pseudo = rmse_mart_x_pseudo
)



# Print RMSE and DM test comparison
print(rmse_df)
