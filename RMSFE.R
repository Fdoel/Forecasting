load("Forecasting results/forecast_MART_results.RData") # Results from MART pseudo and GS
load("Forecasting results/forecast_MARTX_results.RData") # Results from MART X pseudo and GS
load("Forecasting results/forecast_ART_results.RData") # Results from ART GS
load("Forecasting results/forecast_ARTX_results.RData") # Results from ARX GS
load("Forecasting results/forecast_MAR_results.RData") # Results from MAR and AR
load("Forecasting results/forecast_ARX_results.RData") # Results from AR and MAR
load("Forecasting results/pct_default_MART.RData") # Percentage MART models defaulted
load("Forecasting results/pct_default_MART_X.RData") # Percentage MART models defaulted
load("Forecasting results/forecast_results_MAR_SA.RData") # Results for seasonally adjusted AR and MAR
load("Forecasting results/forecast_results_SA.RData") # Results for seasonally adjusted ART and MART
# Load default percentages for MART models

library(ggplot2)
library(tidyr)
library(zoo)

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

dm_marx_vs_artx <- compute_dm_tests(forecast_ARX,forecast_art_x,actual_matrix, h)
dm_marx_vs_martx <- compute_dm_tests(forecast_ARX,forecast_mart_x_grid,actual_matrix, h)
dm_martx_vs_artx <- compute_dm_tests(forecast_mart_x_grid,forecast_art_x,actual_matrix, h)

# For robustness seasonally adjusted
dm_SA_ar_vs_art <- compute_dm_tests(forecast_SA_art, forecast_SA_AR, actual_matrix, h)

# Post 2000
index_2000 <- 238 # The index in the forecasts that corresponds to the year 2000.
n <- nrow(actual_matrix) # Number of observations in the actual matrix
# Combine all observations after 2000 into a RMSFE dataframe
rmse_mar_grid_pseudo_2000 <- rmse(forecast_MAR[index_2000:n,], actual_matrix[index_2000:n,])
rmse_ar_2000 <- rmse(forecast_AR[index_2000:n,], actual_matrix[index_2000:n,])
rmse_mart_grid_2000 <- rmse(forecast_mart_grid[index_2000:n,], actual_matrix[index_2000:n,])
rmse_mart_pseudo_2000 <- rmse(forecast_mart_pseudo[index_2000:n,], actual_matrix[index_2000:n,])
rmse_mart_x_grid_2000 <- rmse(forecast_mart_x_grid[index_2000:n,], actual_matrix[index_2000:n,])
rmse_mart_x_pseudo_2000 <- rmse(forecast_mart_x_pseudo[index_2000:n,], actual_matrix[index_2000:n,])
rmse_art_2000 <- rmse(forecast_art[index_2000:n,], actual_matrix[index_2000:n,])
rmse_art_x_2000 <- rmse(forecast_art_x[index_2000:n,], actual_matrix[index_2000:n,])
rsme_arx_2000 <- rmse(forecast_ARX[index_2000:n,], actual_matrix[index_2000:n,])

# Combine into a tidy data frame
rmse_df_2000 <- data.frame(
  horizon = 1:h,
  RMSE_mar_grid_pseudo_2000 = rmse_mar_grid_pseudo_2000,
  RMSE_ar_2000 = rmse_ar_2000,
  RMSE_mart_grid_2000 = rmse_mart_grid_2000,
  RMSE_mart_pseudo_2000 = rmse_mart_pseudo_2000,
  RMSE_mart_x_grid_2000 = rmse_mart_x_grid_2000,
  RMSE_mart_x_pseudo_2000 = rmse_mart_x_pseudo_2000,
  RMSE_art_2000 = rmse_art_2000,
  RMSE_art_x_2000 = rmse_art_x_2000,
  RMSE_arx_2000 = rsme_arx_2000
)

# Print RMSE and DM test results for post 2000
print(rmse_df_2000)

# For robustness post 2000
dm_SA_ar_vs_mar_2000 <- compute_dm_tests(forecast_MAR[index_2000:n,], forecast_AR[index_2000:n,], actual_matrix[index_2000:n,], h)
dm_SA_ar_vs_marx_2000 <- compute_dm_tests(forecast_ARX[index_2000:n,], forecast_AR[index_2000:n,], actual_matrix[index_2000:n,], h)
dm_SA_mart_vs_art_2000 <- compute_dm_tests(forecast_mart_grid[index_2000:n,], forecast_art[index_2000:n,], actual_matrix[index_2000:n,], h)

# Extract 12-month ahead forecasts and actuals
mart_12 <- forecast_mart_grid[, 12]
art_12 <- forecast_art[, 12]
actual_12 <- actual_matrix[, 12]
mar_12 <-forecast_MAR[,12]
ar_12 <-forecast_AR[,12]
marx_12 <- forecast_ARX[,12]
artx_12 <- forecast_art_x[,12]
martx_12 <- forecast_mart_x_grid[,12]

# Set time index
start_year <- 1959
n_months <- 787 - 12
dates <- seq(ymd(paste0(start_year, "-06-01")), by = "month", length.out = n_months)

# Forecast starts from index 250
forecast_start <- 250
forecast_dates <- as.yearmon(dates[forecast_start:n_months])

# Combine into a data frame
df <- data.frame(
  Date = forecast_dates,
  Actual = actual_12,
  MART_GS = mart_12,
  SETAR = art_12
)

df_long <- pivot_longer(df, cols = c("Actual", "MART_GS", "SETAR"),
                        names_to = "Series", values_to = "Value")

# Plot
ggplot(df_long, aes(x = Date, y = Value, color = Series)) +
  geom_line(linewidth = 1) +
  labs(title = "12-Month Ahead Forecast vs Actual",
       x = "Date", y = "Monthly Inflation", color = "Legend") +
  scale_x_yearmon(format = "%Y-%m", n = 10) +
  scale_color_manual(values = c("Actual" = "black", "MART_GS" = "blue", "SETAR" = "red")) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Combine into a data frame
df_2 <- data.frame(
  Date = forecast_dates,
  Actual = actual_12,
  MAR = mar_12,
  AR = ar_12
)
df_long_2 <- pivot_longer(df_2, cols = c("Actual", "MAR", "AR"),
                        names_to = "Series", values_to = "Value")

# Plot
ggplot(df_long_2, aes(x = Date, y = Value, color = Series)) +
  geom_line(linewidth = 1) +
  labs(title = "12-Month Ahead Forecast vs Actual",
       x = "Date", y = "Monthly Inflation", color = "Legend") +
  scale_x_yearmon(format = "%Y-%m", n = 10) +
  scale_color_manual(values = c("Actual" = "black", "MAR" = "blue", "AR" = "red")) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Combine into a data frame (note the column names match the legend labels)
df_3 <- data.frame(
  Date = forecast_dates,
  Actual = actual_12,
  MARX = marx_12,
  `ART_X` = artx_12,
  `MART_X` = martx_12
)

# Pivot to long format
df_long_3 <- pivot_longer(df_3, 
                          cols = c("Actual", "MARX", "ART_X", "MART_X"),
                          names_to = "Series", values_to = "Value")

# Plot
ggplot(df_long_3, aes(x = Date, y = Value, color = Series)) +
  geom_line(linewidth = 1) +
  labs(title = "12-Month Ahead Forecast vs Actual",
       x = "Date", y = "Monthly Inflation", color = "Legend") +
  scale_x_yearmon(format = "%Y-%m", n = 10) +
  scale_color_manual(values = c("Actual" = "black", 
                                "MARX" = "blue", 
                                "ART_X" = "red", 
                                "MART_X" = "green")) +
  theme_minimal() +
  theme(legend.position = "bottom")
