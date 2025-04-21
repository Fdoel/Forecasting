library(lubridate)

# Horizon of h=1
load("Forecasting results/AR and MAR forecasting with many draws.RData")
Forecast_AR <-  forecast_causal
Forecast_MAR <- forecast_mixed
Forecast_AR_1 <- Forecast_AR[,1]
Forecast_MAR_1 <- Forecast_MAR[,1]

load("Forecasting results/AR-X and MAR-X forecasting.RData")
Forecast_MAR_X <- forecast_causal
Forecast_MAR_X_1 <- Forecast_MAR_X[,1]

identical(Forecast_MAR, Forecast_MAR_X)

load("Forecasting results/forecast_ART_results.RData")
forecast_ART <- forecast_art[,1]

#load("~forecast_x_results.RData")
#forecast_ART_X <- forecast_art_x[,1]
#forecast_MART_X_GS <- forecast_mart_x[,1]
#RMSE_ART_X <- 

load("Forecasting results/forecast_MART_results.RData")
forecast_MART_GS_1 <- forecast_mart_grid[,1]

load("Forecasting results/forecast_x_pseudo_results.RData")
forecast_MART_x_pseudo_1 <- forecast_mart_x_pseudo[,1]
forecast_MART_pseudo_1 <- forecast_mart_pseudo[,1]

load("Forecasting results/forecast_MARTX_results.RData")
forecast_MART_X_GS_1 <- forecast_mart_x_grid[,1]

load("Forecasting results/forecast_ARTX_results.RData")
forecast_ART_x_1 <- forecast_art_x[,1]

load("Forecasting results/forecast_ART_results.RData")
forecast_ART_1 <- forecast_art[,1]


# Setup the actuals
start_year <- 1959
n_months <- 775  # as you defined earlier
dates <- seq(ymd(paste0(start_year, "-06-01")), by = "month", length.out = n_months)
forecast_start <- 250
forecast_dates <- dates[forecast_start:n_months]
actual <- actual_matrix[,1]
actual_forecast <- actual[forecast_start:n_months]

# Function to compute RMSFE from Jan 2000 onward
forecasting_2000 <- function(forecast, name) {
  # Align forecast length with actual_forecast
  stopifnot(length(forecast) == length(actual_forecast))
  
  # Dates corresponding to forecast period
  valid_indices <- which(forecast_dates >= ymd("2000-01-01"))
  forecast_filtered <- forecast[valid_indices]
  actual_filtered <- actual_forecast[valid_indices]
  
  # Calculate RMSFE
  rmsfe <- sqrt(mean((forecast_filtered - actual_filtered)^2, na.rm = TRUE))
  return(c(model = name, RMSFE = rmsfe))
}

# Compute RMSFE for all models
results <- rbind(
  forecasting_2000(Forecast_AR_1, "AR"),
  forecasting_2000(Forecast_MAR_1, "MAR"),
  forecasting_2000(Forecast_MAR_X_1, "MAR-X"),
  forecasting_2000(forecast_ART_1, "ART"),
  forecasting_2000(forecast_MART_pseudo_1, "MART"),
  forecasting_2000(forecast_MART_GS_1, "MART-GS"),
  forecasting_2000(forecast_MART_X_GS_1, "MART-X-GS"),
  forecasting_2000(forecast_MART_x_pseudo_1, "MART-X"),
  forecasting_2000(forecast_ART_x_1, "ART-X")
)

# Turn into data frame
rmsfe_df_1 <- as.data.frame(results)
rmsfe_df_1$RMSFE <- as.numeric(rmsfe_df_1$RMSFE)
print(rmsfe_df_1)


# Horizon of h=6
Forecast_AR_6 <- Forecast_AR[,6]
Forecast_MAR_X_6 <- Forecast_MAR[,6]

Forecast_MAR_6 <- Forecast_MAR_X[,6]

identical(Forecast_MAR_6, Forecast_MAR_6)


forecast_ART_6 <- forecast_art[,6]

#load("~forecast_x_results.RData")
#forecast_ART_X <- forecast_art_x[,1]
#forecast_MART_X_GS <- forecast_mart_x[,1]
#RMSE_ART_X <- 

forecast_MART_GS_6 <- forecast_mart_grid[,6]

forecast_MART_x_pseudo_6 <- forecast_mart_x_pseudo[,6]
forecast_MART_pseudo_6 <- forecast_mart_pseudo[,6]

forecast_MART_X_GS_6 <- forecast_mart_x_grid[,6]

forecast_ART_x_6 <- forecast_art_x[,6]

forecast_ART_6 <- forecast_art[,6]

actual <- actual_matrix[,6]
actual_forecast <- actual[forecast_start:n_months]

# Compute RMSFE for all models
results <- rbind(
  forecasting_2000(Forecast_AR_6, "AR"),
  forecasting_2000(Forecast_MAR_6, "MAR"),
  forecasting_2000(Forecast_MAR_X_6, "MAR-X"),
  forecasting_2000(forecast_ART_6, "ART"),
  forecasting_2000(forecast_MART_pseudo_6, "MART"),
  forecasting_2000(forecast_MART_GS_6, "MART-GS"),
  forecasting_2000(forecast_MART_X_GS_6, "MART-X-GS"),
  forecasting_2000(forecast_MART_x_pseudo_6, "MART-X"),
  forecasting_2000(forecast_ART_x_6, "ART-X")
)

# Turn into data frame
rmsfe_df_6 <- as.data.frame(results)
rmsfe_df_6$RMSFE <- as.numeric(rmsfe_df_6$RMSFE)
print(rmsfe_df_6)


# Horizon of h=12
Forecast_AR_12 <- Forecast_AR[,12]
Forecast_MAR_12 <- Forecast_MAR[,12]

Forecast_MAR_X_12 <- Forecast_MAR_X[,12]

identical(Forecast_MAR_12, Forecast_MAR_12)


forecast_ART_12 <- forecast_art[,12]

#load("~forecast_x_results.RData")
#forecast_ART_X <- forecast_art_x[,1]
#forecast_MART_X_GS <- forecast_mart_x[,1]
#RMSE_ART_X <- 

forecast_MART_GS_12 <- forecast_mart_grid[,12]

forecast_MART_x_pseudo_12 <- forecast_mart_x_pseudo[,12]
forecast_MART_pseudo_12 <- forecast_mart_pseudo[,12]

forecast_MART_X_GS_12 <- forecast_mart_x_grid[,12]

forecast_ART_x_12 <- forecast_art_x[,12]

forecast_ART_12 <- forecast_art[,12]

actual <- actual_matrix[,12]
actual_forecast <- actual[forecast_start:n_months]

# Compute RMSFE for all models
results <- rbind(
  forecasting_2000(Forecast_AR_12, "AR"),
  forecasting_2000(Forecast_MAR_12, "MAR"),
  forecasting_2000(Forecast_MAR_X_12, "MAR-X"),
  forecasting_2000(forecast_ART_12, "ART"),
  forecasting_2000(forecast_MART_pseudo_12, "MART"),
  forecasting_2000(forecast_MART_GS_12, "MART-GS"),
  forecasting_2000(forecast_MART_X_GS_12, "MART-X-GS"),
  forecasting_2000(forecast_MART_x_pseudo_12, "MART-X"),
  forecasting_2000(forecast_ART_x_12, "ART-X")
)

# Turn into data frame
rmsfe_df_12 <- as.data.frame(results)
rmsfe_df_12$RMSFE <- as.numeric(rmsfe_df_12$RMSFE)
print(rmsfe_df_12)

# Comparing best performing method with causal counterpart
compute_dm_tests <- function(forecast1, forecast2, actual, h) {
  p_values <- numeric(h)
  for (i in 1:h) {
    e1 <- actual - forecast1
    e2 <- actual - forecast2
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

DM_1 <- compute_dm_tests(forecast_MART_x_pseudo_1,forecast_ART_x_1,actual_forecast[,1],1)
DM_6 <- compute_dm_tests(forecast_MART_x_pseudo_6,forecast_ART_x_6,actual_forecast[,6],)
DM_12 <- compute_dm_tests(forecast_MART_x_pseudo_12,forecast_ART_x_12,actual_forecast[,12],1)
