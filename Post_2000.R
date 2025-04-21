library(lubridate)


load("~/Documents/GitHub/Forecasting/Forecasting results/AR and MAR forecasting with many draws.RData")
Forecast_AR <-  forecast_causal
Forecast_MAR <- forecast_mixed
Forecast_AR <- Forecast_AR[,1]
Forecast_MAR <- Forecast_MAR[,1]

load("~/Documents/GitHub/Forecasting/Forecasting results/AR-X and MAR-X forecasting.RData")
Forecast_MAR_X <- forecast_causal
Forecast_MAR_X <- Forecast_MAR_X[,1]

identical(Forecast_MAR, Forecast_MAR_X)

load("~/Documents/GitHub/Forecasting/Forecasting results/forecast_ART_results.RData")
forecast_ART <- forecast_art[,1]

#load("~forecast_x_results.RData")
#forecast_ART_X <- forecast_art_x[,1]
#forecast_MART_X_GS <- forecast_mart_x[,1]
#RMSE_ART_X <- 

load("~/Documents/GitHub/Forecasting/Forecasting results/forecast_MART_results.RData")
forecast_MART_GS <- forecast_mart_grid[,1]

load("~/Documents/GitHub/Forecasting/Forecasting results/forecast_x_pseudo_results.RData")
forecast_MART_x_pseudo <- forecast_mart_x_pseudo[,1]
forecast_MART_pseudo <- forecast_mart_pseudo[,1]

load("~/Documents/GitHub/Forecasting/Forecasting results/forecast_MARTX_results.RData")
forecast_MART_X_GS <- forecast_mart_x_grid[,1]

load("~/Documents/GitHub/Forecasting/Forecasting results/forecast_ARTX_results.RData")
forecast_ART_x <- forecast_art_x[,1]

load("~/Documents/GitHub/Forecasting/Forecasting results/forecast_ART_results.RData")
forecast_ART <- forecast_art[,1]


# Setup the actuals
start_year <- 1959
n_months <- 775  # as you defined earlier
dates <- seq(ymd(paste0(start_year, "-06-01")), by = "month", length.out = n_months)
forecast_start <- 250
forecast_dates <- dates[forecast_start:n_months]
actual <- actual_matrix
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
  forecasting_2000(Forecast_AR, "AR"),
  forecasting_2000(Forecast_MAR, "MAR"),
  forecasting_2000(Forecast_MAR_X, "MAR-X"),
  forecasting_2000(forecast_ART, "ART"),
  forecasting_2000(forecast_MART_pseudo, "MART"),
  forecasting_2000(forecast_MART_GS, "MART-GS"),
  forecasting_2000(forecast_MART_X_GS, "MART-X-GS"),
  forecasting_2000(forecast_MART_x_pseudo, "MART-X"),
  forecasting_2000(forecast_ART_x, "ART-X")
)

# Turn into data frame
rmsfe_df <- as.data.frame(results)
rmsfe_df$RMSFE <- as.numeric(rmsfe_df$RMSFE)
print(rmsfe_df)