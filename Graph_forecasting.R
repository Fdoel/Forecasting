load("MART vs ART gerund door max.RData")
# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)

# Ik heb nodig:
# Forecasting AR, MAR ???
# Forecasting ART, MART en gs MART
# Forecasting AR-X, MAR-X
# Forecasting ART-X, MART-X en gs MART-X

# Pakken nu wel de goede vraag ik me af, ik denk aaleen als 

load("AR and MAR forecasting with many draws.RData")
Forecast_AR <-  forecast_causal
Forecast_MAR <- forecast_mixed
RMSE_AR <- rmse_causal
RMSE_MAR <- rmse_mixed

load("AR-X and MAR-X forecasting.RData")
Forecast_MAR_X <- forecast_causal
RMSE_MAR_X <- rmse_causal

identical(Forecast_MAR, Forecast_MAR_X)

load("MART vs ART gerund door max.RData")
forecast_ART <- forecast_art[,1]
RMSE_ART <- rmse_art

#load("~forecast_x_results.RData")
#forecast_ART_X <- forecast_art_x[,1]
#forecast_MART_X_GS <- forecast_mart_x[,1]
#RMSE_ART_X <- 

load("forecast_x_pseudo2_results.RData")
forecast_MART <- forecast_mart_pseudo[,1]
forecast_MART_GS <- forecast_mart_grid[,1]
RMSE_MART <- rmse_mart
RMSE_MART_GS <- rmse_

load("forecast_x_pseudo_results.RData")
forecast_MART_X <- forecast_mart_x_pseudo[,1]
RMSE_MART_X <- rmse_mart

# Set time index
start_year <- 1959
n_months <- 787 -12 
dates <- seq(ymd(paste0(start_year, "-06-01")), by = "month", length.out = n_months)

# Create actuals time series
actuals <- actual_matrix[,12]  # assuming first column contains actuals

# Time series for forecast start
forecast_start <- 250
forecast_dates <- dates[forecast_start:n_months]

# Create a data frame with all forecasts (aligned to same length)
df <- data.frame(
  Date = forecast_dates,
  Actual = actuals,
  #AR = forecast_AR,
  #MAR = forecast_MAR,
  ART = forecast_ART,
  MART_GS = forecast_MART_GS
  #`ART-X` = forecast_ART_X,
  #`MART-X` = forecast_MART_X_GS
)

# Convert to long format for ggplot
df_long <- df %>%
  pivot_longer(cols = -c(Date, Actual), names_to = "Model", values_to = "Forecast")

# Calculate relative error (or change if you prefer another metric)
df_long <- df_long %>%
  mutate(RelativeError = Forecast - Actual)

# Plot
ggplot(df_long, aes(x = Date, y = RelativeError, color = Model)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    title = "Forecast Error Relative to Actual Values",
    subtitle = "From month 250 onward (~1980)",
    y = "Forecast - Actual",
    x = "Date"
  ) +
  theme_minimal() +
  theme(legend.title = element_blank())

# Use MART_GS model residuals
residuals <- actuals - forecast_MART_GS

# Estimate standard deviation of residuals
resid_sd <- sd(residuals, na.rm = TRUE)

# Confidence Interval (for mean prediction) — approximate
ci_upper <- actuals + 1.96 * (resid_sd / sqrt(length(residuals)))
ci_lower <- actuals - 1.96 * (resid_sd / sqrt(length(residuals)))

# Prediction Interval (for individual future obs)
pi_upper <- actuals + 1.96 * resid_sd
pi_lower <- actuals - 1.96 * resid_sd

df_base <- data.frame(
  Date = forecast_dates,
  Actual = actuals,
  CI_Lower = ci_lower,
  CI_Upper = ci_upper,
  PI_Lower = pi_lower,
  PI_Upper = pi_upper
)

df_forecasts <- df %>%
  select(-Actual) %>%
  pivot_longer(cols = -Date, names_to = "Model", values_to = "Forecast")

ggplot() +
  # Prediction interval (wider band)
  geom_ribbon(data = df_base, aes(x = Date, ymin = PI_Lower, ymax = PI_Upper),
              fill = "lightblue", alpha = 0.3, show.legend = FALSE) +
  
  # Confidence interval (narrower band)
  geom_ribbon(data = df_base, aes(x = Date, ymin = CI_Lower, ymax = CI_Upper),
              fill = "blue", alpha = 0.2, show.legend = FALSE) +
  
  # Actual values
  geom_line(data = df_base, aes(x = Date, y = Actual), color = "black", size = 1) +
  
  # Forecasts
  geom_line(data = df_forecasts, aes(x = Date, y = Forecast, color = Model), size = 1) +
  
  labs(
    title = "Forecasts with 95% Confidence and Prediction Intervals",
    x = "Date",
    y = "Value"
  ) +
  theme_minimal() +
  theme(legend.title = element_blank())