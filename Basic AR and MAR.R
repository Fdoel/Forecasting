# Load required library
library(MASS)
source("MARX_functions.R")


#gwn interface volgen
marx(inflation_df_monthly$inflationNonSA, NULL, p_max=6, sig_level = 0.1)

marx_test1 <- mixed(inflation_df_monthly$inflationNonSA, NULL, p_C=0, p_NC=4)
marx_test2 <- mixed(inflation_df_monthly$inflationNonSA, NULL, p_C=1, p_NC=3)
marx_test3 <- mixed(inflation_df_monthly$inflationNonSA, NULL, p_C=2, p_NC=2)
marx_test4 <- mixed(inflation_df_monthly$inflationNonSA, NULL, p_C=3, p_NC=1)
marx_test5 <- mixed(inflation_df_monthly$inflationNonSA, NULL, p_C=4, p_NC=0)

#forecasting test
p <- 4
M <- 50
buffer <- 50
forecasts_size <- length(inflation_df_monthly[,1]) - p - M - buffer
forecast_vector <- vector(mode = "numeric", length = forecasts_size)

for(i in 1:(forecasts_size-buffer)) {
  forecast_data <- inflation_df_monthly[i:forecasts_size, ]
  forecast_vector[i] <- forecast.marx(y=forecast_data$inflationNonSA, p_C=1, p_NC=3, h=1, M=M, N=100)
}

compare_data <- inflation_df_monthly[]







