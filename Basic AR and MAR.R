# Load required library
library(MASS)
source("MARX_functions.R")
source("MART.R")

# Loop over your values and capture the printed output
p_C_max <- 12
p_NC_max <- 12

LL <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)
AIC <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)
BIC <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)
HQ <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)
N <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)

for (i in 0:p_C_max) {
  for (j in 0:p_NC_max) {
    marx_loop <- marx.t(inflation_df_monthly$inflationNonSA, NULL, p_C = i, p_NC = j)
    information <- information.criteria(type="MARX", marx_loop)
    LL[(i+1),(j+1)] <- information$loglikelihood
    AIC[(i+1),(j+1)] <- information$aic
    BIC[(i+1),(j+1)] <- information$bic
    HQ[(i+1),(j+1)] <- information$hq
    N[(i+1),(j+1)] <- information$n
  }
}





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







