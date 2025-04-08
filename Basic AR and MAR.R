# Load required library
library(MASS)
source("MARX_functions.R")
source("MART.R")
library(forecast)

# Loop over your values and capture the printed output
if(FALSE) {
  marx(inflation_df_monthly$inflationNonSA, NULL, p_max = 18, sig_level = 0.1)
  
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
}


#perfom iid test of resids
model_ar12 <- Arima(inflation_df_monthly$inflationNonSA, order = c(12, 0, 0))
resids_ar12 <- model_ar12$residuals

# Step 2: Create the squared residuals
resids_sq <- resids_ar12^2

# Step 3: Regress residuals on lagged squared residuals (you can set m to desired number of lags)
m <- 12
X <- embed(resids_sq, m + 1)
y <- X[, 1]
X_lags <- X[, -1]

# Step 4: Run the regression
model_test <- lm(y ~ X_lags)

# Step 5: Perform the joint significance test (H0: all δ's = 0)
test_statistic <- summary(model_test)$r.squared * length(y)
p_value <- pchisq(test_statistic, df = m, lower.tail = FALSE)

# Step 6: Output
cat("Chi-squared test statistic:", test_statistic, "\n")
cat("p-value:", p_value, "\n")

ar12 <- marx.t(inflation_df_monthly$inflationNonSA, NULL, p_C = 12, p_NC = 0)
ar12_crit <- information.criteria("MARX", ar12)

mar_final <- marx.t(inflation_df_monthly$inflationNonSA, NULL, p_C = 1, p_NC = 11)
mar111_crit <- information.criteria("MARX", mar_final)


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







