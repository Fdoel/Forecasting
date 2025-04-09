# Load required library
library(MASS)
source("MARX_functions.R")
source("MART.R")
library(forecast)
library(pbmcapply)

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
p <- 12
M <- 50
forecasts_size <- length(inflation_df_monthly[,1]) - p
forecast_vector <- vector(mode = "numeric", length = forecasts_size)

for(i in 1:(forecasts_size)) {
  forecast_data <- inflation_df_monthly[i:forecasts_size+i, ]
  forecast_vector[i] <- forecast.marx(y=forecast_data$inflationNonSA, p_C=1, p_NC=11, h=1, M=M, N=1000)
}

library(pbmcapply)
library(pbmcapply)

# Parameters
h <- 12
N <- 1000
M <- 50

# Model specifications
p_C_mixed <- 1;  p_NC_mixed <- 11   # Mixed model
p_C_causal <- 12; p_NC_causal <- 0  # Purely causal model
p_C_mid <- 11; p_NC_mid <- 1        # 11-lag, 1-lead model (new)

# Data
data_series <- inflation_df_monthly$inflationNonSA
start_index <- 100
end_index <- length(data_series) - h
forecast_indices <- start_index:end_index

# Parallel forecast with progress bar
results_list <- pbmclapply(
  X = forecast_indices,
  FUN = function(t) {
    y_window <- data_series[1:t]
    
    forecast_mixed <- forecast.marx(
      y = y_window,
      p_C = p_C_mixed,
      p_NC = p_NC_mixed,
      h = h,
      M = M,
      N = N
    )
    
    forecast_causal <- forecast.marx(
      y = y_window,
      p_C = p_C_causal,
      p_NC = p_NC_causal,
      h = h,
      M = M,
      N = N
    )
    
    forecast_mid <- forecast.marx(
      y = y_window,
      p_C = p_C_mid,
      p_NC = p_NC_mid,
      h = h,
      M = M,
      N = N
    )
    
    actual <- data_series[(t + 1):(t + h)]
    
    return(list(
      mixed = forecast_mixed,
      causal = forecast_causal,
      mid = forecast_mid,
      actual = actual
    ))
  },
  mc.cores = parallel::detectCores() - 1
)

# Convert list results to matrices
forecast_mixed <- do.call(rbind, lapply(results_list, `[[`, "mixed"))
forecast_causal <- do.call(rbind, lapply(results_list, `[[`, "causal"))
forecast_mid <- do.call(rbind, lapply(results_list, `[[`, "mid"))
actual_matrix <- do.call(rbind, lapply(results_list, `[[`, "actual"))

# Label columns
colnames(forecast_mixed) <- paste0("h", 1:h)
colnames(forecast_causal) <- paste0("h", 1:h)
colnames(forecast_mid) <- paste0("h", 1:h)
colnames(actual_matrix) <- paste0("h", 1:h)

# === RMSE computation ===
rmse <- function(forecast, actual) {
  sqrt(colMeans((forecast - actual)^2, na.rm = TRUE))
}

rmse_mixed <- rmse(forecast_mixed, actual_matrix)
rmse_causal <- rmse(forecast_causal, actual_matrix)
rmse_mid <- rmse(forecast_mid, actual_matrix)

rmse_df <- data.frame(
  horizon = 1:h,
  RMSE_mixed = rmse_mixed,
  RMSE_causal = rmse_causal,
  RMSE_mid = rmse_mid
)

# Print RMSE comparison
print(rmse_df)