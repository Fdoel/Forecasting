# Load required library
library(MASS)
source("MARX_functions.R")
source("MART.R")

# Loop over your values and capture the printed output
p_C_max <- 12
p_NC_max <- 12

AIC <- matrix(NA, nrow = p_C_max, ncol = p_NC_max)
BIC <- matrix(NA, nrow = p_C_max, ncol = p_NC_max)
HQ <- matrix(NA, nrow = p_C_max, ncol = p_NC_max)

for (i in 0:p_C_max) {
  
  for (j in 0:p_NC_max) {
    marx_loop <- marx.t(inflation_df_monthly$inflationSA, NULL, p_C = i, p_NC = j)
    
    information <- information_criteria_MARX(marx_loop, inflation_df_monthly$inflationSA, NULL)
    AIC[i,j] <- information[1]
    BIC[i,j] <- information[2]
    HQ[i,j] <- information[3]
  }
}

marx_v2 <- marx.t(inflation_df_monthly$inflationSA, NULL, p_C=3, p_NC=3)

criteria <- information_criteria_MARX(marx_v2, inflation_df_monthly$inflationSA, NULL)

fit <- ar(inflation_df_monthly$inflationSA, order.max=20)


information_criteria_MARX <- function(marx_model, y, x) {
  mse <- sum(marx_model$residuals^2)
  n <- length(marx_model$residuals)
  
  if (length(x) > 1){
    numcol <- NCOL(x)
  }
  else{
    numcol = 0
  }
  
  p <- length(marx_model$coef.c) + length(marx_model$coef.nc)
  
  aic <- -2*mse/n + (2/n)*(p+1+numcol)
  bic <- -2*mse/n + ((log(n))/n)*(p+1+numcol)
  hq <- -2*mse/n + ((2*log(log(n)))/n)*(p+1+numcol)
  return(c(aic = aic, bic = bic, hq = hq))
}





marx(inflation_df_monthly$inflationSA, NULL, p_max=6, sig_level=0.1)
marx(inflation_df_monthly$inflationNonSA, NULL, p_max=6, sig_level=0.1)

yoy_data <- inflation_df_monthly[-(1:12), ]
marx(yoy_data$inflation_yoy_nonSA, NULL, p_max=6, sig_level=0.1)










