# Load required library
library(MASS)
source("MARX_functions.R")
source("MART.R")

# Loop over your values and capture the printed output
p_C_max <- 6
p_NC_max <- 6

AIC <- matrix(NA, nrow = p_C_max, ncol = p_NC_max)
BIC <- matrix(NA, nrow = p_C_max, ncol = p_NC_max)
HQ <- matrix(NA, nrow = p_C_max, ncol = p_NC_max)

for (i in 0:p_C_max) {
  for (j in 0:p_NC_max) {
    marx_loop <- marx.t(inflation_df_monthly$inflationNonSA, NULL, p_C = i, p_NC = j)
    information <- information.criteria(type="MARX", marx_loop)
    AIC[i,j] <- information$aic
    BIC[i,j] <- information$bic
    HQ[i,j] <- information$hq
  }
}

