# Load required library

library(MASS)
source("MARX_functions.R")
source("MART.R")

load("inflation_df_monthly.RData")


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
    marx_loop <- MART(inflation_df_monthly$inflationNonSA, NULL, p_C = i, p_NC = j, 1)
    information <- information.criteria(type="MART", marx_loop)
    LL[(i+1),(j+1)] <- information$loglikelihood
    AIC[(i+1),(j+1)] <- information$aic
    BIC[(i+1),(j+1)] <- information$bic
    HQ[(i+1),(j+1)] <- information$hq
    N[(i+1),(j+1)] <- information$n
  }
  
}