library(MASS)
library(parallel)
install.packages("pbmcapply")
library(pbmcapply)
source("MARX_functions.R")
source("MART.R")

load("inflation_df_monthly.RData")

p_C_max <- 12
p_NC_max <- 12

if (.Platform$OS.type == "windows") {
  n_cores <- 1
} else {
  n_cores <- 2
  RNGkind("L'Ecuyer-CMRG")
}

threshold <- 0.29
gamma <- 20

param_grid <- expand.grid(i = 0:p_C_max, j = 0:p_NC_max)

run_model_info <- function(params) {
  i <- params$i
  j <- params$j
  
  model <- SMART(inflation_df_monthly$inflationNonSA, NULL, p_C = i, p_NC = j, gamma, threshold)
  info <- information.criteria(type = "SMART", model)
  
  return(data.frame(i = i, j = j,
                    loglikelihood = info$loglikelihood,
                    aic = info$aic,
                    bic = info$bic,
                    hq = info$hq,
                    n = info$n))
}

info_results <- pbmclapply(
  1:nrow(param_grid),
  function(idx) run_model_info(param_grid[idx, ]),
  mc.cores = n_cores
)

# Combine results into a data frame
info_df <- do.call(rbind, info_results)

# Now convert to matrices
LL <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)
AIC <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)
BIC <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)
HQ <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)
N <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)

for (row in 1:nrow(info_df)) {
  i <- info_df$i[row] + 1
  j <- info_df$j[row] + 1
  LL[i, j] <- info_df$loglikelihood[row]
  AIC[i, j] <- info_df$aic[row]
  BIC[i, j] <- info_df$bic[row]
  HQ[i, j] <- info_df$hq[row]
  N[i, j] <- info_df$n[row]
}
