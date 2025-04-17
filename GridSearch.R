library(pbmcapply)
source("MART.R")
load("inflation_df_monthly.RData")

# Set model parameters
thresholds <- seq(0.1, 0.6, by = 0.1) 
ds <- seq(1, 6)
gammas <- seq(1,15)
p_C_max <- 6
p_NC_max <- 0

# Set number of cores based on OS
if (.Platform$OS.type == "windows") {
  n_cores <- 1
} else {
  n_cores <- 3
  RNGkind("L'Ecuyer-CMRG")  # Safe parallel RNG
}

# Create a parameter grid of all combinations
param_grid <- expand.grid(
  threshold = thresholds,
  d = ds,
  g = gammas,
  i = 1:p_C_max, 
  j = 0
)

# Function to run MART and get BIC value
run_model <- function(params) {
  t <- params$threshold
  d <- params$d
  gamma <- params$g
  i <- params$i
  j <- params$j
  
  # Run MART model
  SMART_d <- SMART(inflation_df_monthly$inflationNonSA, NULL, i, j, t, gamma, d)
  bic_value <- information.criteria("SMART", SMART_d)
  return(data.frame(threshold = t, d = d, i = i, j = j, bic = bic_value))
}

# Use pbmclapply with a progress bar
bic_results <- pbmclapply(
  1:nrow(param_grid),
  function(idx) run_model(param_grid[idx, ]),
  mc.cores = n_cores
)

# Combine all results into one data frame
bic_mart_d_df <- do.call(rbind, bic_results)

# Save the results
save(bic_sart_tdig_df, file = "bic_sart_t_d_i_g.RData")

View(bic_sart_tdig_df)


