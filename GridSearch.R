# Install pbmcapply if not already installed
# install.packages("pbmcapply")

library(pbmcapply)

# Load your scripts and data
source("MART.R")
load("inflation_df_monthly.RData")

# Set model parameters
thresholds <- seq(0, 1, by = 0.1)
gamma <- seq(1, 10, by=1)
p_C_max <- 12
p_NC_max <- 12

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
  gamma = gammas,
  i = 0:p_C_max,
  j = 0:p_NC_max
)

# Function to run MART and get BIC value
run_model <- function(params) {
  t <- params$threshold
  g <- params$gamma
  i <- params$i
  j <- params$j
  
  # Run MART model
  SMART_c <- SMART(inflation_df_monthly$inflationNonSA, NULL, i, j, t, g)
  bic_value <- information.criteria("SMART", SMART_c)
  return(data.frame(threshold = t, gamma = g, i = i, j = j, bic = bic_value))
}

# Use pbmclapply with a progress bar
bic_results <- pbmclapply(
  1:nrow(param_grid),
  function(idx) run_model(param_grid[idx, ]),
  mc.cores = n_cores
)

# Combine all results into one data frame
bic_smart_df <- do.call(rbind, bic_results)

# Save the results
save(bic_smart_df, file = "bic_smart_df.RData")