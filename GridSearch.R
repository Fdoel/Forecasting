# Install pbmcapply if not already installed
# install.packages("pbmcapply")

library(pbmcapply)

# Load your scripts and data
source("MART.R")
load("inflation_df_monthly.RData")

# Set model parameters
thresholds <- median(inflation_df_monthly$inflationNonSA, na.rm = TRUE)
ds <- seq(1, 12)
p_C_max <- 4
p_NC_max <- 4

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
  i = 0:p_C_max,
  j = 0:p_NC_max
)

# Function to run MART and get BIC value
run_model <- function(params) {
  t <- params$threshold
  d <- params$d
  i <- params$i
  j <- params$j
  
  # Run MART model
  MART_d <- MART(inflation_df_monthly$inflationNonSA, NULL, i, j, t, d)
  bic_value <- information.criteria("MART", MART_d)
  return(data.frame(threshold = t, d=d, i = i, j = j, bic = bic_value))
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
save(bic_mart_d_df, file = "bic_mart_d_df.RData")
