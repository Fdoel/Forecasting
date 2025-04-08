# Install pbmcapply if not already installed
# install.packages("pbmcapply")

library(pbmcapply)

# Load your scripts and data
source("MART.R")
load("inflation_df_monthly.RData")

# Set model parameters
thresholds <- seq(0, 1, by = 0.1)
p_C_max <- 12
p_NC_max <- 12

# Set number of cores based on OS
if (.Platform$OS.type == "windows") {
  n_cores <- 1
} else {
  n_cores <- 2
  RNGkind("L'Ecuyer-CMRG")  # Safe parallel RNG
}

# Create a parameter grid of all combinations
param_grid <- expand.grid(
  threshold = thresholds,
  i = 0:p_C_max,
  j = 0:p_NC_max
)

# Function to run MART and get BIC value
run_model <- function(params) {
  t <- params$threshold
  i <- params$i
  j <- params$j
  
  # Run MART model
  MART_c <- MART(inflation_df_monthly$inflationNonSA, NULL, i, j, t)
  bic_value <- information.criteria("MART", MART_c)
  
  # Return result as a row
  return(data.frame(threshold = t, i = i, j = j, bic = bic_value))
}

# Use pbmclapply with a progress bar
bic_results <- pbmclapply(
  1:nrow(param_grid),
  function(idx) run_model(param_grid[idx, ]),
  mc.cores = n_cores
)

# Combine all results into one data frame
bic_c_df <- do.call(rbind, bic_results)

# Find the best parameter combination (lowest BIC)
min_bic_row <- bic_c_df[which.min(bic_c_df$bic), ]

# Print the best result
print(min_bic_row)
