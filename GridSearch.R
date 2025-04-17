library(pbmcapply)

# Load your scripts and data
source("MART.R")
load("inflation_df_monthly.RData")

# Set model parameters
thresholds <- seq(0.1, 0.6, by = 0.1) #0.6 # median(inflation_df_monthly$inflationNonSA, na.rm = TRUE)
ds <- seq(1, 6)
gammas <- seq(1,10)
p_C_max <- 6
p_NC_max <- 2

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
  x = NULL,
  y = inflation_df_monthly$inflationNonSA,
  d = ds,
  i = 0:p_C_max, 
  j = 0
)

# Function to run MART and get BIC value
run_model <- function(params) {
  t <- params$threshold
  y <- params$y
  x <- params$x
  d <- params$d
  i <- params$i
  j <- params$j
  
  # Run MART model
  MART_d <- MART(y, x, i, j, t, d)
  bic_value <- information.criteria("MART", MART_d)
  return(data.frame(threshold = t, d = d, i = i, j = j, bic = bic_value))
}

# Use pbmclapply with a progress bar
bic_results <- pbmclapply(
  1:nrow(param_grid),
  function(idx) run_model(param_grid[idx, ]),
  mc.cores = n_cores
)

# Combine all results into one data frame
bic_art_x_t_d_df <- do.call(rbind, bic_results)

# Save the results
save(bic_art_x_t_d_df, file = "bic_art_x_t_d_df.RData")

View(bic_art_x_t_d_df)



