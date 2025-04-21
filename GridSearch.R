library(pbmcapply)
source("MART.R")
load("inflation_df_monthly.RData")

# Set model parameters
thresholds <- seq(0.1, 0.6, by = 0.1) 
ds <- seq(1, 6)
gammas <- seq(1,10)
p_C_max <- 6
p_NC_max <- 6

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
  #g = gammas,
  i = 0:p_C_max, 
  j = 0:p_NC_max
)

# Function to run MART and get BIC value
run_model <- function(params) {
  t <- params$threshold
  d <- params$d
  #gamma <- params$g
  i <- params$i
  j <- params$j
  
  # Run MART model
  MART_d <- MART(inflation_df_monthly$inflationSA, NULL, i, j,t,d)
  bic_value <- information.criteria("MART", MART_d)
  return(data.frame(threshold=t,d=d,i = i, j = j, bic = bic_value))
}

# Use pbmclapply with a progress bar
bic_results <- pbmclapply(
  1:nrow(param_grid),
  function(idx) run_model(param_grid[idx, ]),
  mc.cores = n_cores
)

# Combine all results into one data frame
ART_SA <- do.call(rbind, bic_results)

# Voeg de gamma-kolom toe aan het uiteindelijke data frame
#bic_smart_tdig_df$gamma <- param_grid$g


# Save the results
save(ART_SA, file = "ART_SA.RData")

View(ART_SA)


