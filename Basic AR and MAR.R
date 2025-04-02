# Load required library
library(MASS)
source("MARX_functions.R")

# Create an empty list to store outputs
results_list <- list()

# Loop over your values and capture the printed output
for (i in 0:6) {
  for (j in 0:6) {
    output <- capture.output(
      marx_loop <- marx(inflation_df_monthly$inflationSA, NULL, p_max = 20, sig_level = 0.1, p_C = i, p_NC = j)
    )
    # Save the captured output with a meaningful name
    results_list[[paste0("p_C_", i, "_p_NC_", j)]] <- output
  }
}

#example of marx function, simple MAR model
marx_example <- marx(inflation_df_monthly$inflationSA, NULL, p_max = 20, sig_level = 0.1) #user friendly function with interface
