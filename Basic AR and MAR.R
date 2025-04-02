# Load required library
library(MASS)
source("MARX_functions.R")

#loop and store results 
for (i in 0:6) {
  for (j in 0:6) {
    marx_loop <- marx(inflation_df_monthly$inflationSA, NULL, p_max = 20, sig_level = 0.1, p_C=i, p_NC=j)
  }
}

#example of marx function, simple MAR model
marx_example <- marx(inflation_df_monthly$inflationSA, NULL, p_max = 20, sig_level = 0.1) #user friendly function with interface
