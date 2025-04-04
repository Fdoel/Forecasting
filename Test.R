#test

# Load required library
library(MASS)
source("MARX_functions.R")
source("MART.R")

test_data <- inflation_df

# Loop over your values and capture the printed output
p_C_max <- 5
p_NC_max <- 5

AIC <- matrix(NA, nrow = p_C_max, ncol = p_NC_max)
BIC <- matrix(NA, nrow = p_C_max, ncol = p_NC_max)
HQ <- matrix(NA, nrow = p_C_max, ncol = p_NC_max)

for (i in 0:p_C_max) {
  
  for (j in 0:p_NC_max) {
    marx_loop <- marx.t(test_data$inflation, NULL, p_C = i, p_NC = j)
    
    information <- information_criteria_MARX(marx_loop, test_data$inflation, NULL)
    AIC[i,j] <- information[1]
    BIC[i,j] <- information[2]
    HQ[i,j] <- information[3]
  }
}

marx(inflation_df$inflation, NULL, p_max=10, sig_level=0.1)