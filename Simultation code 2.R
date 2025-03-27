set.seed(321) # Set seed for reproducibility

# Load required library
library(MASS)
source("MARX_functions.R")

# STEP 1: series from causal AR(p) model is generated.
ar_process <- function(phi, T) {
  epsilon_t <- rt(T + r, df = 3) * sigma  # t-distributed errors with df=3 en keer sigma?
  v <- rep(0, T + r)  
  for (t in (r + 1):T) {
    v[t] <- phi * v[t - 1] + epsilon_t[t]
  }
  return(v)
}

# STEP 2: series from non causal component is generated
nc_process <- function(psi, v, T) {
  y <- rep(0, T + s)  
  for (t in (T - s):1) {
    y[t] <- psi * y[t + 1] + v[t]
  }

  return(y)
}

# estimation mean
estimate_parameters <- function(y,v) {
  return(c(mean(y), mean(v)))
}
  
# Monte Carlo simulation
monte_carlo_simulation <- function(T, phi, psi, num_simulations = 10000) {
  estimates <- matrix(NA, nrow = num_simulations, ncol = 2)  #Storage for estimates
  models <- rep(NULL,num_simulations)
  
  for (sim in 1:num_simulations) {
    v_t <- ar_process(phi, T)
    y_t <- nc_process(psi, v_t, T)
    
    # Discard initial 100 observations at the beginning and end
    v_t <- v_t[(101):(T - 100)]
    y_t <- y_t[(101):(T - 100)]
    
    estimates[sim, ] <- estimate_parameters(y_t, v_t)
    
    models[sim] <- mixed(y_t,NULL,1,1)
  }
  
  colnames(estimates) <- c("Causal", "nonCausal")
  cbind(estimates, models)
}

#function for getting results from simulations
get_final_estimates <- function(simulation) {
  cols_to_use <- as.data.frame(simulation[, c("Causal", "nonCausal")])
  cols_to_use <- as.data.frame(lapply(cols_to_use, as.numeric))  # force numeric
  mean_estimates <- colMeans(cols_to_use)
  sd_estimates <- sd_estimates <- apply(cols_to_use, 2, sd)
  
  return(data.frame(
    mean = mean_estimates,
    sd = sd_estimates
  ))
}

#Estimate parameters
# Lag and lead 1,1
r <- 1
s <- 1
df <- 3
sigma <- 0.1

# Simulation parameters
n_sim <- 100
sample_sizes <- c(300, 500, 800)
param_combinations <- list(
  c(0.9, 0.9),
  c(0.9, 0.1),
  c(0.1, 0.9)
)

#hier ff mooie loop maken met sample sizes and param combinations
simulation_estimates <- matrix(NA, nrow = length(sample_sizes) * length(param_combinations), ncol = 7)
colnames(simulation_estimates) <- c("Sample_size", "phi", "psi", "lag_est", "lag_sd", "lead_est", "lead_sd")
i=1

for (n in sample_sizes) {
  for (params in param_combinations) {
    phi <- params[1]
    psi <- params[2]
    
    simulation <- monte_carlo_simulation(n, phi, psi, n_sim)
    estimates <- get_final_estimates(simulation)
    
    print(estimates)
    
    simulation_estimates[i,1] = n
    simulation_estimates[i,2] = phi
    simulation_estimates[i,3] = psi
    simulation_estimates[i,4] = estimates[1,1]
    simulation_estimates[i,5] = estimates[1,2]
    simulation_estimates[i,6] = estimates[2,1]
    simulation_estimates[i,7] = estimates[2,2]
    i <- i+1
  }
}






