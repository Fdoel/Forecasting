# Code voor simulatie van model.
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
  
  simulation <- as.data.frame(simulation)
  colnames(simulation) <- c("Causal", "nonCausal", "models")
  cols_to_use <- simulation$models
  cols_to_use <- as.data.frame(simulation$models)
  colnames(cols_to_use) <- seq(1, ncol(cols_to_use))
  
  # Column names
  columns <- c("int", "lag 1", "lead 1", "exo", "df", "scale")
  
  # Creating DataFrame
  df = as.data.frame(cols_to_use, columns=columns)
  
  # Calculate row means
  means <- apply(df, 1, mean)
  sds <- apply(df, 1, sd)
  
  
  mean_estimates <- c(means[2] ,means[3])
  sd_estimates <- c(sds[2], sds[3])
  df_estimates <- means[5]
  scale_estimates <- means[6]
  
  return(list(
    mean = mean_estimates,
    sd = sd_estimates,
    df = df_estimates,
    scale = scale_estimates
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
simulation_estimates <- matrix(NA, nrow = length(sample_sizes) * length(param_combinations), ncol = 9)
colnames(simulation_estimates) <- c("Sample_size", "phi", "psi", "lag_est", "lag_sd", "lead_est", "lead_sd", "df", "scale")
i=1


sim_frama <- monte_carlo_simulation(300, 0.9, 0.9, 100)
get_final_estimates(sim_frama)

for (T in sample_sizes) {
  for (params in param_combinations) {
    phi <- params[1]
    psi <- params[2]
    
    simulation <- monte_carlo_simulation(T, phi, psi, n_sim)
    estimates <- get_final_estimates(simulation)
    
    simulation_estimates[i,1] = T
    simulation_estimates[i,2] = phi
    simulation_estimates[i,3] = psi
    simulation_estimates[i, 4] = estimates$mean[1]  # Lag estimate
    simulation_estimates[i, 5] = estimates$sd[1]  # Lag standard deviation
    simulation_estimates[i, 6] = estimates$mean[2]  # Lead estimate
    simulation_estimates[i, 7] = estimates$sd[2]  # Lead standard deviation
    simulation_estimates[i, 8] = estimates$df  # df value
    simulation_estimates[i, 9] = estimates$scale  # scale value
    
    i <- i+1
  }
}