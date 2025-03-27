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
n_sim <- 10000
sample_sizes <- c(100, 200, 500)
param_combinations <- list(
  c(0.9, 0.9),
  c(0.9, 0.1),
  c(0.1, 0.9)
)

#hier ff mooie loop maken met sample sizes and param combinations

sim_500_9_9 <- monte_carlo_simulation(500,0.9,0.9,1000)

est <- get_final_estimates(sim_500_9_9)


