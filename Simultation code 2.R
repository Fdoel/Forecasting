set.seed(321) # Set seed for reproducibility

# Load required library
library(MASS)

# Simulation parameters
n_sim <- 10000
sample_sizes <- c(100, 200, 500)
param_combinations <- list(
  c(0.9, 0.9),
  c(0.9, 0.1),
  c(0.1, 0.9)
)

# Lag and lead 1,1
r <- 1
s <- 1

# Hoe verwerk ik r and s initial zero???

df <- 3
sigma <- 0.1

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

# estimation mean and sd
estimate_parameters <- function(y,v) {
  return(c(mean(y), mean(v)))
}
  
# Monte Carlo simulation
monte_carlo_simulation <- function(T, phi, psi, num_simulations = 10000) {
  estimates <- matrix(NA, nrow = num_simulations, ncol = 2)  # Opslag voor schattingen
  
  for (sim in 1:num_simulations) {
    v_t <- ar_process(phi, T)
    y_t <- nc_process(psi, v_t, T)
    
    # Discard initial 100 observations at the beginning and end
    v_t <- v_t[(101):(T - 100)]
    y_t <- y_t[(101):(T - 100)]
    
    estimates[sim, ] <- estimate_parameters(y_t, v_t)
  }
  print((estimates))
}

# Er gaat iets mis, mean very very small
sim_500_9_9 <- monte_carlo_simulation(500,0.9,0.9,10000)

# Gemiddelde en standaarddeviatie van de schattingen
mean_estimates <- colMeans(sim_500_9_9)
sd_estimates <- apply(sim_500_9_9, 2, sd)

# Resultaten printen
print(paste("Gemiddelde van y:", mean_estimates[1], " Standaarddeviatie:", sd_estimates[1]))
print(paste("Gemiddelde van v:", mean_estimates[2], " Standaarddeviatie:", sd_estimates[2]))

