# Code voor simulatie van model.
set.seed(64) # Set seed for reproducibility

# Load required library
library(MASS)
source("MARX_functions.R")
source("MART.R")


# STEP 1: series from causal AR(p) model is generated.
ar_process <- function(phi_1, phi_2, T, c = 0.25) {
  d <- 1
  epsilon_t <- rt(T + r, df = 3) * sigma  # t-distributed errors with df=3 en keer sigma
  v <- rep(0, T + r - 1)  
  for (t in (r + 1):T) {
    c <- as.numeric(c)
    # Logistic smooth transition function
    smooth_transition <- 1 / (1 + exp(-gamma * (v[t - d] - theta)))
    
    # Weighted sum of phi_1 and phi_2 based on the logistic function
    phi_t <- (1 - smooth_transition) * phi_1 + smooth_transition * phi_2
    
    # Update the value of v[t] based on the smooth transition of coefficients
    v[t] <- phi_t * v[t - 1] + epsilon_t[t]
  }
  return(v)
}

# STEP 2: series from non causal component is generated
nc_process <- function(psi_1, psi_2, v, T, c = 0.25) {
  d <- 1
  y <- rep(0, T + s - 1)  
  for (t in (T - s):1) {
    # Ik vergelijk nu threshold c met gegenereerde data voor v, moet dat bij y? Of moet de threshold c met iets anders vergeleken worden?
    if (v[t + d] > c) { # plus d en niet min d
      y[t] <- psi_1 * y[t + 1] + v[t]
    }
    else{
      y[t] <- psi_2 * y[t + 1] + v[t]
    }
  }
  return(y)
}

# estimation mean
estimate_parameters <- function(y,v) {
  return(c(mean(y), mean(v)))
}

# Monte Carlo simulation
monte_carlo_simulation <- function(T, phi_1, phi_2, psi_1, psi_2, num_simulations = 10000, c = 0.25) {  
  estimates <- matrix(NA, nrow = num_simulations, ncol = 2)  # Storage for estimates  
  models <- vector("list", num_simulations)  # Initialize models as a list  
  
  for (sim in 1:num_simulations) {
    #message(paste("Simulation:", sim))  # Track simulation number  
    
    # Step 1: Generate AR process  
    v_t <- tryCatch({  
      ar_process(phi_1, phi_2, T, c)  
    }, error = function(e) {  
      stop(paste("Error in ar_process at simulation", sim, ":", e$message))  
    })  
    
    # Step 2: Generate NC process  
    y_t <- tryCatch({  
      nc_process(psi_1, psi_2, v_t, T, c)  
    }, error = function(e) {  
      stop(paste("Error in nc_process at simulation", sim, ":", e$message))  
    })  
    
    # Step 3: Trim observations  
    if (length(v_t) < (T - 200) | length(y_t) < (T - 200)) {  
      stop(paste("Error: Process length too short at simulation", sim))  
    }  
    v_t <- v_t[(101):(T - 100)]  
    y_t <- y_t[(101):(T - 100)]  
    
    # Step 4: Estimate parameters  
    estimates[sim, ] <- tryCatch({  
      estimate_parameters(y_t, v_t)  
    }, error = function(e) {  
      stop(paste("Error in estimate_parameters at simulation", sim, ":", e$message))  
    })  
    
    # Step 5: Perform MART estimation
    estimation <- tryCatch({
      MART(y_t, NULL, 1, 1, 0.25)
    }, error = function(e) {
      message(paste("Error in MART at simulation", sim, ":", e$message))
      return(NULL)  # return NULL to indicate failure
    })
    
    # Skip to next sim if estimation failed
    if (is.null(estimation)) {
      next
    }
    
    # Step 6: Extract coefficients  
    vec <- tryCatch({  
      c(estimation$coef.c1, estimation$coef.c2, estimation$coef.nc1,  
        estimation$coef.nc2, estimation$coef.exo1, estimation$coef.exo2,  
        estimation$coef.int, estimation$scale, estimation$df)  
    }, error = function(e) {  
      stop(paste("Error extracting coefficients at simulation", sim, ":", e$message))  
    })  
    
    models[[sim]] <- as.numeric(vec)  
  }  
  
  colnames(estimates) <- c("Causal", "nonCausal")  
  cbind(estimates, models)  
}


#function for getting results from simulations
get_final_estimates <- function(simulation) {
  simulation <- as.data.frame(simulation)
  colnames(simulation) <- c("Causal", "nonCausal", "models")
  cols_to_use <- simulation[3]
  #print(cols_to_use)
  # Extract the list column (assuming the column is the only column in the data frame)
  cols_to_use_list <- cols_to_use[[1]]
  # Convert the list into a matrix (each list element becomes a row)
  sim_matrix <- do.call(rbind, cols_to_use_list)
  cols_to_use <- as.data.frame(sim_matrix)
  colnames(sim_matrix) <- c("lag 1 regime 1", "lag 1 regime 2", "lead 1 regime 1", "lead 1 regime 2", "exo1", "exo2", "int", "scale", "df")
  #print(sim_matrix)
  # Calculate col means
  means <- apply(sim_matrix, 2, mean)
  sds <- apply(sim_matrix, 2, sd)
  #print(sds)
  
  
  mean_estimates <- c(means[1] ,means[2], means[3], means[4])
  sd_estimates <- c(sds[1], sds[2], sds[3], sds[4])
  df_estimates <- means[9]
  scale_estimates <- means[8]
  
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
n_sim <- 10
sample_sizes <- c(300, 500, 800)
param_combinations <- list(
  c(0.9, 0.9, 0.3, 0.3),
  c(0.9, 0.1, 0.3, 0.7),
  c(0.7, 0.3, 0.1, 0.9)
)

#hier ff mooie loop maken met sample sizes and param combinations
simulation_estimates <- matrix(NA, nrow = length(sample_sizes) * length(param_combinations), ncol = 15)
colnames(simulation_estimates) <- c("Sample_size", "phi_1", "phi_2", "psi_1", "psi_2", "lag_est_1", "lag_sd_1", "lag_est_2", "lag_sd_2", "lead_est_1", "lead_sd_1", "lead_est_2", "lead_sd_2", "df", "scale")
i <- 1

for (T in sample_sizes) {
  for (params in param_combinations) {
    phi_1 <- params[1]
    phi_2 <- params[2]
    psi_1 <- params[3]
    psi_2 <- params[4]
    
    simulation <- monte_carlo_simulation(T, phi_1, phi_2, psi_1, psi_2, n_sim)
    estimates <- get_final_estimates(simulation)
    
    simulation_estimates[i,1] = T
    simulation_estimates[i,2] = phi_1
    simulation_estimates[i,3] = phi_2
    simulation_estimates[i,4] = psi_1
    simulation_estimates[i,5] = psi_2
    simulation_estimates[i,6] = as.numeric(estimates$mean[1]) # Lag estimate regime 1
    simulation_estimates[i,7] = as.numeric(estimates$sd[1])  # Lag standard deviation regime 1
    simulation_estimates[i,8] = as.numeric(estimates$mean[2]) # Lag estimate regime 2
    simulation_estimates[i,9] = as.numeric(estimates$sd[2])  # Lag standard deviation regime 2
    simulation_estimates[i,10] = as.numeric(estimates$mean[3])  # Lead estimate regime 1
    simulation_estimates[i,11] = as.numeric(estimates$sd[3])  # Lead standard deviation regime 1
    simulation_estimates[i,12] = as.numeric(estimates$mean[4])  # Lead estimate regime 2
    simulation_estimates[i,13] = as.numeric(estimates$sd[4])  # Lead standard deviation regime 2
    simulation_estimates[i,14] = as.numeric(estimates$df)  # df value
    simulation_estimates[i,15] = as.numeric(estimates$scale)  # scale value
    
    i <- i+1
  }
}
view(simulation_estimates)
