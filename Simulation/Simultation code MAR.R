# Code voor simulatie van model.
set.seed(321) # Set seed for reproducibility

# Load required library
library(MASS)
library(pbmcapply)
source("MARX_functions.R")
source("MART.R")

d <- 1

# STEP 1: series from causal AR(p) model is generated.
ar_process <- function(phi, T, c = 0.5) {
  epsilon_t <- rt(T + r, df = 3) * sigma  # t-distributed errors with df=3 en keer sigma
  v <- rep(0, T + r - 1)  
  for (t in (r + 1):T) {
    c <- as.numeric(c)
    v[t] <- phi * v[t - 1] + epsilon_t[t]
  }
  return(v) # [-(1:r)] # Remove the first r elements and return the rest
}

# STEP 2: series from non causal component is generated
nc_process <- function(psi, v, T, c = 0.5) {
  y <- rep(0, T + s - 1)  
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
monte_carlo_simulation <- function(T, phi, psi, num_simulations = 10000, c = 0.5) {  
  estimates <- matrix(NA, nrow = num_simulations, ncol = 2)  # Storage for estimates  
  models <- vector("list", num_simulations)  # Initialize models as a list  
  
  for (sim in 1:num_simulations) {  
    #message(paste("Simulation:", sim))  # Track simulation number  
    
    # Step 1: Generate AR process  
    v_t <- tryCatch({  
      ar_process(phi, T, c)  
    }, error = function(e) {  
      stop(paste("Error in ar_process at simulation", sim, ":", e$message))  
    })  
    
    # Step 2: Generate NC process  
    y_t <- tryCatch({  
      nc_process(psi, v_t, T, c)  
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
      marx.t(y_t,NULL,1,1) 
    }, error = function(e) {  
      stop(paste("Error in MART at simulation", sim, ":", e$message))  
    })  
    
    # Step 6: Extract coefficients  
    vec <- tryCatch({  
      c(estimation$coef.c, estimation$coef.nc, estimation$coef.exo,
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
  # Extract the list column (assuming the column is the only column in the data frame)
  cols_to_use_list <- cols_to_use[[1]]
  # Convert the list into a matrix (each list element becomes a row)
  sim_matrix <- do.call(rbind, cols_to_use_list)
  cols_to_use <- as.data.frame(sim_matrix)
  colnames(sim_matrix) <- c("lag", "lead", "exo", "int", "scale", "df")

  # Calculate col means
  means <- apply(sim_matrix, 2, mean)
  sds <- apply(sim_matrix, 2, sd)
  
  mean_estimates <- c(means[1] ,means[2])
  sd_estimates <- c(sds[1], sds[2])
  df_estimates <- means[6]
  scale_estimates <- means[5]
  
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
n_sim <- 10000
sample_sizes <- c(300, 500, 800)
param_combinations <- list(
  c(0.9, 0.9),
  c(0.9, 0.1),
  c(0.1, 0.9)
)
# >>> CHANGED >>>
# Create parameter grid
param_grid <- expand.grid(
  T = sample_sizes,
  phi_psi = param_combinations,
  stringsAsFactors = FALSE
)

# Run the simulation loop with progress bar
results_list <- pbmclapply(1:nrow(param_grid), mc.cores = parallel::detectCores() - 1, FUN = function(i) {
  T <- param_grid$T[i]
  phi <- param_grid$phi_psi[[i]][1]
  psi <- param_grid$phi_psi[[i]][2]
  
  simulation <- monte_carlo_simulation(T, phi, psi, n_sim)
  estimates <- get_final_estimates(simulation)
  
  result_row <- numeric(9)
  result_row[1] <- T
  result_row[2] <- phi
  result_row[3] <- psi
  result_row[4] <- as.numeric(estimates$mean[1])
  result_row[5] <- as.numeric(estimates$sd[1])
  result_row[6] <- as.numeric(estimates$mean[2])
  result_row[7] <- as.numeric(estimates$sd[2])
  result_row[8] <- as.numeric(estimates$df)
  result_row[9] <- as.numeric(estimates$scale)
  
  return(result_row)
})
# <<< CHANGED <<<

# >>> CHANGED >>>
# Bind result matrix
simulation_estimates <- do.call(rbind, results_list)
colnames(simulation_estimates) <- c("Sample_size", "phi", "psi", "lag_est", "lag_sd", "lead_est", "lead_sd", "df", "scale")
# <<< CHANGED <<<

print(dim(simulation_estimates))