
## -----------------------------------------------------------------------------
# Testing if the coefficients differ significantly in the regimes 
# -----------------------------------------------------------------------------

# MART 
# Final model is MART(1,1) with median as threshold and d = 3

model.MART <- MART(inflation_df_monthly$inflationNonSA,NULL,1,1,median(inflation_df_monthly$inflationNonSA),3)
model.ART <- MART(inflation_df_monthly$inflationNonSA,NULL,2,0,median(inflation_df_monthly$inflationNonSA),3)
# Define the function for the MART(1,1)
calculate_z_p_values <- function(model.MART) {
  # Calculate z-values
  z_c1 <- model.MART$coef.c1 / model.MART$se.dist[1]
  z_c2 <- model.MART$coef.c2 / model.MART$se.dist[2]
  z_nc1 <- model.MART$coef.nc1 / model.MART$se.dist[3]
  z_nc2 <- model.MART$coef.nc2 / model.MART$se.dist[4]
  
  # Calculate p-values
  p_c1 <- 2 * (1 - pnorm(abs(z_c1)))
  p_c2 <- 2 * (1 - pnorm(abs(z_c2)))
  p_nc1 <- 2 * (1 - pnorm(abs(z_nc1)))
  p_nc2 <- 2 * (1 - pnorm(abs(z_nc2)))
  
  # Return z-values and p-values as a list
  return(list(
    z_c1 = z_c1,
    z_c2 = z_c2,
    z_nc1 = z_nc1,
    z_nc2 = z_nc2,
    p_c1 = p_c1,
    p_c2 = p_c2,
    p_nc1 = p_nc1,
    p_nc2 = p_nc2
  ))
}

calculate_z_p_values_ar <- function(model.ART) {
  # Calculate z-values
  z_c1_1 <- model.ART$coef.c1[1] / model.ART$se.dist[1]
  z_c2_1 <- model.ART$coef.c2[1] / model.ART$se.dist[3]
  z_c1_2 <- model.ART$coef.c1[2] / model.ART$se.dist[2]
  z_c2_2 <- model.ART$coef.c2[2] / model.ART$se.dist[4]

  
  # Calculate p-values
  p_c1_1 <- 2 * (1 - pnorm(abs(z_c1_1)))
  p_c2_1 <- 2 * (1 - pnorm(abs(z_c2_1))) 
  p_c1_2 <- 2 * (1 - pnorm(abs(z_c1_2)))
  p_c2_2 <- 2 * (1 - pnorm(abs(z_c2_2)))

  
  # Return z-values and p-values as a list
  return(list(
    z_c1_1 = z_c1_1,
    z_c2_1 = z_c2_1,
    z_c1_2 = z_c1_2,
    z_c2_2 = z_c2_2,
    p_c1_1 = p_c1_1,
    p_c2_1 = p_c2_1,
    p_c1_2 = p_c1_2,
    p_c2_2 = p_c2_2
  ))
}
# Call function for AR(2)
result_ar <- calculate_z_p_values_ar(model.ART)
result_ar

# Call the function with the model object
result <- calculate_z_p_values(model.MART)

# View the result
result

# Testing for causual part
diff_c <- model.MART$coef.c1 - (model.MART$coef.c2)
se_c <- sqrt(model.MART$se.dist[1]^2 + model.MART$se.dist[2]^2) #- 2*model.MART$var.cov[1,2]
z_c <- diff_c / se_c
p_value_c <- 2 * (1 - pnorm(abs(z_c)))
cat("Causal z:", z_c, "\nCausal p-value:", p_value_c, "\n")

# Testing for non causual part
diff_nc <- model.MART$coef.nc1 - (model.MART$coef.nc2)
se_nc <- sqrt(model.MART$se.dist[3]^2 + model.MART$se.dist[4]^2) #- 2*model.MART$var.cov[3,4]
z_nc <- diff_nc / se_nc
p_value_nc <- 2 * (1 - pnorm(abs(z_nc)))
cat("Non Causal z:", z_nc, "\n Non Causal p-value:", p_value_nc, "\n")

wald_test <- function(beta1, beta2, varcov_matrix, idx1, idx2) {
  # Compute difference vector
  d <- beta1 - beta2
  
  # Construct variance of difference using covariances
  var_d <- varcov_matrix[idx1, idx1] +
    varcov_matrix[idx2, idx2] -
    varcov_matrix[idx1, idx2] -
    varcov_matrix[idx2, idx1]
  
  # If d is vector, result is quadratic form (Wald stat)
  W <- t(d) %*% solve(var_d) %*% d
  
  # Degrees of freedom = length of vector
  df <- length(d)
  
  # p-value from chi-squared distribution
  p_value <- 1 - pchisq(W, df)
  
  return(list(statistic = as.numeric(W), df = df, p_value = p_value))
}

# Testing for causal part in AR(2)
b_1 <- c(model.ART$coef.c1[1],model.ART$coef.c1[2])
b_2 <- c(model.ART$coef.c2[1],model.ART$coef.c2[2])
cov <- model.ART$var.cov

ar_wald <- wald_test(b_1, b_2, cov, idx1 = c(1, 2), idx2 = c(3, 4))
ar_wald

