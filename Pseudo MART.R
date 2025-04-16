# MART pseudo function

# Load required libraries
library(MASS)          # For statistical distributions and matrix functions
source("MARX_functions.R")  # Custom functions for MAR model estimation
source("MART.R")            # MART model training and forecasting routines (includes information criteria calculations)
library(forecast)      # For ARIMA modeling and forecast tools
library(pbmcapply)     # For parallel processing with progress bar

#' @title The model selection for pseudo-ARX function
#' @description This function allows you to calculate AIC, BIC, HQ for pseudo-ARX models.
#' @param y Data vector of time series observations.
#' @param x Matrix of data (every column represents one time series). Specify NULL or "not" if not wanted.
#' @param p_max Maximum number of autoregressive terms to be included.
#' @keywords selection
#' @keywords pseudo-causal
#' @return \item{bic}{Vector containing values BIC for p=0 up to p_max.}
#' @return \item{aic}{Vector containing values AIC for p=0 up to p_max.}
#' @return \item{hq}{vector containing values HQ for p=0 up to p_max.}
#' @author Sean Telg
#' @export
#' @examples
#' data <- sim.marx(c('t',1,1), c('t',1,1),100,0.5,0.4,0.3)
#' selection.lag(data$y,data$x,8)

selection.lag_t <- function(y,x,p_max,c,d=3){
  if (is.null(x)){
    x <- "not"
  }
  
  bic_results <- bic(y,x,p_max,c,d=3)
  aic_results <- aic(y,x,p_max,c,d=3)
  hq_results <- hq(y,x,p_max,c,d=3)
  
  bic_vec <- bic_results[[2]]
  colnames(bic_vec) <- paste('p =', 0:p_max)
  aic_vec <- aic_results[[2]]
  colnames(aic_vec) <- paste('p =', 0:p_max)
  hq_vec <- hq_results[[2]]
  colnames(hq_vec) <- paste('p =', 0:p_max)
  
  cat('Order Selection Criteria Pseudo Causal Model:', "\n")
  cat(' ', "\n")
  cat('BAYESIAN INFORMATION CRITERION', "\n")
  print(bic_vec)
  cat(' ', "\n")
  cat('Minimum value attained at p = ')
  cat(which.min(bic_vec) -1)
  cat(' ',  "\n")
  cat(' ',  "\n")
  cat('AKAIKE INFORMATION CRITERION', "\n")
  print(aic_vec)
  cat(' ', "\n")
  cat('Minimum value attained at p = ')
  cat(which.min(aic_vec) - 1)
  cat(' ',  "\n")
  cat(' ',  "\n")
  cat('HANNAN-QUINN INFORMATION CRITERION', "\n")
  print(hq_vec)
  cat(' ', "\n")
  cat('Minimum value attained at p = ')
  cat(which.min(hq_vec) - 1)
  cat(' ', "\n")
  cat(' ', "\n")
  
  return(list(bic = bic_vec,aic=aic_vec,hq=hq_vec))
}

#' @title The Bayesian/Schwarz information criterion (BIC) function
#' @description This function allows you to calculate the Bayesian/Schwarz information criteria (BIC) for ARX models.
#' @param y Data vector of time series observations.
#' @param x Matrix of data (every column represents one time series). Specify NULL or "not" if not wanted.
#' @param p_max Maximum number of autoregressive terms to be included.
#' @keywords selection
#' @return \item{p}{Lag order chosen by BIC.}
#' @return \item{values}{Vector containing values BIc for p = 0 up to p_max.}
#' @author Sean Telg
#' @export
#' @examples
#' data <- sim.marx(c('t',1,1), c('t',1,1),100,0.5,0.4,0.3)
#' bic(data$y, data$x,8)

bic <- function(y,x,p_max,c,d=1){
  
  if (is.null(x)){
    x <- "not"
  }
  
  y <- fBasics::vec(y)
  y <- y - mean(y)
  # Demean c om de threshold consistent te houden
  c <- c - mean(y)
  n <- length(y) - max(p_max,d)
  
  if (length(x) > 1){
    numcol <- NCOL(x)
  }
  else{
    numcol = 0
  }
  
  crit <- matrix(data=NA, nrow=(p_max+1), ncol=1)
  
  for (p in 0:p_max){
    
    arx.ls_T_results <- arx.ls_T(y,x,p,c,d)
    n <- length(arx.ls_T_results[[5]])
    Cov <- arx.ls_T_results[[6]]
    crit[(p+1)] <- -2*Cov/n + ((log(n))/n)*(2*p+1+2*numcol)
  }
  
  p_bic <- which.min(crit) - 1
  
  crit <- t(crit)
  colnames(crit) <- paste('p =', 0:p_max)
  
  return(list(p = p_bic, values= crit))
}


#' @title The Akaike information criterion (AIC) function
#' @description This function allows you to calculate the Akaike information criteria (AIC) for ARX models.
#' @param y Data vector of time series observations.
#' @param x Matrix of data (every column represents one time series). Specify NULL or "not" if not wanted.
#' @param p_max Maximum number of autoregressive terms to be included.
#' @keywords selection
#' @return \item{p}{Lag order chosen by AIC.}
#' @return \item{values}{Vector containing values AIC for p = 0 up to p_max.}
#' @author Sean Telg
#' @export
#' @examples
#' data <- sim.marx(c('t',1,1), c('t',1,1),100,0.5,0.4,0.3)
#' aic(data$y, data$x,8)

aic <- function(y,x,p_max,c,d=3){
  
  if (is.null(x)){
    x <- "not"
  }
  
  y <- fBasics::vec(y)
  y <- y - mean(y)
  # Demean c om de threshold consistent te houden
  c <- c - mean(y)
  n <- length(y) - max(p_max,d)
  
  if (length(x) > 1){
    numcol <- NCOL(x)
  }
  else{
    numcol = 0
  }
  
  crit <- matrix(data=NA, nrow=(p_max+1), ncol=1)
  
  for (p in 0:p_max){
    
    arx.ls_T_results <- arx.ls_T(y,x,p,c,d)
    n <- length(arx.ls_T_results[[5]])
    Cov <- arx.ls_T_results[[6]]
    crit[(p+1)] <- -2*Cov/n + (2/n)*(2*p+1+2*numcol)
  }
  
  p_aic <- which.min(crit) - 1
  
  crit <- t(crit)
  colnames(crit) <- paste('p =', 0:p_max)
  
  return(list(p = p_aic,values = crit))
  
}

#' @title The Hannan-Quinn (HQ) information criterion function
#' @description This function allows you to calculate the Hannan-Quinn (HQ) information criteria for ARX models.
#' @param y       Data vector of time series observations.
#' @param x       Matrix of data (every column represents one time series). Specify NULL or "not" if not wanted.
#' @param p_max   Maximum number of autoregressive terms to be included.
#' @keywords selection
#' @return \item{p}{Lag order chosen by HQ.}
#' @return \item{values}{Vector containing values HQ for p = 0 up to p_max.}
#' @author Sean Telg
#' @export
#' @examples
#' data <- sim.marx(c('t',1,1), c('t',1,1),100,0.5,0.4,0.3)
#' hq(data$y, data$x,8)

hq <- function(y,x,p_max,c,d=3){
  
  if (is.null(x)){
    x <- "not"
  }
  
  y <- fBasics::vec(y)
  y <- y - mean(y)
  # Demean c om de threshold consistent te houden
  c <- c - mean(y)
  n <- length(y) - max(p_max,d)
  
  if (length(x) > 1){
    numcol <- NCOL(x)
  }
  else{
    numcol = 0
  }
  
  crit <- matrix(data=NA, nrow=(p_max+1), ncol=1)
  
  for (p in 0:p_max){
    
    arx.ls_T_results <- arx.ls_T(y,x,p,c,d)
    n <- length(arx.ls_T_results[[5]])
    Cov <- arx.ls_T_results[[6]]
    crit[(p+1)] <- -2*Cov/n + ((2*log(log(n)))/n)*(2*p+1+2*numcol)
  }
  
  p_hq <- which.min(crit) - 1
  
  crit <- t(crit)
  colnames(crit) <- paste('p =', 0:p_max)
  
  return(list(p = p_hq, values = crit))
  
}

#' #' @title The lag-lead model selection for MARX function
#' #' @description This function allows you to determine the MARX model (for p = r + s) that maximizes the t-log-likelihood.
#' #' @param y Data vector of time series observations.
#' #' @param x Matrix of data (every column represents one time series). Specify NULL or "not" if not wanted.
#' #' @param p_pseudo Number of autoregressive terms to be included in the pseudo-causal model.
#' #' @keywords selection
#' #' @keywords causal-noncausal
#' #' @return \item{p.C}{The number of lags selected.}
#' #' @return \item{p.NC}{The number of leads selected.}
#' #' @return \item{loglikelihood}{The value of the loglikelihood for all models with p = r + s.}
#' #' @author Sean Telg
#' #' @export
#' #' @examples
#' #' data <- sim.marx(c('t',3,1), c('t',3,1),100,0.5,0.4,0.3)
#' #' selection.lag.lead(data$y,data$x,2)
#' 
selection.lag.lead_T <- function(y, x, p_pseudo, c, d = 1) {
  y <- as.numeric(y)
  # Check if x is NULL and set it to 'not' if true
  if (is.null(x)) {
    x <- "not"
  }

  P_C <- as.numeric(seq(length = (p_pseudo + 1), from = 0, by = 1))
  P_C <- as.numeric(fBasics::vec(P_C))
  print(P_C)
  P_NC <- as.numeric(rev(P_C))
  P_NC <- as.numeric(fBasics::vec(P_NC))
  print(P_NC)

  
  loglik <- c()

  for (i in 1:(p_pseudo + 1)) {

    # Using tryCatch to handle potential errors during the MART call and the subsequent operations
    tryCatch({
      MART_results <- MART(y, x, P_C[i], P_NC[i], c, d)
      sig <- as.numeric(MART_results[[8]])
      df  <- as.numeric(MART_results[[9]])
      E   <- MART_results[[10]]
      n <- length(E)

      # Check if the components are numeric
      if (!is.numeric(E)) stop("E is not numeric.")
      if (!is.numeric(df)) stop("df is not numeric.")
      if (!is.numeric(sig)) stop("sig is not numeric.")

      loglik[i] <- (n*lgamma((df+1)/2) - n*log(sqrt(df*pi*sig^2)) - n*lgamma(df/2) - ((df+1)/2)*log(1+(E/sig)^2/df) %*% matlab::ones(n,1))


    }, error = function(e) {
      # Catch errors during this iteration and provide informative feedback
      cat("Error in iteration", i, "\n")
      cat("Error message: ", e$message, "\n")
      cat("Skipping to next iteration...\n")
      loglik[i] <- NA  # Assign NA to loglik in case of error
    })
  }

  # After the loop, check for a valid max log-likelihood index
  if (all(is.na(loglik))) {
    stop("All log-likelihoods failed. Cannot proceed.")
  }

  maxloglik <- which.max(loglik)
  print(maxloglik)

  if (!is.numeric(maxloglik)) stop("maxloglik is not numeric.")

  P <- cbind(P_C, P_NC)
  print(P[maxloglik, ])
  P <- fBasics::vec(P[maxloglik, ])

  if (!is.numeric(P)) stop("P is not numeric.")

  p_C  <- P[1]
  p_NC <- P[2]

  if (!is.numeric(p_C)) stop("p_C is not numeric.")
  if (!is.numeric(p_NC)) stop("p_NC is not numeric.")
  if (!is.numeric(rev(loglik))) stop("rev(loglik) is not numeric.")

  return(list(p.C = p_C, p.NC = p_NC, loglikelihood = rev(loglik)))
}


selection.lag_t(inflation_df_monthly$inflationNonSA,inflation_df_monthly$UNRATE,12,median(inflation_df_monthly$inflationNonSA),d=3)

p_pseudo <- readline(prompt = "Choose lag order for pseudo causal model: ")
p_pseudo <- as.numeric(p_pseudo)

pseudo <- arx.ls_T(inflation_df_monthly$inflationNonSA,GS1,p_pseudo,median(inflation_df_monthly$inflationNonSA),d=3)
Cov_pseudo <- pseudo[[4]]
U_pseudo <- pseudo[[5]]
test_cdf_pseudo <- cbind(U_pseudo, stats::pnorm(U_pseudo,0,Cov_pseudo))

kstest_results <- stats::ks.test(test_cdf_pseudo[,1],"pnorm",0,Cov_pseudo)
jarquebera     <- tseries::jarque.bera.test(U_pseudo)

if (kstest_results$p.value < 0.05){
  hh_pseudo = 1 			## reject
}else{
  hh_pseudo = 0 			## not reject
}


if (jarquebera$p.value < 0.05){
  jarque_check = 1
}else{
  jarque_check = 0
}


if (hh_pseudo == 0){
  
  cat(' ', "\n")
  cat(' ', "\n")
  cat('THE KS-TEST FAILS TO REJECT THE NULL OF NORMALLY DISTRIBUTED RESIDUALS OF THE PURELY CAUSAL ARX MODEL', "\n")
  cat('p-value:')
  cat(kstest_results$p.value, "\n")
  cat('WARNING: MIxED ARX MODEL MIGHT NOT BE IDENTIFIABLE!', "\n")
  cat(' ', "\n")
  cat(' ', "\n")
}else{
  
  cat(' ', "\n")
  cat(' ', "\n")
  cat('THE KS-TEST REJECTS THE NULL OF NORMALLY DISTRIBUTED RESIDUALS OF THE PURELY CAUSAL ARX MODEL', "\n")
  cat('p-value:')
  cat(kstest_results$p.value, "\n")
  cat(' ', "\n")
  cat(' ', "\n")
}


if (jarque_check == 0){
  
  cat(' ', "\n")
  cat(' ', "\n")
  cat('THE JB-TEST FAILS TO REJECT THE NULL OF NORMALLY DISTRIBUTED RESIDUALS OF THE PURELY CAUSAL ARX MODEL', "\n")
  cat('p-value:')
  cat(jarquebera$p.value, "\n")
  cat('WARNING: MIxED ARX MODEL MIGHT NOT BE IDENTIFIABLE!', "\n")
  cat(' ', "\n")
  cat(' ', "\n")
}else{
  
  cat(' ', "\n")
  cat(' ', "\n")
  cat('THE JB-TEST REJECTS THE NULL OF NORMALLY DISTRIBUTED RESIDUALS OF THE PURELY CAUSAL ARX MODEL', "\n")
  cat('p-value:')
  cat(jarquebera$p.value, "\n")
  cat(' ', "\n")
  cat(' ', "\n")
}

stats::qqnorm(U_pseudo, main="Normal Probability Plot of Residuals")
stats::qqline(U_pseudo)

selection.lag.lead_results <- selection.lag.lead_T(inflation_df_monthly$inflationNonSA,GS1,p_pseudo,median(inflation_df_monthly$inflationNonSA),d=3)
p_C <- selection.lag.lead_results[[1]]
p_NC <- selection.lag.lead_results[[2]]


## -----------------------------------------------------------------------------
# Residual diagnostics: test for independence of AR(p) residuals (Hecq et al. 2016) and test for no serial correlation (MARX package paper of HEcq et al.)
# -----------------------------------------------------------------------------

# Fit a 12-lag AR model to the inflation series
model_ar2 <- Arima(inflation_df_monthly$inflationNonSA, order = c(2, 0, 0))
resids_ar2 <- model_ar2$residuals  # Extract residuals

# Step 2: Square the residuals for use as regressors
resids_sq <- resids_ar2^2

# Step 3: Create lag matrix manually
m <- 2
n <- length(resids_ar2)

# Create the response variable y (residuals from t = m+1 to n)
y <- resids_ar2[(m + 1):n]

# Create lagged squared residuals matrix
X_lags <- matrix(NA, nrow = n - m, ncol = m)
for (i in 1:m) {
  X_lags[, i] <- resids_sq[(m + 1 - i):(n - i)]
}

# Step 4: Regress current residual on lagged squared residuals
model_test <- lm(y ~ X_lags)

# Step 5: Test for joint significance of lag coefficients (H0: residuals are i.i.d.)
test_statistic <- summary(model_test)$r.squared * length(y)
p_value <- pchisq(test_statistic, df = m, lower.tail = FALSE)

# Step 6: Output results of the chi-squared test
cat("Chi-squared test statistic:", test_statistic, "\n")
cat("p-value:", p_value, "\n")

# Step 1: Perfome Ljung-Box test
Box.test(resids_ar2, lag = 6, type = "Ljung-Box")

