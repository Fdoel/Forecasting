#' Mixed causal-noncausal autoregressions with exogenous regressors.
#'
#' @title The MARX function
#' @description This interface-based function allows you to perform model selection for MARX models based on information criteria.
#' @param y Data vector of time series observations.
#' @param x Matrix of data (every column represents one time series). Specify NULL or "not" if not wanted.
#' @param p_max Maximum number of autoregressive parameters (leads + lags) to be included.
#' @param sig_level Significance level for the construction of inference.
#' @param p_C Number of lags (if not specified by the user a model selection procedure is used to determine the number of lags).
#' @param p_NC Number of leads (if not specified by the user a model selection procedure is used to determine the number of leads).
#' @keywords estimation
#' @keywords selection
#' @return The function returns the values of the information criteria for the pseudo-causal models. The user is asked to choose a value for "p". Extensive output for the MARX(r,s,q) model (with p = r + s) which maximizes the log-likelihood is reported.
#' @author Sean Telg
#' @export
#' @examples
#' data <- sim.marx(c('t',1,1), c('t',1,1),100,0.5,0.4,0.3)
#' p_max <- 8
#' sig_level <- 0.05
#' marx(data$y, data$x, p_max, sig_level,1,1) ## p_C and p_NC chosen to be 1: MARX(1,1,1) output.
#' marx(data$y, NULL, p_max,sig_level,1,1) ## MAR(1,1), no exogenous variable specified.

marx <- function(y,x,p_max,sig_level,p_C,p_NC){

  if (is.null(x)){
    x <- "not"
  }

  print(match.call())
  nargin <- length(as.list(match.call())) -1

  if (nargin == 4){
    p = 1
  }
  else{
    p = 0
  }

  if (p == 1){

    selection.lag(y,x,p_max)
    p_pseudo <- readline(prompt = "Choose lag order for pseudo causal model: ")
    p_pseudo <- as.numeric(p_pseudo)

    pseudo <- arx.ls(y,x,p_pseudo)
    Cov_pseudo <- pseudo[[4]]
    U_pseudo <- pseudo[[5]]
    test_cdf_pseudo <- cbind(U_pseudo, stats::pnorm(U_pseudo,0,Cov_pseudo))

    kstest_results <- stats::ks.test(test_cdf_pseudo[,1],"pnorm",0,Cov_pseudo)
    jarquebera     <- tseries::jarque.bera.test(U_pseudo)

    if (kstest_results$p.value < 0.05){
      hh_pseudo = 1 			## reject
    }
    else{
      hh_pseudo = 0 			## not reject
    }


    if (jarquebera$p.value < 0.05){
      jarque_check = 1
    }
    else{
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
    }
    else{

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
    }
    else{

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

    selection.lag.lead_results <- selection.lag.lead(y,x,p_pseudo)
    p_C <- selection.lag.lead_results[[1]]
    p_NC <- selection.lag.lead_results[[2]]

  }

  marx.t_results <- marx.t(y,x,p_C,p_NC);
  B_C  <- marx.t_results[[1]]
  B_NC <- marx.t_results[[2]]
  B_x  <- marx.t_results[[3]]
  IC   <- marx.t_results[[4]]
  sig  <- marx.t_results[[5]]
  df   <- marx.t_results[[6]]

  inference_results <- inference(y,x,B_C,B_NC,B_x,IC,sig,df,sig_level)
  BC     <- inference_results[[1]]
  BNC    <- inference_results[[2]]
  Bx     <- inference_results[[3]]
  ICC    <- inference_results[[4]]
  std_c  <- inference_results[[5]]
  std_nc <- inference_results[[6]]
  std_x  <- inference_results[[7]]
  std_ic <- inference_results[[8]]

  cat(' ', "\n")
  cat(' ', "\n")
  cat('Estimated Causal Parameters: ', "\n")
  print(B_C)
  cat(' ', "\n")
  cat('Corresponding Standard Errors: ', "\n")
  print(std_c)
  cat(' ', "\n")
  cat('Corresponding Confidence Intervals: ', "\n")
  print(BC)
  cat(' ', "\n")
  cat(' ', "\n")

  cat('Estimated Noncausal Parameters: ', "\n")
  print(B_NC)
  cat(' ', "\n")
  cat('Corresponding Standard Errors: ', "\n")
  print(std_nc)
  cat(' ', "\n")
  cat('Corresponding Confidence Intervals: ', "\n")
  print(BNC)
  cat(' ', "\n")
  cat(' ', "\n")

  cat('Estimated Parameters Exogenous Variables: ', "\n")
  print(B_x)
  cat(' ', "\n")
  cat('Corresponding Standard Errors: ', "\n")
  print(std_x)
  cat(' ', "\n")
  cat('Corresponding Confidence Intervals: ', "\n")
  print(Bx)
  cat(' ', "\n")
  cat(' ', "\n")


  cat('Estimated Intercept: ', "\n")
  print(IC)
  cat(' ', "\n")
  cat('Corresponding Standard Errors: ', "\n")
  print(std_ic)
  cat(' ', "\n")
  cat('Corresponding Confidence Intervals: ', "\n")
  print(ICC)
  cat(' ', "\n")
  cat(' ', "\n")

  cat('Estimated Distributional Parameters (df,sig): ', "\n")
  print(c(df,sig))
  cat(' ', "\n")
  cat(' ', "\n")
}

#' @title The regressor matrix function
#' @description This function allows you to create a regressor matrix.
#' @param y   Data vector of time series observations.
#' @param x   Matrix of data (every column represents one time series). Specify NULL or "not" if not wanted.
#' @param p   Number of autoregressive terms to be included.
#' @keywords estimation
#' @return    \item{Z}{Regressor matrix}
#' @author Sean Telg
#' @export
#' @examples
#' data <- sim.marx(c('t',3,1),c('t',1,1),100,0.5,0.4,0.3)
#' regressor.matrix(data$y, data$x, 2)

regressor.matrix <- function(y,x,p){

  if (is.null(x)){
    x <- "not"
  }

  y <- fBasics::vec(y)

  n <- length(y)

  if (p==1){
    k <-1
  }
  else{
    k <- NCOL(y)
  }


  if (p > 0){
    Z <- matlab::zeros(n,k*p)

    for (i in 1:p){
      Z[(1+i):n,((i-1)*k+1):(i*k)] <- y[1:(n-i)]
    }

    Z <- Z[(1+p):n,]

  }
  else{
    Z <- matrix(,nrow=n,ncol=0)
  }

  if (x == "not" && length(x) == 1){
    Z <- Z
  }


  if (NCOL(x) == 1 && x != "not"){
    Z <- cbind(Z,x[(1+p):n])
  }
  else if (NCOL(x) > 1 && x != "not"){
    Z <- cbind(Z,x[(1+p):n,])
  }


  return(matrix = Z)
}

#' @title The ARX estimation by OLS function
#' @description This function allows you to estimate ARX models by ordinary least squares (OLS).
#' @param y Data vector of time series observations.
#' @param x Matrix of data (every column represents one time series). Specify NULL or "not" if not wanted.
#' @param p Number of autoregressive terms to be included.
#' @keywords estimation
#' @keywords pseudo-causal
#' @return \item{coefficients}{Vector of estimated coefficients.}
#' @return \item{coef.auto}{Vector of estimated autoregressive parameters.}
#' @return \item{coef.exo}{Vector of estimated exogenous parameters.}
#' @return \item{mse}{Mean squared error.}
#' @return \item{residuals}{Residuals.}
#' @return \item{loglikelihood}{Value of the loglikelihood.}
#' @return \item{fitted.values}{Fitted values.}
#' @return \item{df}{Degrees of freedom.}
#' @return \item{vcov}{Variance-covariance matrix of residuals.}
#' @author Sean Telg
#' @export
#' @examples
#' data <- sim.marx(c('t',3,1),c('t',1,1),100,0.5,0.4,0.3)
#' arx.ls(data$y,data$x,2)

arx.ls <- function(y,x,p){

  if (is.null(x)){
    x <- "not"
  }

  n <- length(y) - p

  Y <- y[(p+1):length(y)]
  int <- rep(1,(length(y)-p))
  Z <- regressor.matrix(y,x,p)
  Z <- cbind(int,Z)

  df <- nrow(Z) - NCOL(Z)

  B <- solve(t(Z) %*% Z) %*% (t(Z) %*% Y)

  if (p > 0){
    if (length(x) > 1){
      rownames(B) <- c('int', paste('lag', 1:p), paste('exo', 1:NCOL(x)))
    }
    else{
      rownames(B) <- c('int', paste('lag', 1:p))
    }
  }
  else{
    if (length(x) > 1){
      rownames(B) <- c('int', paste('exo', 1:NCOL(x)))
    }
    else{
      rownames(B) <- 'int'
    }
  }

  FV <- Z %*% B
  U <- Y - FV

  sig <- (t(U) %*% U)
  sig <- as.numeric(sig)

  Cov <- (1/n)*sig
  Cov <- as.numeric(Cov)

  sigma2 <- sum((Y - Z %*% B)^2)/df
  qz <- qr(Z)
  vcov <- sigma2*chol2inv(qz$qr)
  colnames(vcov) <- rownames(vcov) <- colnames(Z)

  Loglik <- -(n/2)*(1 + log(2*pi)+log(Cov))

  if (p == 0){
    B_auto <- 0
  }
  else{
    B_auto <- B[2:(p+1)]
  }

  if (length(x) > 1){
    B_x <- B[(p+2):length(B)]
  }
  else{
    B_x <- 0
  }

  return(list(coefficients = B, coef.auto = B_auto, coef.exo = B_x, mse = Cov, residuals = U, loglikelihood = Loglik, fitted.values = FV, df = df,vcov=vcov))
}

#' @title The estimation of the MARX model by t-MLE function
#' @description This function allows you to estimate the MARX model by t-MLE.
#' @param y Data vector of time series observations.
#' @param x Matrix of data (every column represents one time series). Specify NULL or "not" if not wanted.
#' @param p_C Number of lags.
#' @param p_NC Number of leads.
#' @param params0 Starting values for the parameters to be estimated (both model and distributional parameters).
#' @keywords estimation
#' @keywords causal-noncausal
#' @return \item{coef.c}{Estimated causal coefficients.}
#' @return \item{coef.nc}{Estimated noncausal coefficients.}
#' @return \item{coef.exo}{Estimated exogenous coefficients.}
#' @return \item{coef.int}{Estimated intercept.}
#' @return \item{scale}{Estimated scale parameter.}
#' @return \item{df}{Estimated degrees of freedom.}
#' @return \item{residuals}{Residuals.}
#' @return \item{se.dist}{Standard errors of the distributional parameters.}
#' @author Sean Telg
#' @export
#' @examples
#' data <- sim.marx(c('t',3,1),c('t',3,1),100,0.5,0.4,0.3)
#' marx.t(data$y,data$x,1,1)

marx.t <- function(y,x,p_C,p_NC,params0){

  #print(match.call())
  nargin <- length(as.list(match.call())) - 1

  if (is.null(x)){
    x <- "not"
  }

  if (length(x) == 1){
    numcol <- 0
  }
  else{
    numcol <- NCOL(x)
  }

  if(numcol > 1){
    x.rev <- matrix(data=NA,nrow=length(x[,1]),ncol=numcol)
    for (i in 1:numcol){
      x.rev[,i] <- rev(x[,i])
    }
  }
  else{
    x.rev <- matrix(data=NA,nrow=length(x),ncol=numcol)
    x.rev <- rev(x)
  }


  if (nargin < 5){
    y    <- fBasics::vec(y)
    z    <- rev(y)
    z    <- fBasics::vec(z)
    BC0  <- arx.ls(y,x,p_C)[[2]]
    Bx0  <- arx.ls(y,x,p_C)[[3]]
    BNC0 <- arx.ls(z,x.rev,p_NC)[[2]]
    IC0  <- 0
    df0  <- 20
    sig0 <- 2

    BC0 <- fBasics::vec(BC0)
    BNC0 <- fBasics::vec(BNC0)
    Bx0 <- fBasics::vec(Bx0)

  if (length(x) > 1){
      if (p_C > 0 && p_NC > 0){
        params0 <- rbind(BC0,BNC0,Bx0,IC0,sig0,df0)
      }
      else if (p_NC > 0 && p_C == 0){
        params0 <- rbind(BNC0,Bx0,IC0,sig0,df0)
      }
      else if (p_C > 0 && p_NC == 0){
        params0 <- rbind(BC0,Bx0,IC0,sig0,df0)
      }
      else if (p_C == 0 && p_NC == 0){
        params0 <- rbind(Bx0,IC0,sig0,df0)
      }
  }
  else{
    if (p_C > 0 && p_NC > 0){
      params0 <- rbind(BC0,BNC0,IC0,sig0,df0)
    }
    else if (p_NC > 0 && p_C == 0){
      params0 <- rbind(BNC0,IC0,sig0,df0)
    }
    else if (p_C > 0 && p_NC == 0){
      params0 <- rbind(BC0,IC0,sig0,df0)
    }
    else if (p_C == 0 && p_NC == 0){
      params0 <- rbind(IC0,sig0,df0)
    }
  }
}

  optimization_results <- stats::optim(params0,ll.max,gr=NULL,y=fBasics::vec(y),p_C=p_C,p_NC=p_NC,x=x,method="BFGS",hessian=TRUE)
  PARAMS <- optimization_results$par

  if (length(x) > 1){
    numcol <- NCOL(x)

    if (p_C > 0 && p_NC > 0){
        B_C  <- PARAMS[1:p_C]
        B_NC <- PARAMS[(p_C+1):(p_C + p_NC)]
        B_x  <- PARAMS[(p_C + p_NC + 1):(p_C + p_NC + numcol)]
        IC   <- PARAMS[(p_C + p_NC + numcol + 1)]
        sig  <- PARAMS[(p_C + p_NC + numcol + 2)]
        df   <- PARAMS[(p_C + p_NC + numcol + 3)]
      }
      else if (p_NC > 0 && p_C == 0){
        B_C  <- 0
        B_NC <- PARAMS[1:p_NC]
        B_x  <- PARAMS[(p_NC + 1):(p_NC + numcol)]
        IC   <- PARAMS[(p_NC + numcol + 1)]
        sig  <- PARAMS[(p_NC + numcol + 2)]
        df   <- PARAMS[(p_NC + numcol + 3)]
      }
      else if (p_C > 0 && p_NC == 0){
        B_NC <- 0
        B_C  <- PARAMS[1:p_C]
        B_x  <- PARAMS[(p_C + 1):(p_C + numcol)]
        IC   <- PARAMS[(p_C + numcol + 1)]
        sig  <- PARAMS[(p_C + numcol + 2)]
        df   <- PARAMS[(p_C + numcol + 3)]
      }
      else if (p_C == 0 && p_NC == 0){
        B_NC  <- 0
        B_C   <- 0
        B_x   <- PARAMS[(p_C + 3):(p_C + 2 + numcol)]
        IC    <- PARAMS[(p_C + numcol + 3)]
        sig   <- PARAMS[(p_C + numcol + 4)]
        df    <- PARAMS[(p_C + numcol + 5)]
      }
  }
  else{
    numcol <- 0
    B_x <- 0
        if (p_C > 0 && p_NC > 0){
          B_C  <- PARAMS[1:p_C]
          B_NC <- PARAMS[(p_C+1):(p_C + p_NC)]
          IC   <- PARAMS[(p_C + p_NC + 1)]
          sig  <- PARAMS[(p_C + p_NC + 2)]
          df   <- PARAMS[(p_C + p_NC + 3)]
        }
        else if (p_NC > 0 && p_C == 0){
          B_C  <- 0
          B_NC <- PARAMS[1:p_NC]
          IC   <- PARAMS[(p_NC + 1)]
          sig  <- PARAMS[(p_NC + 2)]
          df   <- PARAMS[(p_NC + 3)]
        }
        else if (p_C > 0 && p_NC == 0){
          B_NC <- 0
          B_C  <- PARAMS[1:p_C]
          IC   <- PARAMS[(p_C + 1)]
          sig  <- PARAMS[(p_C + 2)]
          df   <- PARAMS[(p_C + 3)]
        }
        else if (p_C == 0 && p_NC == 0){
          B_NC  <- 0
          B_C   <- 0
          IC    <- PARAMS[1]
          sig   <- PARAMS[2]
          df    <- PARAMS[3]
        }
  }

  ZC1 <- y[(p_C+1):length(y)]
  ZC1 <- fBasics::vec(ZC1)
  ZC2 <- regressor.matrix(y,"not",p_C)

  if (p_C == 1){
    ZC2 <- fBasics::vec(ZC2)
  }

  if (p_C > 0){
    V <- ZC1 - ZC2 %*% B_C
  }
  else{
    V <- ZC1
  }

  U <- rev(V)
  U <- fBasics::vec(U)

  ZNC1 <- U[(p_NC + 1):length(U)]
  ZNC1 <- fBasics::vec(ZNC1)
  ZNC2 <- regressor.matrix(U,"not",p_NC)

  if(numcol > 1){
    for (i in 1:numcol){
      x[,i] <- rev(x[,i])
    }
  }
  else{
    x <- rev(x)
  }

  if(length(x) > 1){
      if (numcol > 1 ){
        x <- x[(p_NC +1):length(U),]
      }
      else{
        x <- x[(p_NC +1):length(U)]
        x <- fBasics::vec(x)
      }
  }
  else{
    x <- "not"
  }


  if (p_NC == 1){
    ZNC2 <- fBasics::vec(ZNC2)
  }

  if (length(x) > 1){
    if (p_NC > 0){
      E <- rev(ZNC1 - (ZNC2 %*% B_NC) - IC - (x %*% B_x))
    }
    else{
      E <- rev(ZNC1 - IC - (x %*% B_x))
    }
  }
  else{
    if (p_NC > 0){
      E <- rev(ZNC1 - (ZNC2 %*% B_NC) - IC)
    }
    else{
      E <- rev(ZNC1 - IC)
    }

  }

  se <- sqrt(diag(solve(optimization_results$hessian)))
  se.dist <- se[(length(se)-1):length(se)]
  se.dist <- rev(se.dist)

  return(list(coef.c = B_C, coef.nc = B_NC, coef.exo = B_x, coef.int = IC, scale = sig,df = df,residuals = E, se.dist = se.dist))
}

#' @title Asymptotic inference for the MARX function
#' @description This function allows you to calculate standard errors and confidence intervals for parameters of the MARX model.
#' @param y Data vector of time series observations.
#' @param x Matrix of data (every column represents one time series). Specify NULL or "not" if not wanted.
#' @param B_C Estimated causal parameters of the MARX.
#' @param B_NC Estimated noncausal parameters of the MARX.
#' @param B_x Estimated parameters of the exogenous variables in the MARX.
#' @param IC Estimated intercept.
#' @param sig Estimated scale parameter of the assumed underlying Student-t distribution of the residuals.
#' @param df Estimated degrees of freedom of the assumed underlying Student-t distribution of the residuals.
#' @param sig_level Significance level for the construction of inference.
#' @keywords inference
#' @keywords causal-noncausal
#' @return \item{CI.c}{Confidence intervals for causal parameters.}
#' @return \item{CI.nc}{Confidence intervals for noncausal parameters.}
#' @return \item{CI.exo}{Confidence intervals for exogenous parameters.}
#' @return \item{CI.int}{Confidence interval for intercept.}
#' @return \item{se.c}{Standard errors of causal parameters.}
#' @return \item{se.nc}{Standard errors of noncausal parameters.}
#' @return \item{se.exo}{Standard errors of exogenous parameters.}
#' @return \item{se.int}{Standard error of intercept.}
#' @author Sean Telg
#' @export
#' @examples
#' data <- sim.marx(c('t',1,1), c('t',1,1),100,0.5,0.4,0.3)
#' y <- data$y
#' x <- data$x
#' res <- marx.t(y,x,1,1)
#' inference(y,x,res$coef.c,res$coef.nc,res$coef.exo,res$coef.int,res$scale,res$df,0.05)

inference <- function(y,x,B_C,B_NC,B_x,IC,sig,df,sig_level){

  if (is.null(x)){
    x <- "not"
  }

  p_C = length(B_C)
  p_NC = length(B_NC)
  p_x = length(B_x)

  if (length(B_C) == 1 && B_C == 0){
    p_C = 0
  }

  if (length(B_NC) == 1 && B_NC == 0){
    p_NC = 0
  }

  obser = length(y) - p_C - p_NC

  regressor_C <- regressor.matrix(y,"not",p_C)
  regressand_C <- fBasics::vec(y[(p_C+1):length(y)])

  if (p_C == 1){
    regressor_C <- fBasics::vec(regressor_C)
  }

  regressor_NC  <- regressor.matrix(fBasics::vec(rev(y)),"not",p_NC)

  if (p_NC > 1){
    for(i in 1:NCOL(regressor_NC)){
      regressor_NC[,i] <- fBasics::vec(rev(regressor_NC[,i]))
    }
  }

  if (p_NC == 1){
    regressor_NC <- fBasics::vec(rev(fBasics::vec(regressor_NC)))
  }

  #if (p_NC > 1){
  #  for (i in 1:length(regressor_NC[1,])){
  #    regressor_NC[,i] <- rev(regressor_NC[,i])
  #  }
  #}

  regressand_NC <- fBasics::vec(y[1:(length(y)-p_NC)])

  ## Causal part:

  if (p_C > 0){
    if (p_NC > 0){
      U <- regressand_NC - regressor_NC %*% B_NC
    }
    else{
      U <- regressand_NC
    }

    U_regressor_C <- regressor.matrix(U,"not",p_C)
    Gam_C <- (t(U_regressor_C) %*% U_regressor_C)/obser
    Cov_C <- (df+3)/(df+1)*sig^2*solve(Gam_C)

    std_c <- matrix(data=NA,nrow=length(Cov_C[,1]), ncol=1)
    BC <- matrix(data=NA,nrow=length(Cov_C[,1]), ncol=2)

    for(i in 1:length(Cov_C[,1])){
      std_c[i] = sqrt(Cov_C[i,i]/obser)
      BC[i,] = c((B_C[i] - stats::qnorm(1-(sig_level/2))*std_c[i]),(B_C[i] + stats::qnorm(1-(sig_level/2))*std_c[i]))
    }
  }
  else if (p_C == 0){
    BC = 0
    std_c = 0
  }

  ## Noncausal part:

  if (p_NC > 0){

    if (p_C > 0){
      V <- regressand_C - regressor_C %*% B_C
    }
    else{
      V <- regressand_C
    }

    V_regressor_NC <- regressor.matrix(fBasics::vec(rev(V)),"not",p_NC)


    if (p_NC > 1){
      for(i in 1:NCOL(V_regressor_NC)){
        V_regressor_NC[,i] <- fBasics::vec(rev(V_regressor_NC[,i]))
      }
    }

    if (p_NC == 1){
      V_regressor_NC <- fBasics::vec(rev(fBasics::vec(V_regressor_NC)))
    }

    #for (i in 1:length(V_regressor_NC[1,])){
    #  V_regressor_NC[,i] <- rev(V_regressor_NC[,i])
    #}

    Gam_NC <- (1/obser)* (t(V_regressor_NC) %*% V_regressor_NC)
    Cov_NC <- (df+3)/(df+1)*sig^2*solve(Gam_NC)

    std_nc <- matrix(data=NA,nrow=length(Cov_NC[,1]), ncol=1)
    BNC <- matrix(data=NA,nrow=length(Cov_NC[,1]), ncol=2)

    for(i in 1:length(Cov_NC[,1])){
      std_nc[i] = sqrt(Cov_NC[i,i]/obser)
      BNC[i,] = c((B_NC[i] - stats::qnorm(1-(sig_level/2))*std_nc[i]),(B_NC[i] + stats::qnorm(1-(sig_level/2))*std_nc[i]))
    }
  }
  else if (p_NC == 0){
    BNC = 0
    std_nc = 0
  }

  ## Intercept

  Cov_IC <- (df+3)/(df+1)*sig^2
  std_ic <- sqrt(Cov_IC/obser)
  BIC <- c((IC - stats::qnorm(1-(sig_level/2))*std_ic),(IC + stats::qnorm(1-(sig_level/2))*std_ic))

  ## Coefficients of exogenous variables
  if (length(x) > 1){
      Gam_x <- (1/obser) * (t(x) %*% x)
      Cov_x <- (df+3)/(df+1)*sig^2*solve(Gam_x)

      std_x <- matrix(data=NA,nrow=length(Cov_x[,1]), ncol=1)
      Bx <- matrix(data=NA,nrow=length(Cov_x[,1]), ncol=2)

    for(i in 1:length(Cov_x[,1])){
      std_x[i] = sqrt(Cov_x[i,i]/obser)
      Bx[i,] = c((B_x[i] - stats::qnorm(1-(sig_level/2))*std_x[i]),(B_x[i] + stats::qnorm(1-(sig_level/2))*std_x[i]))
    }
  }
  else{
    Bx <- 0
    std_x <- 0
  }

  return(list(CI.c = BC, CI.nc = BNC, CI.exo = Bx, CI.int = BIC, se.c = std_c, se.nc = std_nc, se.exo = std_x, se.int = std_ic))
}

#' @title The value of the t-log-likelihood for MARX function
#' @description This function allows you to determine the value of the t-log-likelihood for the MARX model.
#' @param params List of parameters.
#' @param y Data vector of time series observations.
#' @param x Matrix of data (every column represents one time series). Specify NULL or "not" if not wanted.
#' @param p_C Number of lags.
#' @param p_NC Number of leads.
#' @keywords optimization
#' @return \item{neg.loglikelihood}{Minus the loglikelihood.}
#' @author Sean Telg
#' @export
#' @examples
#' data <- sim.marx(c('t',1,1), c('t',1,1),100,0.5,0.4,0.3)
#' y <- data$y
#' x <- data$x
#' p_C <- 1
#' p_NC <- 1
#' params <- c(0.5,0.4,0.3,0,1,1)
#' ll.max(params,y,x,p_C,p_NC)

ll.max <- function(params,y,x,p_C,p_NC){

  if (is.null(x)){
    x <- "not"
  }

  y <- fBasics::vec(y)
  if (length(x) > 1){
    colnum <- NCOL(x)

    if (p_C > 0 && p_NC > 0){
      BC1  <- params[1:p_C]
      BNC1 <- params[(p_C+1):(p_C + p_NC)]
      Bx1  <- params[(p_C+ p_NC + 1):(p_C + p_NC + colnum)]
      IC1  <- params[(p_C + p_NC + colnum + 1)]
      sig1 <- params[(p_C + p_NC + colnum + 2)]
      df1  <- params[(p_C + p_NC + colnum + 3)]
    }
    else if (p_NC > 0 && p_C == 0){
      BC1  <- 0
      BNC1 <- params[1:p_NC]
      Bx1  <- params[(p_NC+1):(p_NC + colnum)]
      IC1  <- params[(p_NC + colnum + 1)]
      sig1 <- params[(p_NC + colnum + 2)]
      df1  <- params[(p_NC + colnum + 3)]
    }
    else if (p_C > 0 && p_NC == 0){
      BNC1 <- 0
      BC1  <- params[1:p_C]
      Bx1  <- params[(p_C + 1): (p_C + colnum)]
      IC1  <- params[(p_C + colnum + 1)]
      sig1 <- params[(p_C + colnum + 2)]
      df1  <- params[(p_C + colnum + 3)]
    }
    else if (p_C == 0 && p_NC == 0){
      BNC1  <- 0
      BC1   <- 0
      Bx1   <- params[(1:colnum)]
      IC1   <- params[(colnum + 1)]
      sig1  <- params[(colnum + 2)]
      df1   <- params[(colnum + 3)]
    }
  }
  else{
    colnum <- 0

    if (p_C > 0 && p_NC > 0){
      BC1  <- params[1:p_C]
      BNC1 <- params[(p_C+1):(p_C + p_NC)]
      IC1  <- params[(p_C + p_NC + 1)]
      sig1 <- params[(p_C + p_NC + 2)]
      df1  <- params[(p_C + p_NC + 3)]
    }
    else if (p_NC > 0 && p_C == 0){
      BC1  <- 0
      BNC1 <- params[1:p_NC]
      IC1  <- params[(p_NC + 1)]
      sig1 <- params[(p_NC + 2)]
      df1  <- params[(p_NC + 3)]
    }
    else if (p_C > 0 && p_NC == 0){
      BNC1 <- 0
      BC1  <- params[1:p_C]
      IC1  <- params[(p_C + 1)]
      sig1 <- params[(p_C + 2)]
      df1  <- params[(p_C + 3)]
    }
    else if (p_C == 0 && p_NC == 0){
      BNC1  <- 0
      BC1   <- 0
      IC1   <- params[1]
      sig1  <- params[2]
      df1   <- params[3]
    }

  }

  ZC1 <- y[(p_C+1):length(y)]
  ZC1 <- fBasics::vec(ZC1)
  ZC2 <- regressor.matrix(y,"not",p_C)

  if (p_C == 1){
    ZC2 <- fBasics::vec(ZC2)
  }

  if (p_C > 0){
    V <- ZC1 - (ZC2 %*% BC1)
  }
  else{
    V <- ZC1
  }

  U <- rev(V)
  U <- fBasics::vec(U)

  ZNC1 <- U[(p_NC + 1):length(U)]
  ZNC1 <- fBasics::vec(ZNC1)
  ZNC2 <- regressor.matrix(U,"not",p_NC)

  if(colnum > 1){
    for (i in 1:colnum){
      x[,i] <- rev(x[,i])
    }
  }
  else{
    x <- rev(x)
  }

  if (length(x) > 1){
    if (colnum > 1){
      x <- x[(p_NC +1):length(U),]
    }
    else{
      x <- x[(p_NC + 1):length(U)]
      x <- fBasics::vec(x)
    }
  }
  else{
    x = "not"
  }

  if (p_NC == 1){
    ZNC2 <- fBasics::vec(ZNC2)
  }

  if (length(x) > 1){
    if (p_NC > 0){
      E <- rev(ZNC1 - (ZNC2 %*% BNC1) - IC1 - (x %*% Bx1))
    }
    else{
      E <- rev(ZNC1 - IC1 - (x %*% Bx1))
    }
  }
  else{
    if (p_NC > 0){
      E <- rev(ZNC1 - (ZNC2 %*% BNC1) - IC1)
    }
    else{
      E <- rev(ZNC1 - IC1)
    }
  }

  n <- length(E)

  loglik_eval <- -(n*lgamma((df1+1)/2) - n*log(sqrt(df1*pi*sig1^2)) - n*lgamma(df1/2) - ((df1+1)/2)*log(1+(E/sig1)^2/df1) %*% matlab::ones(n,1))

  return(neg.loglikelihood = loglik_eval)
}

#' @title The lag-lead model selection for MARX function
#' @description This function allows you to determine the MARX model (for p = r + s) that maximizes the t-log-likelihood.
#' @param y Data vector of time series observations.
#' @param x Matrix of data (every column represents one time series). Specify NULL or "not" if not wanted.
#' @param p_pseudo Number of autoregressive terms to be included in the pseudo-causal model.
#' @keywords selection
#' @keywords causal-noncausal
#' @return \item{p.C}{The number of lags selected.}
#' @return \item{p.NC}{The number of leads selected.}
#' @return \item{loglikelihood}{The value of the loglikelihood for all models with p = r + s.}
#' @author Sean Telg
#' @export
#' @examples
#' data <- sim.marx(c('t',3,1), c('t',3,1),100,0.5,0.4,0.3)
#' selection.lag.lead(data$y,data$x,2)

selection.lag.lead <- function(y,x,p_pseudo){

  if (is.null(x)){
    x <- "not"
  }

  P_C <- seq(length=(p_pseudo+1), from=0, by=1)
  P_C <- fBasics::vec(P_C)
  P_NC <- rev(P_C)
  P_NC <- fBasics::vec(P_NC)

  n <- length(y) - p_pseudo

  loglik <- c()

  for (i in 1:(p_pseudo+1)){

    marx.t_results <- marx.t(y,x,P_C[i],P_NC[i]);
    sig <- marx.t_results[[5]]
    df  <- marx.t_results[[6]]
    E   <- marx.t_results[[7]]

    loglik[i] <- (n*lgamma((df+1)/2) - n*log(sqrt(df*pi*sig^2)) - n*lgamma(df/2) - ((df+1)/2)*log(1+(E/sig)^2/df) %*% matlab::ones(n,1))

  }

  maxloglik <- which.max(loglik)

  P = cbind(P_C,P_NC)

  P = fBasics::vec(P[maxloglik,])

  p_C  <- P[1]
  p_NC <- P[2]

  return(list(p.C = p_C, p.NC = p_NC,loglikelihood = rev(loglik)))
}

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

selection.lag <- function(y,x,p_max){

  if (is.null(x)){
    x <- "not"
  }

  bic_results <- bic(y,x,p_max)
  aic_results <- aic(y,x,p_max)
  hq_results <- hq(y,x,p_max)

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

  return(list(bic = bic_vec, aic = aic_vec, hq = hq_vec))
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

bic <- function(y,x,p_max){

  if (is.null(x)){
    x <- "not"
  }

  y <- fBasics::vec(y)
  # y <- y - mean(y)
  # n <- length(y) - p_max

  if (length(x) > 1){
    numcol <- NCOL(x)
  }
  else{
    numcol = 0
  }

  crit <- matrix(data=NA, nrow=(p_max+1), ncol=1)

  for (p in 0:p_max){

    arx.ls_results <- arx.ls(fBasics::vec(y),x,p)
    n <- length(arx.ls_results[[5]])
    Cov <- arx.ls_results[[6]]
    crit[(p+1)] <- -2*Cov/n + ((log(n))/n)*(p+1+numcol)
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

aic <- function(y,x,p_max){

  if (is.null(x)){
    x <- "not"
  }

  y <- fBasics::vec(y)
  # y <- y - mean(y)
  #n <- length(y) - p_max

  if (length(x) > 1){
    numcol <- NCOL(x)
  }
  else{
    numcol = 0
  }

  crit <- matrix(data=NA, nrow=(p_max+1), ncol=1)

  for (p in 0:p_max){

    arx.ls_results <- arx.ls(fBasics::vec(y),x,p)
    n <- length(arx.ls_results[[5]])
    Cov <- arx.ls_results[[6]]
    crit[(p+1)] <- -2*Cov/n + (2/n)*(p+1+numcol)
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

hq <- function(y,x,p_max){

  if (is.null(x)){
    x <- "not"
  }

  y <- fBasics::vec(y)
  # y <- y - mean(y)
  # n <- length(y) - p_max

  if (length(x) > 1){
    numcol <- NCOL(x)
  }
  else{
    numcol = 0
  }

  crit <- matrix(data=NA, nrow=(p_max+1), ncol=1)

  for (p in 0:p_max){

    arx.ls_results <- arx.ls(fBasics::vec(y),x,p)
    n <- length(arx.ls_results[[5]])
    Cov <- arx.ls_results[[6]]
    crit[(p+1)] <- -2*Cov/n + ((2*log(log(n)))/n)*(p+1+numcol)
  }

  p_hq <- which.min(crit) - 1

  crit <- t(crit)
  colnames(crit) <- paste('p =', 0:p_max)

  return(list(p = p_hq, values = crit))

}

#' @title The simulation of MARX processes
#' @description This function allows you to simulate MARX processes based on different underlying distribution.
#' @param dist.eps vector containing the error distribution and its parameters (options: t, normal, stable).
#' @param dist.x   vector containing the distribution of x and its parameters (options: t, normal, stable). Specify NULL or "not" if not wanted.
#' @param obs      Number of observations for simulated process.
#' @param c_par    vector of causal parameters.
#' @param nc_par   vector of noncausal parameters.
#' @param exo_par  Parameter of the exogenous variable.
#' @keywords simulation
#' @return \item{y}{Simulated data y.}
#' @return \item{x}{Simulated data x (exogenous variable).}
#' @author Sean Telg
#' @export
#' @examples
#' dist.eps <- c('t',1,1) ## t-distributed errors with 1 degree of freedom and scale parameter 1
#' dist.x   <- c('normal',0,1) ## standard normally distributed x variable
#' obs <- 100
#' c_par <- c(0.2,0.4)
#' nc_par <- 0.8
#' exo_par <- 0.5
#' sim.marx(dist.eps,dist.x,obs,c_par,nc_par,exo_par) ## Simulates a MARX(2,1,1) process

sim.marx <- function(dist.eps, dist.x, obs, c_par, nc_par, exo_par){

  if (is.null(dist.x)){
    exo_par = 0
    dist.x = 0
  }

  if (length(dist.x) == 1){
    exo_par = 0
  }

  M <- 5*obs + 300
  M.star <- 500   #0.05*M

  if (dist.eps[1] == 't'){
    eps <- as.numeric(dist.eps[3])*(stats::rt((3*(obs+200)), df = as.numeric(dist.eps[2]), ncp = 0))
  } else if (dist.eps[1] == 'cauchy'){
    eps <- stats::rcauchy((3*(obs+200)), location= as.numeric(dist.eps[2]), scale=as.numeric(dist.eps[3]))
  } else if (dist.eps[1] == 'normal') {
    eps <- stats::rnorm((3*(obs+200)), mean = as.numeric(dist.eps[2]), sd=as.numeric(dist.eps[3]))
  } else if (dist.eps[1] == 'stable') {
    eps <- stabledist::rstable((3*(obs+200)), alpha = as.numeric(dist.eps[2]), beta=as.numeric(dist.eps[3]), gamma=as.numeric(dist.eps[4]), delta=as.numeric(dist.eps[5]))
  }

  if (dist.x[1] == 't'){
    x <- as.numeric(dist.x[3])*(stats::rt((3*(M+200)), df = as.numeric(dist.x[2]), ncp = 0))
  } else if (dist.x[1] == 'cauchy'){
    x <- stats::rcauchy((3*(M+200)), location= as.numeric(dist.x[2]), scale=as.numeric(dist.x[3]))
  } else if (dist.x[1] == 'normal') {
    x <- stats::rnorm((3*(M+200)), mean = as.numeric(dist.x[2]), sd=as.numeric(dist.x[3]))
  } else if (dist.x[1] == 'stable') {
    x <- stabledist::rstable((3*(M+200)), alpha = as.numeric(dist.x[2]), beta=as.numeric(dist.x[3]), gamma=as.numeric(dist.x[4]), delta=as.numeric(dist.x[5]))
  }

  eps <- eps[(2*(obs+200)):(3*(obs+200))]

  r <- length(c_par)
  s <- length(nc_par)

  if (length(dist.x) != 1){
    BC <- c(c_par, rep(0,M.star - length(c_par)))
    Phi.matrix <- matrix(data=0, nrow=M.star, ncol=M.star)
    Psi <- c()
    Psi[1] <- 1

    for (i in 1:M.star){
      for (j in 1:i){
        Phi.matrix[i,j] = Psi[i-j+1]*BC[j]
      }
      Psi[(i+1)] <- sum(Phi.matrix[i,])
    }


    ## Generating linear filter over x_t "similar to MA representation for AR"
    MA_rep <- matrix(data=NA, nrow=(obs+200), ncol=M.star)
    for (i in 1:(obs+200)){
      for (j in 1:M.star){
        MA_rep[i,j] <- Psi[j]*x[M-(j-1)-(i-1)]
      }
    }
    exo_var <- c()
    exo_var <- MA_rep[,1]
    MA_rep <- rowSums(MA_rep)
    #MA_rep <- rev(MA_rep)
  }
  else{
    MA_rep <- matrix(data=0, nrow=(obs+200),ncol=M.star)
    exo_var <- c()
  }


  u <- matrix(data=NA, nrow=(obs+100+s), ncol=1)
  U.mat <- u[(201+obs):(200+obs+s)] <- rep(1,s)
  for (i in (obs+200):1){
    if (length(nc_par) == 1){
      AR <- nc_par * U.mat
    }
    else{
      AR <- sum(t(nc_par*U.mat))
    }
    if (i == (obs+200)){AR <- 0}

    u[i] <- AR + eps[i]

    if (s == 1){
      U.mat <- u[i]
    }
    else{
      U.mat <- c(u[i], t(U.mat[1:(length(U.mat)-1)]))
      U.mat <- t(U.mat)
    }
  }
  y <- matrix(data=0, nrow=(obs+r+200), ncol=1)
  Y <- y[1:r] <- rep(0,r)

  for (i in 1:(obs+200)){
    if (r == 1){
      AR2 <- c_par * Y
    }
    else{
      AR2 <- sum(t(c_par*Y))
    }
    if (i == 1){AR2 <- 0}

    y[i+r] <- AR2 + u[i] + exo_par * MA_rep[i+r]

    if (r == 1){
      Y <- y[i+r]
    }
    else{
      Y <- c(y[i+r], t(Y[1:(length(Y)-1)]))
      Y <- t(Y)
    }
  }
  u <- u[1:obs]
  y <- y[(101+r):(obs+100+r)]
  exo_var <- exo_var[(101+r):(obs+100+r)]
  return(list(y = y, x = exo_var))
}

#' @title Companion form function
#' @description   This function allows you to compute a companion form matrix in order to check the stability of causal and noncausal part of the ARX model.
#' @param pol     Coefficient vector. If polynomial is 1 - ax - bx^2, coefficient vector is c(a, b).
#' @keywords stability, stationarity
#' @return \item{C}{Companion matrix C.}
#' @author Sean Telg
#' @export
#' @examples
#' pol <- c(0.3,0.4)
#' C <- companion.form(pol)

companion.form <- function(pol){
  r <- length(pol)
  b <- cbind(diag((r-1)),rep(0,(r-1)))
  C <- rbind(t(pol),b)

  return(C)
}

#' @title Coefficients of the moving average representation function
#' @description   This function allows you to invert a polynomial (either the causal or the noncausal one) and output the corresponding coefficients of the moving average representation.
#' @param pol     Coefficient vector. If polynomial is 1 - ax - bx^2, coefficient vector is c(a, b).
#' @param M       Truncation value M (how many MA coefficients should be computed?).
#' @keywords stability, stationarity
#' @return \item{psi}{Vector containing coefficients of the moving average representation.}
#' @author Sean Telg
#' @export
#' @examples
#' pol <- c(0.3,0.4)
#' psi <- companion.form(pol)

compute.MA <- function(pol,M){
  r <- length(pol)
  str <- c(rep(0,(r-1)),1)
  psi <- c(1)

  for (j in 1:M){
    psi[j+1] <- t(pol) %*% rev(str)

    if (r > 1){
      str <- c(str[2:length(str)],psi[j+1])
    }
    else{
      str <- psi[j+1]
    }
  }
  return(psi)
}

#' @title Forecasting function for the MARX model
#' @description   This function allows you to forecast with the mixed causal-noncausal model with possibly exogenous regressors.
#' @param y       Data vector y.
#' @param X       (optional) Matrix with data (column represent a series).
#' @param p_C     Number of lags (causal order).
#' @param p_NC    Number of leads (noncausal order).
#' @param X.for   (optional) Matrix with forecasted values for X (column represents series).
#' @param h       Forecast horizon h.
#' @param M       (optional) Truncation value M for MA representation. Default value: 50.
#' @param N       (optional) Number of simulations to forecast noncausal component. Default: 10,000.
#' @keywords forecasting
#' @return \item{y.for}{Vector containing forecasted values for y.}
#' @author Sean Telg
#' @export
#' @examples
#' ## Forecasting MAR(0,1) model 4-periods ahead for lnbev (from dataset)
#' data <- MARX::dataset[,2]
#' y.for <- forecast.marx(y=data, p_C=0, p_NC=1, h=4, M=50, N=1000)

forecast.marx <- function(y,X,p_C,p_NC,X.for,h,M,N){

  set.seed(9999)
  if (missing(X) == TRUE){
    X = NULL
  }

  if (missing(N) == TRUE){
    N = 10000
  }

  object <- mixed(y,X,p_C,p_NC)
  obs <- length(y)

  ## Check whether there are exogenous variables and whether truncation M is known

  if (missing(X.for) == TRUE && missing(M) == TRUE){
    X.for = NULL
    M = 50
  }
  else if(missing(X.for) == TRUE && missing(M) == FALSE){
    X.for = NULL
    M = M
  }
  else if(missing(X.for) == FALSE && missing(M) == TRUE){
    if (NCOL(X.for) == 1){
      if(is.null(X.for) == TRUE){
        M = 50
      }
      else{
        M = length(X.for)
      }
    }
    else{
      M = length(X.for[,1])
    }
  }
  else if(missing(X.for) == FALSE && missing(M) == FALSE){
    if (NCOL(X.for) == 1){
      if(is.null(X.for) == TRUE){
        M = M
      }
      else{
        M = min(length(X.for), M)
      }
    }
    else{
      M = min(length(X.for[,1]),M)
    }
  }

  coef.caus <- c()
  if (object$order[1] == 0){
    r = 1
    coef.caus <- object$coefficients[(r+1)]
  }
  else{
    r = object$order[1]
    coef.caus <- object$coefficients[2:(r+1)]
  }

  coef.noncaus <- c()
  if (object$order[2] == 0){
    s = 1
    coef.noncaus <- object$coefficients[(r+1+s)]
  }
  else{
    s = object$order[2]
    coef.noncaus <- object$coefficients[(r+2):(r+1+s)]
  }

  coef.exo <- c()
  if (object$order[3] == 0){
    q = 1
    coef.exo <- object$coefficients[(r+1+s+q)]
  }
  else{
    q = object$order[3]
    coef.exo <- object$coefficients[(r+1+s+1):(r+s+1+q)]
  }

  ## Simulate future epsilon and use forecasted X
  hve <- c()
  hve2 <- matrix(data=0, nrow=N,ncol=h)

  for (iter in 1:N){

    eps.sim <- object$coefficients["scale",]*stats::rt(M,object$coefficients["df",])

    z2 <- c()
    for (i in 1:M){
      if(is.null(X.for) == TRUE){
        z2[i] <- eps.sim[i]
      }
      else{
        if(NCOL(X.for) > 1){
          z2[i] <- eps.sim[i] +  coef.exo %*% t(X.for[i,])
        }
        else{
          z2[i] <- eps.sim[i] + coef.exo * X.for[i]
        }
      }
    }

    ## Compute filtered values u = phi(L)y and moving average values
    phi <- c(1,coef.caus)

    u <- c()
    for (i in (r+1):obs){
      u[i] <- phi %*% y[i:(i-r)]
    }
    w <- c(u[(obs-s+1):obs],z2)

    C <- matrix(data=0, nrow=(M+s), ncol=(M+s))
    C[1,] <- compute.MA(coef.noncaus,(M+s-1))

    if (s > 1){
      for (i in 2:s){
        C[i,] <- c(0, C[(i-1),1:(length(C[(i-1),])-1)])
      }
    }

    for (i in (s+1):(M+s)){
      C[i,] <- c(rep(0,(i-1)),1,rep(0,(M+s-i)))
    }

    D = solve(C)

    e <- D %*% w

    h1 <- c()

    for (i in 1:s){
      h1[i] <- metRology::dt.scaled(e[i], df=object$coefficients["df",], sd=object$coefficients["scale",])
    }

    hve[iter] = prod(h1)

    for (j in 1:h){
      mov.av <-  C[1,1:(M-j+1)] %*% z2[j:M]
      hve2[iter,j] <- mov.av * hve[iter]

    }
  }

  y.star <- y[(obs-r+1):obs]
  y.for <- c()
  exp <- c()

  for (j in 1:h){
    exp[j] = ((1/N)*sum(hve2[,j]))/((1/N)*sum(hve))

    if(length(coef.caus) == 1){
      y.for[j] <-  object$coefficients[1]/(1-sum(coef.noncaus)) + coef.caus * y.star + exp[j]
    }
    else{
      y.for[j] <-  object$coefficients[1]/(1-sum(coef.noncaus)) + t(coef.caus) %*% y.star + exp[j]
    }

    y.star <- c(y.for[j], y.star[1:(length(y.star)-1)])
  }

  return(y.for)
}

mixed.combine <- function(y,x,p_C,p_NC){

  q <- NCOL(x)

  if (is.null(x) == TRUE){
    x <- "not"
    q <- 0
  }

  sig_level = 0.05
  est <- marx.t(y,x,p_C,p_NC)
  residuals <- est$residuals
  fitted.values <- y[(1+p_C):(length(y)-p_NC)] - residuals

  inf <- inference(y,x,est$coef.c,est$coef.nc,est$coef.exo,est$coef.int,est$scale,est$df,sig_level)

  coefficients <- c(est$coef.int,est$coef.c,est$coef.nc,est$coef.exo, est$df, est$scale)
  coefficients <- fBasics::vec(coefficients)

  if(length(x) > 1){
    if (p_C > 0 && p_NC > 0){
      rownames(coefficients) <- c('int', paste('lag', 1:p_C), paste('lead', 1:p_NC), paste('exo', 1:NCOL(x)), 'df', 'scale')
    }
    else if (p_C > 0 && p_NC == 0){
      rownames(coefficients) <- c('int', paste('lag', 1:p_C), 'lead', paste('exo', 1:NCOL(x)), 'df', 'scale')
    }
    else if (p_C == 0 && p_NC > 0){
      rownames(coefficients) <- c('int', 'lag', paste('lead', 1:p_NC), paste('exo', 1:NCOL(x)), 'df', 'scale')
    }
    else if (p_C == 0 && p_NC == 0){
      rownames(coefficients) <- c('int', 'lag', 'lead', paste('exo', 1:NCOL(x)), 'df', 'scale')
    }
  }
  else{
    if (p_C > 0 && p_NC > 0){
      rownames(coefficients) <- c('int', paste('lag', 1:p_C), paste('lead', 1:p_NC), 'exo', 'df', 'scale')
    }
    else if (p_C > 0 && p_NC == 0){
      rownames(coefficients) <- c('int', paste('lag', 1:p_C), 'lead', 'exo', 'df', 'scale')
    }
    else if (p_C == 0 && p_NC > 0){
      rownames(coefficients) <- c('int', 'lag', paste('lead', 1:p_NC), 'exo', 'df', 'scale')
    }
    else if (p_C == 0 && p_NC == 0){
      rownames(coefficients) <- c('int', 'lag', 'lead', 'exo', 'df', 'scale')
    }
  }

  se <- c(inf$se.int, inf$se.c, inf$se.nc, inf$se.exo, est$se.dist[1], est$se.dist[2])

  if (length(x) > 1){
    degree <- NROW(est$residuals) - (p_C + p_NC + NCOL(x) + 1)
  }
  else{
    degree <- NROW(est$residuals) - (p_C + p_NC + 1)
  }

  return(list(coefficients = coefficients, se = se, df.residual = degree, residuals= residuals, fitted.values = fitted.values, order=c(p_C,p_NC,q)))
}


#' @title  The pseudo-causal model function
#' @description This function allows you to estimate pseudo-causal ARX models by OLS (compatible with most functions in lm() class).
#' @aliases   pseudo
#'            pseudo.default
#'            print.pseudo
#'            summary.pseudo
#' @param y      Data vector of time series observations.
#' @param x      Matrix of data (every column represents one time series). Specify NULL or "not" if not wanted.
#' @param p      Number of lags to be included.
#' @param object An object of the class "pseudo"
#' @param ...    Other arguments
#' @keywords pseudo-causal
#' @keywords estimation
#' @examples
#' data <- sim.marx(c('t',1,1), c('t',1,1),100,0.5,0.4,0.3)
#' object <- pseudo(data$y, data$x, 2)
#' class(object) <- "pseudo"
#' summary(object)
#' @rdname pseudo
#' @export pseudo
pseudo <- function(y,x,p){ UseMethod("pseudo")}
#' @return An object of class \code{"pseudo"} is a list containing the following components:
#' @return \item{coefficients}{Vector of estimated coefficients.}
#' @return \item{coef.auto}{Vector of estimated autoregressive parameters.}
#' @return \item{coef.exo}{Vector of estimated exogenous parameters.}
#' @return \item{mse}{Mean squared error.}
#' @return \item{residuals}{Residuals.}
#' @return \item{loglikelihood}{Value of the loglikelihood.}
#' @return \item{fitted.values}{Fitted values.}
#' @return \item{df}{Degrees of freedom.}
#' @return \item{vcov}{Variance-covariance matrix of residuals.}
#'
#' @rdname pseudo
#' @method pseudo default
#' @S3method  pseudo default
pseudo.default <- function(y,x,p){

  if (is.null(x)){
    x <- "not"
  }

  est <- arx.ls(y,x,p)
  est$call <- match.call();
  class(est) <- "pseudo"
  est
}
#' @rdname pseudo
#' @method print pseudo
#' @S3method print pseudo
print.pseudo <- function(x,...){

  cat("Call:\n")
  print(x$call)
  cat("\n Coefficients:\n")
  print(x$coefficients)
}
#' @rdname pseudo
#' @method summary pseudo
#' @S3method summary pseudo
summary.pseudo <- function(object,...){

  se <- sqrt(diag(object$vcov))
  tval <- stats::coef(object) / se

  TAB <- cbind(Estimate = stats::coef(object),
               StdErr = se,
               t.value = tval,
               p.value = 2*stats::pt(-abs(tval),df=object$df))

  colnames(TAB) <- c("Estimate", "Std.Err", "t value", "p value")

  res <- list(call=object$call, coefficients = TAB)

  class(res) <- "summary.pseudo"
  res
}


pseudo <- function(y,x,p){ UseMethod("pseudo")}

pseudo.default <- function(y,x,p){

  if (is.null(x)){
    x <- "not"
  }

  est <- arx.ls(y,x,p)
  est$call <- match.call();

  class(est) <- "pseudo"
  est
}

print.pseudo <- function(x,...){

  cat("Call:\n")
  print(x$call)
  cat("\n Coefficients:\n")
  print(x$coefficients)

}

summary.pseudo <- function(object,...){

  se <- sqrt(diag(object$vcov))
  tval <- stats::coef(object) / se

  TAB <- cbind(Estimate = stats::coef(object),
               StdErr = se,
               t.value = tval,
               p.value = 2*stats::pt(-abs(tval),df=object$df))

  colnames(TAB) <- c("Estimate", "Std.Err", "t value", "p value")

  res <- list(call=object$call, coefficients = TAB)

  class(res) <- "summary.pseudo"
  res
}

#' @title The MARX estimation function
#' @description This function allows you to estimate mixed causal-noncausal MARX models by t-MLE (compatible with most functions in lm() class).
#' @aliases    mixed
#'             mixed.default
#'             print.mixed
#'             summary.mixed
#' @param y      Data vector of time series observations.
#' @param x      Matrix of data (every column represents one time series). Specify NULL or "not" if not wanted.
#' @param p_C    Number of lags to be included.
#' @param p_NC   Number of leads to be included.
#' @param object An object of the class "mixed".
#' @param ...    Other parameters.
#' @keywords causal-noncausal
#' @keywords estimation
#' @examples
#' data <- sim.marx(c('t',1,1), c('t',1,1),100,0.5,0.4,0.3)
#' object <- mixed(data$y, data$x, 1, 1)
#' class(object) <- "mixed"
#' summary(object)
#' @rdname mixed
#' @export mixed
mixed <- function(y,x,p_C,p_NC){ UseMethod("mixed")}
#' @return An object of class \code{"mixed"} is a list containing the following components:
#' @return \item{coefficients}{Vector of estimated coefficients.}
#' @return \item{se}{Standard errors of estimated coefficients.}
#' @return \item{df.residual}{Degrees of freedom residuals.}
#' @return \item{residuals}{Residuals.}
#' @return \item{fitted.values}{Fitted values.}
#' @return \item{order}{Vector containing (r,s,q), i.e. causal order r, noncausal order s, number of exogenous regressors q.}
#' @rdname mixed
#' @method mixed default
#' @S3method mixed default
mixed.default <- function(y,x,p_C,p_NC){

  est <- mixed.combine(y,x,p_C,p_NC);
  est$call <- match.call();
  class(est) <- "mixed"
  est
}
#' @rdname mixed
#' @method print mixed
#' @s3method print mixed
print.mixed <- function(x,...){

  cat("Call:\n")
  print(x$call)
  cat("\n Coefficients:\n")
  print(x$coefficients)
}
#' @rdname mixed
#' @method summary mixed
#' @S3method summary mixed
summary.mixed <- function(object,...){

  tval <- stats::coef(object) / object$se

  TAB <- cbind(Estimate = stats::coef(object),
               StdErr = object$se,
               t.value = tval,
               p.value = 2*stats::pt(-abs(tval),df=object$df.residual))

  colnames(TAB) <- c("Estimate", "Std.Err", "t value", "p value")

  res <- list(call=object$call, coefficients = TAB)

  class(res) <- "summary.mixed"
  res
}

mixed <- function(y,x,p_C,p_NC){ UseMethod("mixed")}

mixed.default <- function(y,x,p_C,p_NC){

  est <- mixed.combine(y,x,p_C,p_NC);
  est$call <- match.call();
  class(est) <- "mixed"
  est
}

print.mixed <- function(x,...){

  cat("Call:\n")
  print(x$call)
  cat("\n Coefficients:\n")
  print(x$coefficients)

}

summary.mixed <- function(object,...){

  tval <- stats::coef(object) / object$se

  TAB <- cbind(Estimate = stats::coef(object),
               StdErr = object$se,
               t.value = tval,
               p.value = 2*stats::pt(-abs(tval),df=object$df.residual))

  colnames(TAB) <- c("Estimate", "Std.Err", "t value", "p value")

  res <- list(call=object$call, coefficients = TAB)

  class(res) <- "summary.mixed"
  res
}
