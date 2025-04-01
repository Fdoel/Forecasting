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
regressor.matrix_T <- function(y, x, p, c) {
  # Handle NULL x case
  if (missing(x) || is.null(x)) {
    x <- "not"
  }
  
  y <- fBasics::vec(y)
  n <- length(y)
  
  if (p == 1) {
    k <- 1
  } else {
    k <- NCOL(y)
  }
  
  # Create Z matrix
  if (p > 0) {
    Z <- matlab::zeros(n, k * p)
    for (i in 1:p) {
      Z[(1 + i):n, ((i - 1) * k + 1):(i * k)] <- y[1:(n - i)]
    }
    Z <- Z[(1 + p):n, ]
  } else {
    Z <- matrix(, nrow = n, ncol = 0)
  }
  
  if (identical(x, "not")) {
  } else if (NCOL(x) == 1) {
    Z <- cbind(Z, x[(1 + p):n])
  } else if (NCOL(x) > 1) {
    Z <- cbind(Z, x[(1 + p):n, ])
  }
  
  # Create thresholded ZT matrix
  if(identical(x, "not")) {
    if(p == 1) {
      nT <- length(Z)
      Z_c <- Z
    } else if(p ==0) {
      nT <- length(Z)
      Z_c <- 0
    } else {
      nT <- nrow(Z)
      Z_c <- Z[,1:p]
    }
  } else {
    nT <- nrow(Z)
    Z_c <- Z[,1:p]
  }
  if(p != 0) {
    Z_c <- cbind(Z_c, Z_c)
    if (!identical(x, "not")) {
      Z_x <- Z[,(p + 1):ncol(Z)]
      Z_x <- cbind(Z_x, Z_x)
      mX <- ncol(Z_x)
      for(i in 1:nT) {
        if(Z_c[i, 1] > c) {
          Z_x[i, 1:(mX/2)] <- 0
        } else {
          Z_x[i, (mX/2 + 1):mX] <- 0
        }
      }
    }
    mC <- ncol(Z_c)
    for(i in 1:nT) {
      if(Z_c[i, 1] > c) {
        Z_c[i, 1:(mC/2)] <- 0
      } else {
        Z_c[i, (mC/2 + 1):mC] <- 0
      }
    }
    if (!identical(x, "not")) {
      ZT <- cbind(Z_c, Z_x)
    } else {
      ZT <- Z_c
    }
  } else {
    ZT <- Z
  }
  return(matrix = ZT)
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

arx.ls_T <- function(y,x,p,c){
  
  if (is.null(x)){
    x <- "not"
  }
  
  n <- length(y) - p
  
  Y <- y[(p+1):length(y)]
  int <- rep(1,(length(y)-p))
  ZT <- regressor.matrix_T(y,x,p,c)
  ZT <- cbind(int,ZT)
  
  df <- nrow(ZT) - NCOL(ZT)
  
  B <- solve(t(ZT) %*% ZT) %*% (t(ZT) %*% Y)
  
  if (p > 0){
    if (length(x) > 1){
      rownames(B) <- c('int', paste('lag', 1:(2*p)), paste('exo', 1:(2*NCOL(x))))
    }
    else{
      rownames(B) <- c('int', paste('lag', 1:(2*p)))
    }
  }
  else{
    if (length(x) > 1){
      rownames(B) <- c('int', paste('exo', 1:(2*NCOL(x))))
    }
    else{
      rownames(B) <- 'int'
    }
  }
  
  FV <- ZT %*% B
  U <- Y - FV
  
  sig <- (t(U) %*% U)
  sig <- as.numeric(sig)
  
  Cov <- (1/n)*sig
  Cov <- as.numeric(Cov)
  
  sigma2 <- sum((Y - ZT %*% B)^2)/df
  qz <- qr(ZT)
  vcov <- sigma2*chol2inv(qz$qr)
  colnames(vcov) <- rownames(vcov) <- colnames(ZT)
  
  Loglik <- -(n/2)*(1 + log(2*pi)+log(Cov))
  
  if (p == 0){
    B_auto <- 0
  }
  else{
    B_auto <- B[2:(2*p+1)]
  }
  
  if (length(x) > 1){
    B_x <- B[(2*p+2):length(B)]
  }
  else{
    B_x <- 0
  }
  
  return(list(coefficients = B, coef.auto = B_auto, coef.exo = B_x, mse = Cov, residuals = U, loglikelihood = Loglik, fitted.values = FV, df = df,vcov=vcov))
}

# Dit verwacht params in een bepaalde vorm die nu niet zo uit ARX_T komt
ll.MART.Z <- function(params,y,x,p_C,p_NC,c){
  p_CT <- p_C*2
  p_NCT <- p_NC*2
  if (is.null(x)){
    x <- "not"
  }
  
  y <- fBasics::vec(y)
  if (length(x) > 1){
    colnum <- NCOL(as.matrix(x))
    colnumT <- colnum*2
    
    if (p_C > 0 && p_NC > 0){
      BC1  <- params[1:p_CT]
      BNC1 <- params[(p_CT+1):(p_CT + p_NCT)]
      Bx1  <- params[((p_CT + p_NCT)+ 1):(p_CT + p_NCT + colnumT)]
      IC1  <- params[(p_CT + p_NCT + colnumT + 1)]
      sig1 <- params[(p_CT + p_NCT + colnumT + 2)]
      df1  <- params[(p_CT + p_NCT + colnumT + 3)]
    } else if (p_NC > 0 && p_C == 0){
      BC1  <- 0
      BNC1 <- params[1:(p_NCT)]
      Bx1  <- params[((p_NCT)+1):((p_NCT) + colnumT)]
      IC1  <- params[(p_NCT + colnumT + 1)]
      sig1 <- params[(p_NCT + colnumT + 2)]
      df1  <- params[(p_NCT + colnumT + 3)]
    } else if (p_C > 0 && p_NC == 0){
      BNC1 <- 0
      BC1  <- params[1:(p_CT)]
      Bx1  <- params[(p_CT + 1):(p_CT + colnumT)]
      IC1  <- params[(p_CT + colnumT + 1)]
      sig1 <- params[(p_CT + colnumT + 2)]
      df1  <- params[(p_CT + colnumT + 3)]
    } else if (p_C == 0 && p_NC == 0){
      BNC1  <- 0
      BC1   <- 0
      Bx1   <- params[(1:(colnumT))]
      IC1   <- params[(colnumT + 1)]
      sig1  <- params[(colnumT + 2)]
      df1   <- params[(colnumT + 3)]
    }
  } else{
    colnum <- 0
    colnumT <- 0
    if (p_C > 0 && p_NC > 0){
      BC1  <- params[1:(p_CT)]
      BNC1 <- params[((p_CT)+1):(p_CT + p_NCT)]
      IC1  <- params[(p_CT + p_NCT + 1)]
      sig1 <- params[(p_CT + p_NCT + 2)]
      df1  <- params[(p_CT + p_NCT + 3)]
    } else if (p_NC > 0 && p_C == 0){
      BC1  <- 0
      BNC1 <- params[1:(p_NCT)]
      IC1  <- params[(p_NCT + 1)]
      sig1 <- params[(p_NCT + 2)]
      df1  <- params[(p_NCT + 3)]
    } else if (p_C > 0 && p_NC == 0){
      BNC1 <- 0
      BC1  <- params[1:(p_CT)]
      IC1  <- params[(p_CT + 1)]
      sig1 <- params[(p_CT + 2)]
      df1  <- params[(p_CT + 3)]
    } else if (p_C == 0 && p_NC == 0){
      BNC1  <- 0
      BC1   <- 0
      IC1   <- params[1]
      sig1  <- params[2]
      df1   <- params[3]
    }
  }
  
  ZC1 <- y[(p_C+1):length(y)]
  ZC1 <- fBasics::vec(ZC1)
  ZC2 <- regressor.matrix_T(y,"not",p_C, c)
  
  if (p_C > 0){
    V <- ZC1 - (ZC2 %*% BC1)
  } else{
    V <- ZC1
  }
  U <- rev(V)
  U <- fBasics::vec(U)
  
  ZNC1 <- U[(p_NC + 1):length(U)]
  ZNC1 <- fBasics::vec(ZNC1)
  ZNC2 <- regressor.matrix_T(U,"not",p_NC, c)
  if((colnumT) > 1){
    ZX <- regressor.matrix_T(y, x, 1, c)[, 3:(2+2*ncol(as.matrix(x)))]
    for (i in 1:(colnumT)) {
      ZX[,i] <- rev(ZX[,i])
    }
  }
  if (length(x) > 1){
    if ((colnumT) > 1){
      ZX <- ZX[(p_NC +1):length(U),]
    } else {
      ZX <- ZX[(p_NC + 1):length(U)]
    }
  } else {
    x = "not"
  }
  if (length(x) > 1){
    if (p_NC > 0){
      E <- rev(ZNC1 - (ZNC2 %*% BNC1) - IC1 - (ZX %*% Bx1))
    } else{
      E <- rev(ZNC1 - IC1 - (ZX %*% Bx1))
    }
  } else {
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

# DEZE FUNCTIE MOET NOG AANGEPAST WORDEN
MART <- function(y, x, p_C, p_NC, c) {
  p_CT <- 2*p_C
  p_NCT <- 2*p_NC
  nargin <- length(as.list(match.call())) - 2
  if (is.null(x)){
    x <- "not"
  }
  
  if (length(x) == 1) {
    numcol <- 0
    numcolT <- 0
  } else {
    numcol <- ncol(as.matrix(x))
    numcolT <- numcol*2
  }
  if(numcol > 1){
    x.rev <- matrix(data=NA,nrow=length(x[,1]),ncol=numcol)
    for (i in 1:numcol){
      x.rev[,i] <- rev(x[,i])
    }
  }else{
    x.rev <- matrix(data=NA,nrow=length(x),ncol=numcol)
    x.rev <- rev(x)
  }
  
  
  if (nargin < 5){
    y    <- fBasics::vec(y)
    z    <- rev(y)
    # Hier specificeer je startwaardes voor de parameters voor optimalisatie
    z    <- fBasics::vec(z) # Z hier is basically de toekomst
    BC0  <- arx.ls_T(y,x,p_C,c)[[2]] # Fit een AR model en pak de phi's
    Bx0  <- arx.ls_T(y,x,p_C,c)[[3]] # Fit een AR model en pak de beta's
    BNC0 <- arx.ls_T(z,x.rev,p_NC,c)[[2]] # Fir een AR model op de omgedraaide volgorde, dus basically de toekomst
    IC0  <- 0
    df0  <- 20
    sig0 <- 2
    
    BC0 <- fBasics::vec(BC0)
    BNC0 <- fBasics::vec(BNC0)
    Bx0 <- fBasics::vec(Bx0)
    
    if (length(x) > 1){
      if (p_C > 0 & p_NC > 0){
        params0 <- rbind(BC0,BNC0,Bx0,IC0,sig0,df0)
      } else if (p_NC > 0 & p_C == 0){
        params0 <- rbind(BNC0,Bx0,IC0,sig0,df0)
      } else if (p_C > 0 & p_NC == 0){
        params0 <- rbind(BC0,Bx0,IC0,sig0,df0)
      } else if (p_C == 0 & p_NC == 0){
        params0 <- rbind(Bx0,IC0,sig0,df0)
      }
    } else {
      if (p_C > 0 & p_NC > 0){
        params0 <- rbind(BC0,BNC0,IC0,sig0,df0)
      } else if (p_NC > 0 & p_C == 0){
        params0 <- rbind(BNC0,IC0,sig0,df0)
      } else if (p_C > 0 & p_NC == 0){
        params0 <- rbind(BC0,IC0,sig0,df0)
      } else if (p_C == 0 & p_NC == 0){
        params0 <- rbind(IC0,sig0,df0)
      }
    }
  }
  
  optimization_results <- stats::optim(params0,ll.MART.Z,gr=NULL,y=fBasics::vec(y),p_C=p_C,p_NC=p_NC,x=x,c=c,method="BFGS",hessian=TRUE)
  PARAMS <- optimization_results$par
  
  if (length(x) > 1){
    numcol <- ncol(as.matrix(x))
    ZX <- regressor.matrix_T(y, x, 1, c)[, 3:(2+2*ncol(as.matrix(x)))]
    if (p_C > 0 && p_NC > 0){
      B_C  <- PARAMS[1:(p_CT)]
      B_NC <- PARAMS[(p_CT+1):(p_CT + p_NCT)]
      B_x  <- PARAMS[(p_CT + p_NCT + 1):(p_CT + p_NCT + numcolT)]
      IC   <- PARAMS[(p_CT + p_NCT + numcolT + 1)]
      sig  <- PARAMS[(p_CT + p_NCT + numcolT + 2)]
      df   <- PARAMS[(p_CT + p_NCT + numcol + 3)]
    }
    else if (p_NC > 0 && p_C == 0){
      B_C  <- 0
      B_NC <- PARAMS[1:(p_NCT)]
      B_x  <- PARAMS[(p_NCT + 1):(p_NCT + numcolT)]
      IC   <- PARAMS[(p_NCT + numcolT + 1)]
      sig  <- PARAMS[(p_NCT + numcolT + 2)]
      df   <- PARAMS[(p_NCT + numcolT + 3)]
    }
    else if (p_C > 0 && p_NC == 0){
      B_NC <- 0
      B_C  <- PARAMS[1:(p_CT)]
      B_x  <- PARAMS[(p_CT + 1):(p_CT + numcolT)]
      IC   <- PARAMS[(p_CT + numcolT + 1)]
      sig  <- PARAMS[(p_CT + numcolT + 2)]
      df   <- PARAMS[(p_CT + numcolT + 3)]
    }
    else if(p_C == 0 && p_NC == 0){
      B_NC  <- 0
      B_C   <- 0
      B_x   <- PARAMS[(p_CT + 3):(p_CT + 2 + numcolT)]
      IC    <- PARAMS[(p_CT + numcolT + 3)]
      sig   <- PARAMS[(p_CT + numcolT + 4)]
      df    <- PARAMS[(p_CT + numcolT + 5)]
    }
  } else{
    numcol <- 0
    B_x <- 0
    if (p_C > 0 && p_NC > 0){
      B_C  <- PARAMS[1:(p_CT)]
      B_NC <- PARAMS[(p_CT+1):(p_CT + p_NCT)]
      IC   <- PARAMS[(p_CT + p_NCT + 1)]
      sig  <- PARAMS[(p_CT + p_NCT + 2)]
      df   <- PARAMS[(p_CT + p_NCT + 3)]
    } else if (p_NC > 0 && p_C == 0){
      B_C  <- 0
      B_NC <- PARAMS[1:(p_NCT)]
      IC   <- PARAMS[(2*(p_NC) + 1)]
      sig  <- PARAMS[(2*(p_NC) + 2)]
      df   <- PARAMS[(2*(p_NC) + 3)]
    } else if (p_C > 0 && p_NC == 0){
      B_NC <- 0
      B_C  <- PARAMS[1:(p_CT)]
      IC   <- PARAMS[(p_CT + 1)]
      sig  <- PARAMS[(p_CT + 2)]
      df   <- PARAMS[(p_CT + 3)]
    } else if (p_C == 0 && p_NC == 0){
      B_NC  <- 0
      B_C   <- 0
      IC    <- PARAMS[1]
      sig   <- PARAMS[2]
      df    <- PARAMS[3]
    }
  }
  
  ZC1 <- y[(p_C+1):length(y)]
  ZC1 <- fBasics::vec(ZC1)
  ZC2 <- regressor.matrix_T(y,"not",p_C,c)
  
  if (p_C > 0){
    V <- ZC1 - ZC2 %*% B_C
  } else{
    V <- ZC1
  }
  
  U <- rev(V)
  U <- fBasics::vec(U)
  
  ZNC1 <- U[(p_NC + 1):length(U)]
  ZNC1 <- fBasics::vec(ZNC1)
  ZNC2 <- regressor.matrix_T(U,"not",p_NC,c)
  
  if(numcolT > 1){
    for (i in 1:numcolT){
      ZX[,i] <- rev(ZX[,i])
    }
  }
  
  if(length(x) > 1){
    if (numcolT > 1 ){
      ZX <- ZX[(p_NC +1):length(U),]
    }
    else{
      ZX <- ZX[(p_NC +1):length(U)]
    }
  } else{
    x <- "not"
  }
  
  
  if (length(x) > 1){
    if (p_NC > 0){
      E <- rev(ZNC1 - (ZNC2 %*% B_NC) - IC - (ZX %*% B_x))
    }
    else{
      E <- rev(ZNC1 - IC - (ZX %*% B_x))
    }
  } else{
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
  
  return(list(coef.c1 = B_C[1:p_C], coef.c2 = B_C[(p_C+1):(p_CT)], coef.nc1 = B_NC[1:p_NC], coef.nc2 = B_NC[(p_NC+1):(p_NCT)], coef.exo1 = B_x[1:(length(B_x)/2)], coef.exo2 = B_x[(length(B_x)/2 +1): length(B_x)], coef.int = IC, scale = sig,df = df,residuals = E, se.dist = se.dist))
}

logistic.smooth <- function(y, gamma, c) {
  return(1/(1 + exp(-gamma*(y - c))))
}

#' @title The regressor matrix function for smooth threshold models
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
regressor.matrix_ST <- function(y, x, p, c, gamma) {
  # Handle NULL x case
  if (missing(x) || is.null(x)) {
    x <- "not"
  }
  
  y <- fBasics::vec(y)
  n <- length(y)
  
  if (p == 1) {
    k <- 1
  } else {
    k <- NCOL(y)
  }
  
  # Create Z matrix
  if (p > 0) {
    Z <- matlab::zeros(n, k * p)
    for (i in 1:p) {
      Z[(1 + i):n, ((i - 1) * k + 1):(i * k)] <- y[1:(n - i)]
    }
    Z <- Z[(1 + p):n, ]
  } else {
    Z <- matrix(, nrow = n, ncol = 0)
  }
  
  if (identical(x, "not")) {
  } else if (NCOL(x) == 1) {
    Z <- cbind(Z, x[(1 + p):n])
  } else if (NCOL(x) > 1) {
    Z <- cbind(Z, x[(1 + p):n, ])
  }
  
  # Create thresholded ZT matrix
  if(identical(x, "not")) {
    if(p == 1) {
      nT <- length(Z)
      Z_c <- Z
    } else if(p == 0) {
      nT <- length(Z)
      Z_c <- 0
    } else {
      nT <- nrow(Z)
      Z_c <- Z[,1:p]
    }
  } else {
    nT <- nrow(Z)
    Z_c <- Z[,1:p]
  }
  if(p != 0) {
    Z_c <- cbind(Z_c, Z_c)
    if (!identical(x, "not")) {
      Z_x <- Z[,(p + 1):ncol(Z)]
      Z_x <- cbind(Z_x, Z_x)
      mX <- ncol(Z_x)
      for(i in 1:nT) {
        Z_x_old <- Z_x[i,]
        Z_x[i,1:(mX/2)] <- logistic.smooth(y[i], gamma, c)*Z_x_old[1:(mX/2)]
        Z_x[i, (mX/2 + 1):mX] <- (1-logistic.smooth(y[i], gamma, c))*Z_x_old[(mX/2 + 1):mX]
      }
    }
    mC <- ncol(Z_c)
    for(i in 1:nT) {
      Z_c_old <- Z_c[i,]
      Z_c[i,1:(mC/2)] <- logistic.smooth(y[i], gamma, c)*Z_c_old[1:(mC/2)]
      Z_c[i, (mC/2 + 1):mC] <- (1-logistic.smooth(y[i], gamma, c))*Z_c_old[(mC/2 + 1):mC]
    }
    if (!identical(x, "not")) {
      ZT <- cbind(Z_c, Z_x)
    } else {
      ZT <- Z_c
    }
  } else {
    ZT <- Z
  }
  return(matrix = ZT)
}

#' @title The Smooth transition ARX estimation by OLS function
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

arx.ls_ST <- function(y,x,p,c, gamma){
  
  if (is.null(x)){
    x <- "not"
  }
  
  n <- length(y) - p
  
  Y <- y[(p+1):length(y)]
  int <- rep(1,(length(y)-p))
  ZT <- regressor.matrix_ST(y,x,p,c,gamma)
  ZT <- cbind(int,ZT)
  
  df <- nrow(ZT) - NCOL(ZT)
  
  B <- solve(t(ZT) %*% ZT) %*% (t(ZT) %*% Y)
  
  if (p > 0){
    if (length(x) > 1){
      rownames(B) <- c('int', paste('lag', 1:(2*p)), paste('exo', 1:(2*NCOL(x))))
    }
    else{
      rownames(B) <- c('int', paste('lag', 1:(2*p)))
    }
  }
  else{
    if (length(x) > 1){
      rownames(B) <- c('int', paste('exo', 1:(2*NCOL(x))))
    }
    else{
      rownames(B) <- 'int'
    }
  }
  
  FV <- ZT %*% B
  U <- Y - FV
  
  sig <- (t(U) %*% U)
  sig <- as.numeric(sig)
  
  Cov <- (1/n)*sig
  Cov <- as.numeric(Cov)
  
  sigma2 <- sum((Y - ZT %*% B)^2)/df
  qz <- qr(ZT)
  vcov <- sigma2*chol2inv(qz$qr)
  colnames(vcov) <- rownames(vcov) <- colnames(ZT)
  
  Loglik <- -(n/2)*(1 + log(2*pi)+log(Cov))
  
  if (p == 0){
    B_auto <- 0
  }
  else{
    B_auto <- B[2:(2*p+1)]
  }
  
  if (length(x) > 1){
    B_x <- B[(2*p+2):length(B)]
  }
  else{
    B_x <- 0
  }
  
  return(list(coefficients = B, coef.auto = B_auto, coef.exo = B_x, mse = Cov, residuals = U, loglikelihood = Loglik, fitted.values = FV, df = df,vcov=vcov))
}

ll.SMART.Z <- function(params,y,x,p_C,p_NC,c,gamma) {
  p_CT <- p_C*2
  p_NCT <- p_NC*2
  if (is.null(x)){
    x <- "not"
  }
  
  y <- fBasics::vec(y)
  if (length(x) > 1){
    colnum <- NCOL(as.matrix(x))
    colnumT <- colnum*2
    
    if (p_C > 0 && p_NC > 0){
      BC1  <- params[1:p_CT]
      BNC1 <- params[(p_CT+1):(p_CT + p_NCT)]
      Bx1  <- params[((p_CT + p_NCT)+ 1):(p_CT + p_NCT + colnumT)]
      IC1  <- params[(p_CT + p_NCT + colnumT + 1)]
      sig1 <- params[(p_CT + p_NCT + colnumT + 2)]
      df1  <- params[(p_CT + p_NCT + colnumT + 3)]
    } else if (p_NC > 0 && p_C == 0){
      BC1  <- 0
      BNC1 <- params[1:(p_NCT)]
      Bx1  <- params[((p_NCT)+1):((p_NCT) + colnumT)]
      IC1  <- params[(p_NCT + colnumT + 1)]
      sig1 <- params[(p_NCT + colnumT + 2)]
      df1  <- params[(p_NCT + colnumT + 3)]
    } else if (p_C > 0 && p_NC == 0){
      BNC1 <- 0
      BC1  <- params[1:(p_CT)]
      Bx1  <- params[(p_CT + 1):(p_CT + colnumT)]
      IC1  <- params[(p_CT + colnumT + 1)]
      sig1 <- params[(p_CT + colnumT + 2)]
      df1  <- params[(p_CT + colnumT + 3)]
    } else if (p_C == 0 && p_NC == 0){
      BNC1  <- 0
      BC1   <- 0
      Bx1   <- params[(1:(colnumT))]
      IC1   <- params[(colnumT + 1)]
      sig1  <- params[(colnumT + 2)]
      df1   <- params[(colnumT + 3)]
    }
  } else{
    colnum <- 0
    colnumT <- 0
    if (p_C > 0 && p_NC > 0){
      BC1  <- params[1:(p_CT)]
      BNC1 <- params[((p_CT)+1):(p_CT + p_NCT)]
      IC1  <- params[(p_CT + p_NCT + 1)]
      sig1 <- params[(p_CT + p_NCT + 2)]
      df1  <- params[(p_CT + p_NCT + 3)]
    } else if (p_NC > 0 && p_C == 0){
      BC1  <- 0
      BNC1 <- params[1:(p_NCT)]
      IC1  <- params[(p_NCT + 1)]
      sig1 <- params[(p_NCT + 2)]
      df1  <- params[(p_NCT + 3)]
    } else if (p_C > 0 && p_NC == 0){
      BNC1 <- 0
      BC1  <- params[1:(p_CT)]
      IC1  <- params[(p_CT + 1)]
      sig1 <- params[(p_CT + 2)]
      df1  <- params[(p_CT + 3)]
    } else if (p_C == 0 && p_NC == 0){
      BNC1  <- 0
      BC1   <- 0
      IC1   <- params[1]
      sig1  <- params[2]
      df1   <- params[3]
    }
  }
  
  ZC1 <- y[(p_C+1):length(y)]
  ZC1 <- fBasics::vec(ZC1)
  ZC2 <- regressor.matrix_ST(y,"not",p_C, c, gamma)
  
  if (p_C > 0){
    V <- ZC1 - (ZC2 %*% BC1)
  } else{
    V <- ZC1
  }
  U <- rev(V)
  U <- fBasics::vec(U)
  
  ZNC1 <- U[(p_NC + 1):length(U)]
  ZNC1 <- fBasics::vec(ZNC1)
  ZNC2 <- regressor.matrix_ST(U,"not",p_NC, c, gamma)
  if((colnumT) > 1){
    ZX <- regressor.matrix_ST(y, x, 1, c, gamma)[, 3:(2+2*ncol(as.matrix(x)))]
    for (i in 1:(colnumT)) {
      ZX[,i] <- rev(ZX[,i])
    }
  }
  if (length(x) > 1){
    if ((colnumT) > 1){
      ZX <- ZX[(p_NC +1):length(U),]
    } else {
      ZX <- ZX[(p_NC + 1):length(U)]
    }
  } else {
    x = "not"
  }
  if (length(x) > 1){
    if (p_NC > 0){
      E <- rev(ZNC1 - (ZNC2 %*% BNC1) - IC1 - (ZX %*% Bx1))
    } else{
      E <- rev(ZNC1 - IC1 - (ZX %*% Bx1))
    }
  } else {
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


SMART <- function(y, x, p_C, p_NC, c, gamma) {
  p_CT <- 2*p_C
  p_NCT <- 2*p_NC
  nargin <- length(as.list(match.call())) - 2
  if (is.null(x)){
    x <- "not"
  }
  
  if (length(x) == 1) {
    numcol <- 0
    numcolT <- 0
  } else {
    numcol <- ncol(as.matrix(x))
    numcolT <- numcol*2
  }
  if(numcol > 1){
    x.rev <- matrix(data=NA,nrow=length(x[,1]),ncol=numcol)
    for (i in 1:numcol){
      x.rev[,i] <- rev(x[,i])
    }
  }else{
    x.rev <- matrix(data=NA,nrow=length(x),ncol=numcol)
    x.rev <- rev(x)
  }
  
  
  if (nargin < 6){
    y    <- fBasics::vec(y)
    z    <- rev(y)
    # Hier specificeer je startwaardes voor de parameters voor optimalisatie
    z    <- fBasics::vec(z) # Z hier is basically de toekomst
    BC0  <- arx.ls_ST(y,x,p_C,c, gamma)[[2]] # Fit een AR model en pak de phi's
    Bx0  <- arx.ls_ST(y,x,p_C,c, gamma)[[3]] # Fit een AR model en pak de beta's
    BNC0 <- arx.ls_ST(z,x.rev,p_NC,c, gamma)[[2]] # Fir een AR model op de omgedraaide volgorde, dus basically de toekomst
    IC0  <- 0
    df0  <- 20
    sig0 <- 2
    
    BC0 <- fBasics::vec(BC0)
    BNC0 <- fBasics::vec(BNC0)
    Bx0 <- fBasics::vec(Bx0)
    
    if (length(x) > 1){
      if (p_C > 0 & p_NC > 0){
        params0 <- rbind(BC0,BNC0,Bx0,IC0,sig0,df0)
      } else if (p_NC > 0 & p_C == 0){
        params0 <- rbind(BNC0,Bx0,IC0,sig0,df0)
      } else if (p_C > 0 & p_NC == 0){
        params0 <- rbind(BC0,Bx0,IC0,sig0,df0)
      } else if (p_C == 0 & p_NC == 0){
        params0 <- rbind(Bx0,IC0,sig0,df0)
      }
    } else {
      if (p_C > 0 & p_NC > 0){
        params0 <- rbind(BC0,BNC0,IC0,sig0,df0)
      } else if (p_NC > 0 & p_C == 0){
        params0 <- rbind(BNC0,IC0,sig0,df0)
      } else if (p_C > 0 & p_NC == 0){
        params0 <- rbind(BC0,IC0,sig0,df0)
      } else if (p_C == 0 & p_NC == 0){
        params0 <- rbind(IC0,sig0,df0)
      }
    }
  }
  
  optimization_results <- stats::optim(params0,ll.SMART.Z,gr=NULL,y=fBasics::vec(y),p_C=p_C,p_NC=p_NC,x=x,c=c,gamma=gamma,method="BFGS",hessian=TRUE)
  PARAMS <- optimization_results$par
  
  if (length(x) > 1){
    numcol <- ncol(as.matrix(x))
    ZX <- regressor.matrix_ST(y, x, 1, c, gamma)[, 3:(2+2*ncol(as.matrix(x)))]
    if (p_C > 0 && p_NC > 0){
      B_C  <- PARAMS[1:(p_CT)]
      B_NC <- PARAMS[(p_CT+1):(p_CT + p_NCT)]
      B_x  <- PARAMS[(p_CT + p_NCT + 1):(p_CT + p_NCT + numcolT)]
      IC   <- PARAMS[(p_CT + p_NCT + numcolT + 1)]
      sig  <- PARAMS[(p_CT + p_NCT + numcolT + 2)]
      df   <- PARAMS[(p_CT + p_NCT + numcol + 3)]
    }
    else if (p_NC > 0 && p_C == 0){
      B_C  <- 0
      B_NC <- PARAMS[1:(p_NCT)]
      B_x  <- PARAMS[(p_NCT + 1):(p_NCT + numcolT)]
      IC   <- PARAMS[(p_NCT + numcolT + 1)]
      sig  <- PARAMS[(p_NCT + numcolT + 2)]
      df   <- PARAMS[(p_NCT + numcolT + 3)]
    }
    else if (p_C > 0 && p_NC == 0){
      B_NC <- 0
      B_C  <- PARAMS[1:(p_CT)]
      B_x  <- PARAMS[(p_CT + 1):(p_CT + numcolT)]
      IC   <- PARAMS[(p_CT + numcolT + 1)]
      sig  <- PARAMS[(p_CT + numcolT + 2)]
      df   <- PARAMS[(p_CT + numcolT + 3)]
    }
    else if(p_C == 0 && p_NC == 0){
      B_NC  <- 0
      B_C   <- 0
      B_x   <- PARAMS[(p_CT + 3):(p_CT + 2 + numcolT)]
      IC    <- PARAMS[(p_CT + numcolT + 3)]
      sig   <- PARAMS[(p_CT + numcolT + 4)]
      df    <- PARAMS[(p_CT + numcolT + 5)]
    }
  } else{
    numcol <- 0
    B_x <- 0
    if (p_C > 0 && p_NC > 0){
      B_C  <- PARAMS[1:(p_CT)]
      B_NC <- PARAMS[(p_CT+1):(p_CT + p_NCT)]
      IC   <- PARAMS[(p_CT + p_NCT + 1)]
      sig  <- PARAMS[(p_CT + p_NCT + 2)]
      df   <- PARAMS[(p_CT + p_NCT + 3)]
    } else if (p_NC > 0 && p_C == 0){
      B_C  <- 0
      B_NC <- PARAMS[1:(p_NCT)]
      IC   <- PARAMS[(2*(p_NC) + 1)]
      sig  <- PARAMS[(2*(p_NC) + 2)]
      df   <- PARAMS[(2*(p_NC) + 3)]
    } else if (p_C > 0 && p_NC == 0){
      B_NC <- 0
      B_C  <- PARAMS[1:(p_CT)]
      IC   <- PARAMS[(p_CT + 1)]
      sig  <- PARAMS[(p_CT + 2)]
      df   <- PARAMS[(p_CT + 3)]
    } else if (p_C == 0 && p_NC == 0){
      B_NC  <- 0
      B_C   <- 0
      IC    <- PARAMS[1]
      sig   <- PARAMS[2]
      df    <- PARAMS[3]
    }
  }
  
  ZC1 <- y[(p_C+1):length(y)]
  ZC1 <- fBasics::vec(ZC1)
  ZC2 <- regressor.matrix_ST(y,"not",p_C,c, gamma)
  
  if (p_C > 0){
    V <- ZC1 - ZC2 %*% B_C
  } else{
    V <- ZC1
  }
  
  U <- rev(V)
  U <- fBasics::vec(U)
  
  ZNC1 <- U[(p_NC + 1):length(U)]
  ZNC1 <- fBasics::vec(ZNC1)
  ZNC2 <- regressor.matrix_ST(U,"not",p_NC,c, gamma)
  
  if(numcolT > 1){
    for (i in 1:numcolT){
      ZX[,i] <- rev(ZX[,i])
    }
  }
  
  if(length(x) > 1){
    if (numcolT > 1 ){
      ZX <- ZX[(p_NC +1):length(U),]
    }
    else{
      ZX <- ZX[(p_NC +1):length(U)]
    }
  } else{
    x <- "not"
  }
  
  
  if (length(x) > 1){
    if (p_NC > 0){
      E <- rev(ZNC1 - (ZNC2 %*% B_NC) - IC - (ZX %*% B_x))
    }
    else{
      E <- rev(ZNC1 - IC - (ZX %*% B_x))
    }
  } else{
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
  
  return(list(coef.c1 = B_C[1:p_C], coef.c2 = B_C[(p_C+1):(p_CT)], coef.nc1 = B_NC[1:p_NC], coef.nc2 = B_NC[(p_NC+1):(p_NCT)], coef.exo1 = B_x[1:(length(B_x)/2)], coef.exo2 = B_x[(length(B_x)/2 +1): length(B_x)], coef.int = IC, scale = sig,df = df,residuals = E, se.dist = se.dist))
}