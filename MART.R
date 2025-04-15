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
regressor.matrix_T <- function(y, x, p, c, d=1) {
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
    Z <- Z[(1 + max(p,d)):n, ]
  } else {
    Z <- matrix(, nrow = n-d, ncol = 0)
  }
  
  if (identical(x, "not")) {
  } else if (NCOL(x) == 1) {
    Z <- cbind(Z, x[(1 + max(p,d)):n])
  } else if (NCOL(x) > 1) {
    Z <- cbind(Z, x[(1 + max(p,d)):n, ])
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
        if(y[max(p,d)+i-d] > c) {
          Z_x[i, 1:(mX/2)] <- 0
        } else {
          Z_x[i, (mX/2 + 1):mX] <- 0
        }
      }
    }
    mC <- ncol(Z_c)
    for(i in 1:nT) {
      if(y[max(p,d)+i-d] > c) {
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

arx.ls_T <- function(y,x,p,c,d=1){
  
  if (is.null(x)){
    x <- "not"
  }
  
  n <- length(y) - p
  
  Y <- y[(max(p,d)+1):length(y)]
  int <- rep(1,(length(y)-max(p,d)))
  ZT <- regressor.matrix_T(y,x,p,c,d)
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
ll.MART.Z <- function(params,y,x,p_C,p_NC,c,d=1){
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
  
  ZC1 <- y[(max(p_C,d)+1):length(y)]
  ZC1 <- fBasics::vec(ZC1)
  ZC2 <- regressor.matrix_T(y,"not",p_C, c, d)
  
  if (p_C > 0){
    V <- ZC1 - (ZC2 %*% BC1)
  } else{
    V <- ZC1
  }
  U <- rev(V)
  U <- fBasics::vec(U)
  
  ZNC1 <- U[(max(p_NC,d) + 1):length(U)]
  ZNC1 <- fBasics::vec(ZNC1)
  ZNC2 <- regressor.matrix_T(U,"not",p_NC, c,d)
  if((colnumT) > 1){
    ZX <- regressor.matrix_T(y, x, 1, c,d)[, 3:(2+2*ncol(as.matrix(x)))]
    for (i in 1:(colnumT)) {
      ZX[,i] <- rev(ZX[,i])
    }
  }
  if (length(x) > 1){
    if ((colnumT) > 1){
      ZX <- ZX[(max(p_NC,d) +1):length(U),]
    } else {
      ZX <- ZX[(max(p_NC,d) + 1):length(U)]
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
MART <- function(y, x, p_C, p_NC, c, d=1) {
  p_CT <- 2*p_C
  p_NCT <- 2*p_NC
  nargin <- length(as.list(match.call())) - 3
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
  
  
  if (nargin < 7){
    y    <- fBasics::vec(y)
    z    <- rev(y)
    z    <- fBasics::vec(z)
    BC0  <- arx.ls_T(y,x,p_C,c,d)[[2]] 
    Bx0  <- arx.ls_T(y,x,p_C,c,d)[[3]]
    BNC0 <- arx.ls_T(z,x.rev,p_NC,c,d)[[2]]
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
    ZX <- regressor.matrix_T(y, x, 1, c,d)[, 3:(2+2*ncol(as.matrix(x)))]
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
  
  ZC1 <- y[(max(p_C,d)+1):length(y)]
  ZC1 <- fBasics::vec(ZC1)
  ZC2 <- regressor.matrix_T(y,"not",p_C,c,d)
  
  if (p_C > 0){
    V <- ZC1 - ZC2 %*% B_C
  } else{
    V <- ZC1
  }
  
  U <- rev(V)
  U <- fBasics::vec(U)
  
  ZNC1 <- U[(max(p_NC,d) + 1):length(U)]
  ZNC1 <- fBasics::vec(ZNC1)
  ZNC2 <- regressor.matrix_T(U,"not",p_NC,c,d)
  
  if(numcolT > 1){
    for (i in 1:numcolT){
      ZX[,i] <- rev(ZX[,i])
    }
  }
  
  if(length(x) > 1){
    if (numcolT > 1 ){
      ZX <- ZX[(max(p_NC,d) +1):length(U),]
    }
    else{
      ZX <- ZX[(max(p_NC,d) +1):length(U)]
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
  se.dist <- se
  var.cov <- solve(optimization_results$hessian)
  
  
  return(list(coef.c1 = B_C[1:p_C], coef.c2 = B_C[(p_C+1):(p_CT)], coef.nc1 = B_NC[1:p_NC], coef.nc2 = B_NC[(p_NC+1):(p_NCT)], coef.exo1 = B_x[1:(length(B_x)/2)], coef.exo2 = B_x[(length(B_x)/2 +1): length(B_x)], coef.int = IC, scale = sig,df = df,residuals = E, se.dist = se.dist, var.cov = var.cov))
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
regressor.matrix_ST <- function(y, x, p, c, gamma,d=1) {
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
    Z <- Z[(1 + max(p,d)):n, ]
  } else {
    Z <- matrix(, nrow = n-d, ncol = 0)
  }
  
  if (identical(x, "not")) {
  } else if (NCOL(x) == 1) {
    Z <- cbind(Z, x[(1 + max(p,d)):n])
  } else if (NCOL(x) > 1) {
    Z <- cbind(Z, x[(1 + max(p,d)):n, ])
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
      Z_x <- Z[,(max(p,d) + 1):ncol(Z)]
      Z_x <- cbind(Z_x, Z_x)
      mX <- ncol(Z_x)
      for(i in 1:nT) {
        Z_x_old <- Z_x[i,]
        Z_x[i,1:(mX/2)] <- logistic.smooth(y[max(p,d)+i-d], gamma, c)*Z_x_old[1:(mX/2)]
        Z_x[i, (mX/2 + 1):mX] <- (1-logistic.smooth(y[max(p,d)+i-d], gamma, c))*Z_x_old[(mX/2 + 1):mX]
      }
    }
    mC <- ncol(Z_c)
    for(i in 1:nT) {
      Z_c_old <- Z_c[i,]
      Z_c[i,1:(mC/2)] <- logistic.smooth(y[max(p,d)+i-d], gamma, c)*Z_c_old[1:(mC/2)]
      Z_c[i, (mC/2 + 1):mC] <- (1-logistic.smooth(y[max(p,d)+i-d], gamma, c))*Z_c_old[(mC/2 + 1):mC]
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

arx.ls_ST <- function(y,x,p,c, gamma,d=1){
  
  if (is.null(x)){
    x <- "not"
  }
  
  n <- length(y) - p
  
  Y <- y[(max(p,d)+1):length(y)]
  int <- rep(1,(length(y)-max(p,d)))
  ZT <- regressor.matrix_ST(y,x,p,c,gamma,d)
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

ll.SMART.Z <- function(params,y,x,p_C,p_NC,c,gamma,d=1) {
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
  ZC2 <- regressor.matrix_ST(y,"not",p_C, c, gamma,d)
  
  if (p_C > 0){
    V <- ZC1 - (ZC2 %*% BC1)
  } else{
    V <- ZC1
  }
  U <- rev(V)
  U <- fBasics::vec(U)
  
  ZNC1 <- U[(max(p_NC,d) + 1):length(U)]
  ZNC1 <- fBasics::vec(ZNC1)
  ZNC2 <- regressor.matrix_ST(U,"not",p_NC, c, gamma,d)
  if((colnumT) > 1){
    ZX <- regressor.matrix_ST(y, x, 1, c, gamma,d)[, 3:(2+2*ncol(as.matrix(x)))]
    for (i in 1:(colnumT)) {
      ZX[,i] <- rev(ZX[,i])
    }
  }
  if (length(x) > 1){
    if ((colnumT) > 1){
      ZX <- ZX[(max(p_NC,d) +1):length(U),]
    } else {
      ZX <- ZX[(max(p_NC,d) + 1):length(U)]
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


SMART <- function(y, x, p_C, p_NC, c, gamma,d=1) {
  p_CT <- 2*p_C
  p_NCT <- 2*p_NC
  nargin <- length(as.list(match.call())) - 3
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
  
  
  if (nargin < 7){
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
    ZX <- regressor.matrix_ST(y, x, 1, c, gamma,d)[, 3:(2+2*ncol(as.matrix(x)))]
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
  
  ZC1 <- y[(max(p_C,d)+1):length(y)]
  ZC1 <- fBasics::vec(ZC1)
  ZC2 <- regressor.matrix_ST(y,"not",p_C,c, gamma,d)
  
  if (p_C > 0){
    V <- ZC1 - ZC2 %*% B_C
  } else{
    V <- ZC1
  }
  
  U <- rev(V)
  U <- fBasics::vec(U)
  
  ZNC1 <- U[(max(p_NC,d) + 1):length(U)]
  ZNC1 <- fBasics::vec(ZNC1)
  ZNC2 <- regressor.matrix_ST(U,"not",p_NC,c, gamma,d)
  
  if(numcolT > 1){
    for (i in 1:numcolT){
      ZX[,i] <- rev(ZX[,i])
    }
  }
  
  if(length(x) > 1){
    if (numcolT > 1 ){
      ZX <- ZX[(max(p_NC,d) +1):length(U),]
    }
    else{
      ZX <- ZX[(max(p_NC,d) +1):length(U)]
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
  se.dist <- se
  
  return(list(coef.c1 = B_C[1:p_C], coef.c2 = B_C[(p_C+1):(p_CT)], coef.nc1 = B_NC[1:p_NC], coef.nc2 = B_NC[(p_NC+1):(p_NCT)], coef.exo1 = B_x[1:(length(B_x)/2)], coef.exo2 = B_x[(length(B_x)/2 +1): length(B_x)], coef.int = IC, scale = sig,df = df,residuals = E, se.dist = se.dist))
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

computePath <- function(y, p, h, d, c, c1, c2, nc1, nc2, intercept, v1, v2) {
  path <- rep(NA, h + d)
  p <- length(which(y > c))/length(c)
  n <- length(y)
  r <- length(c1)
  s <- length(c2)
  y_path <- rep(NA, h + d)
  for(i in range(1:(h+d))) {
    path[i] <- binomial(rbinom(1, 1, (1-p))) + 1
    if(y[n + i - d] > c) {
      y_part <- intercept + c1 * y[(n + i -1):(n + i - r)] 
    } else {
      y_part <- intercept + c2 * y[(n + i -1):(n + i - r)] 
    }
    if(path[i] == 1) {
      y_path[i] <- y_part + v1
    } else {
      y_path[i] <- y_part + v2 
      }
    # Check if the series is consistent so far
    if (sum(!is.na(x)) > d) {
      checks <- sum(!is.na(x))  - d
    }
  }
}

#' @title Forecasting function for the MART model
#' @description   This function allows you to forecast with the mixed causal-noncausal model with possibly exogenous regressors.
#' @param y       Data vector y.
#' @param X       (optional) Matrix with data (column represent a series).
#' @param p_C     Number of lags (causal order).
#' @param p_NC    Number of leads (noncausal order).
#' @param X.for   (optional) Matrix with forecasted values for X (column represents series).
#' @param h       Forecast horizon h.
#' @param M       (optional) Truncation value M for MA representation. Default value: 50.
#' @param N       (optional) Number of simulations to forecast noncausal component. Default: 10,000.
#' @param Seed    (optional) Seed for random number generation. Default: 9999
#' @keywords forecasting
#' @return \item{y.for}{Vector containing forecasted values for y.}
#' @author Sean Telg, Floris van den Doel
#' @export
#' @examples


# Treshold verandert niet fundamenteel de simulatie van toekomstige errors, je moet alleen uitkijken dat dimensies enzo kloppen
forecast.MART <- function(y,X,p_C,p_NC,c,d,X.for,h,M,N,seed=20240402) {
  
  set.seed(seed)
  if (missing(X) == TRUE){
    X = NULL
  }
  
  if (missing(N) == TRUE){
    N = 1000
  }
  
  model <- MART(y, X, p_C, p_NC, c, d)
  obs <- length(y)
  
  ## Check whether there are exogenous variables and whether truncation M is known
  
  if (missing(X.for) == TRUE && missing(M) == TRUE){
    X.for = NULL
    M = 50
  } else if(missing(X.for) == TRUE && missing(M) == FALSE){
    X.for = NULL
    M = M
  } else if(missing(X.for) == FALSE && missing(M) == TRUE){
    if (NCOL(X.for) == 1){
      if(is.null(X.for) == TRUE){
        M = 50
      } else{
        M = length(X.for)
      }
    } else{
      M = length(X.for[,1])
    }
  } else if(missing(X.for) == FALSE && missing(M) == FALSE){
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
  
  r <- length(model$coef.c1)
  s <- length(model$coef.nc1)
  ## Simulate future epsilon and use forecasted X
  hve_reg1 <- c()
  hve_reg2 <- c()
  hve21 <- matrix(data=0, nrow=N,ncol=h)
  hve22 <- matrix(data=0, nrow=N,ncol=h)
  
  for (iter in 1:N){
    
    eps.sim <- model$scale*stats::rt(M,model$df)
    z2 <- c()
    for (i in 1:M){
      if(is.null(X.for) == TRUE){
        z2[i] <- eps.sim[i]
      }
      else{
        # Dit op de een of andere manier omzetten naar threshold maar zie niet 123 hoe
        if(NCOL(X.for) > 1){
          z2[i] <- eps.sim[i] +  coef.exo %*% t(X.for[i,])
        }
        else{
          z2[i] <- eps.sim[i] + coef.exo * X.for[i]
        }
      }
    }
    
    ## Compute filtered values u = phi(L)y and moving average values
    # split in phi1 en phi2 voor twee regimes.
    phi1 <- c(1,model$coef.c1)
    phi2 <- c(1,model$coef.c2)
    
    u <- c()
    for (i in (max(r,d)+1):obs){
      if(y[i-d] >c) {
        u[i] <- phi1 %*% y[i:(i-r)]
      } else {
        u[i] <- phi2 %*% y[i:(i-r)]
      }
    }
    w <- c(u[(obs-s+1):obs],z2)
    # C is a function of the non-causal polynomial. So we need to split it up for now
    C1 <- matrix(data=0, nrow=(M+s), ncol=(M+s))
    C2 <- matrix(data=0, nrow=(M+s), ncol=(M+s))
    C1[1,] <- compute.MA(model$coef.nc1,(M+s-1))
    C2[1,] <- compute.MA(model$coef.nc2,(M+s-1))
    if (C1[1,M+s] > 1 | C2[1,M+s] > 1) {
      fallback_forecast <- forecast.MART(y, X, (p_C + p_NC), 0, c, d, X.for, h, M, N, seed)
      return(list(forecast = fallback_forecast$forecast, defaulted = TRUE))
    }
    if (s > 1){
      for (i in 2:s){
        C1[i,] <- c(0, C1[(i-1),1:(length(C1[(i-1),])-1)])
        C2[i,] <- c(0, C2[(i-1),1:(length(C2[(i-1),])-1)])
      }
    }
    
    for (i in (s+1):(M+s)){
      C1[i,] <- c(rep(0,(i-1)),1,rep(0,(M+s-i)))
      C2[i,] <- c(rep(0,(i-1)),1,rep(0,(M+s-i)))
    }
    
    D1 = solve(C1)
    D2 = solve(C2)
    
    e1 <- D1 %*% w
    e2 <- D2 %*% w
    
    h1 <- c()
    h2 <- c()
    
    for (i in 1:s){
      h1[i] <- metRology::dt.scaled(e1[i], df=model$df, sd=model$scale)
      h2[i] <- metRology::dt.scaled(e2[i], df=model$df, sd=model$scale)
    }
    
    hve_reg1[iter] = prod(h1)
    hve_reg2[iter] = prod(h2)
    
    # At this point all necesarry info is there, only regime path needs to be established
    # Take each regime as a bernoulli distribution
    r <- rbinom((d+h), 1, (1-p)) + 1
    for (j in 1:h){
      mov.av1 <-  C1[1,1:(M-j+1)] %*% z2[j:M]
      mov.av2 <- C2[1,1:(M-j+1)] %*% z2[j:M]
      
      hve21[iter,j] <- mov.av1 * hve_reg1[iter]
      hve22[iter,j] <- mov.av2 * hve_reg2[iter]
      
    }
    
    y.star <- y[(obs-r+1):obs]
    y.for <- c()
    exp1 <- c()
    exp2 <- c()
    
    for (j in 1:h){
      exp1[j] = ((1/N)*sum(hve21[,j]))/((1/N)*sum(hve_reg1))
      exp2[j] = ((1/N)*sum(hve22[,j]))/((1/N)*sum(hve_reg2))
      
      if(y[(obs - d + j)] > c) {
        p <- length(which(y > c))/length(y)
        if(length(model$coef.c1) == 1){
          y.for[j] <- model$coef.c1 * y.star + p *((model$coef.int/(1-sum(model$coef.nc1))  + exp1[j])) + (1-p)*((model$coef.int/(1-sum(model$coef.nc2))  + exp2[j]))
        } else{
          y.for[j] <-  t(model$coef.c1) %*% y.star + p *((model$coef.int/(1-sum(model$coef.nc1))  + exp1[j])) + (1-p)*((model$coef.int/(1-sum(model$coef.nc2))  + exp2[j]))
        }
      } else {
        p <- length(which(y > c))/length(y)
        if(length(model$coef.c1) == 1){
          y.for[j] <- model$coef.c2 * y.star + p *((model$coef.int/(1-sum(model$coef.nc1))  + exp1[j])) + (1-p)*((model$coef.int/(1-sum(model$coef.nc2))  + exp2[j]))
        } else{
          y.for[j] <-  t(model$coef.c2) %*% y.star + p *((model$coef.int/(1-sum(model$coef.nc1))  + exp1[j])) + (1-p)*((model$coef.int/(1-sum(model$coef.nc2))  + exp2[j]))
        }
      }
      y.star <- c(y.for[j], y.star[1:(length(y.star)-1)])
    }
  }
  return(list(forecast = y.for, defaulted = FALSE))
}

information.criteria <- function(type = c("MARX", "MART", "SMART"), model) {
  n <- length(model$residuals)
  
  k <- if (type == "MARX") {
    sum(lengths(list(model$coef.c, model$coef.nc, model$coef.exo, model$coef.int))) + 2
  } else {
    2 * sum(lengths(list(model$coef.c1, model$coef.nc1, model$coef.exo1))) + 2
  }
  
  df <- model$df
  sig <- model$scale
  loglikelihood <- (n*lgamma((df+1)/2) - n*log(sqrt(df*pi*sig^2)) - n*lgamma(df/2) - ((df+1)/2)*sum(log(1+(model$residuals/sig)^2/df)))
  
  aic <- (-2 * loglikelihood + 2 * k) / n
  bic <- (-2 * loglikelihood + log(n) * k) / n
  hq  <- (-2 * loglikelihood + 2 * log(log(n)) * k) / n
  
  return(list(aic = aic, bic = bic, hq = hq, loglikelihood = loglikelihood, k = k, n = n, df = df, sig = sig))
}

# Treshold verandert niet fundamenteel de simulatie van toekomstige errors, je moet alleen uitkijken dat dimensies enzo kloppen
forecast.SMART <- function(y,X,p_C,p_NC,c,gamma,d,X.for,h,M,N,seed=20240402) {
  
  set.seed(seed)
  if (missing(X) == TRUE){
    X = NULL
  }
  
  if (missing(N) == TRUE){
    N = 10000
  }
  
  model <- SMART(y, X, p_C, p_NC, c, gamma, d)
  obs <- length(y)
  
  ## Check whether there are exogenous variables and whether truncation M is known
  
  if (missing(X.for) == TRUE && missing(M) == TRUE){
    X.for = NULL
    M = 50
  } else if(missing(X.for) == TRUE && missing(M) == FALSE){
    X.for = NULL
    M = M
  } else if(missing(X.for) == FALSE && missing(M) == TRUE){
    if (NCOL(X.for) == 1){
      if(is.null(X.for) == TRUE){
        M = 50
      } else{
        M = length(X.for)
      }
    } else{
      M = length(X.for[,1])
    }
  } else if(missing(X.for) == FALSE && missing(M) == FALSE){
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
  
  r <- length(model$coef.c1)
  s <- length(model$coef.nc1)
  ## Simulate future epsilon and use forecasted X
  hve_reg1 <- c()
  hve_reg2 <- c()
  hve21 <- matrix(data=0, nrow=N,ncol=h)
  hve22 <- matrix(data=0, nrow=N,ncol=h)
  
  for (iter in 1:N){
    
    eps.sim <- model$scale*stats::rt(M,model$df)
    z2 <- c()
    for (i in 1:M){
      if(is.null(X.for) == TRUE){
        z2[i] <- eps.sim[i]
      }
      else{
        # Dit op de een of andere manier omzetten naar threshold maar zie niet 123 hoe
        if(NCOL(X.for) > 1){
          z2[i] <- eps.sim[i] +  coef.exo %*% t(X.for[i,])
        }
        else{
          z2[i] <- eps.sim[i] + coef.exo * X.for[i]
        }
      }
    }
    
    ## Compute filtered values u = phi(L)y and moving average values
    # split in phi1 en phi2 voor twee regimes.
    phi1 <- c(1,model$coef.c1)
    phi2 <- c(1,model$coef.c2)
    
    u <- c()
    for (i in (max(r,d)+1):obs){
      if(y[i-d] >c) {
        u[i] <- phi1 %*% y[i:(i-r)]
      } else {
        u[i] <- phi2 %*% y[i:(i-r)]
      }
    }
    w <- c(u[(obs-s+1):obs],z2)
    # C is a function of the non-causal polynomial. So we need to split it up for now
    C1 <- matrix(data=0, nrow=(M+s), ncol=(M+s))
    C2 <- matrix(data=0, nrow=(M+s), ncol=(M+s))
    C1[1,] <- compute.MA(model$coef.nc1,(M+s-1))
    C2[1,] <- compute.MA(model$coef.nc2,(M+s-1))
    
    #if (C1[1,M+s] > 1 | C2[1,M+s] > 1) {
    #  fallback_forecast <- forecast.MART(y, X, (p_C + p_NC), 0, c, d, X.for, h, M, N, seed)
    #  return(list(forecast = fallback_forecast$forecast, defaulted_to_ar = TRUE))
    #}
    
  
    if (s > 1){
      for (i in 2:s){
        C1[i,] <- c(0, C1[(i-1),1:(length(C1[(i-1),])-1)])
        C2[i,] <- c(0, C2[(i-1),1:(length(C2[(i-1),])-1)])
      }
    }
    
    for (i in (s+1):(M+s)){
      C1[i,] <- c(rep(0,(i-1)),1,rep(0,(M+s-i)))
      C2[i,] <- c(rep(0,(i-1)),1,rep(0,(M+s-i)))
    }
    
    D1 = solve(C1)
    D2 = solve(C2)
    
    e1 <- D1 %*% w
    e2 <- D2 %*% w
    
    h1 <- c()
    h2 <- c()
    
    for (i in 1:s){
      h1[i] <- metRology::dt.scaled(e1[i], df=model$df, sd=model$scale)
      h2[i] <- metRology::dt.scaled(e2[i], df=model$df, sd=model$scale)
    }
    
    hve_reg1[iter] = prod(h1)
    hve_reg2[iter] = prod(h2)
    
    for (j in 1:h){
      mov.av1 <-  C1[1,1:(M-j+1)] %*% z2[j:M]
      mov.av2 <- C2[1,1:(M-j+1)] %*% z2[j:M]
      
      hve21[iter,j] <- mov.av1 * hve_reg1[iter]
      hve22[iter,j] <- mov.av2 * hve_reg2[iter]
      
    }
  }
  
  y.star <- y[(obs-r+1):obs]
  y.for <- c()
  exp1 <- c()
  exp2 <- c()
  
  for (j in 1:h){
    exp1[j] = ((1/N)*sum(hve21[,j]))/((1/N)*sum(hve_reg1))
    exp2[j] = ((1/N)*sum(hve22[,j]))/((1/N)*sum(hve_reg2))
    
    p <- logistic.smooth(mean(y), gamma, c)
    if(length(model$coef.c1) == 1){
      y.for[j] <- logistic.smooth(y[(obs - d + j)], c, gamma) * model$coef.c1 * y.star + (1-logistic.smooth(y[(obs - d + j)], c, gamma)) * model$coef.c2 * y.star + p *((model$coef.int/(1-sum(model$coef.nc1))  + exp1[j])) + (1-p)*((model$coef.int/(1-sum(model$coef.nc2))  + exp2[j]))
    } else{
      y.for[j] <-  logistic.smooth(y[(obs - d + j)], c, gamma) * t(model$coef.c1) %*% y.star + (1 - logistic.smooth(y[(obs - d + j)], c, gamma)) * t(model$coef.c2) %*% y.star  + p *((model$coef.int/(1-sum(model$coef.nc1))  + exp1[j])) + (1-p)*((model$coef.int/(1-sum(model$coef.nc2))  + exp2[j]))
    }

    y.star <- c(y.for[j], y.star[1:(length(y.star)-1)])
  }
  return(list(forecast = y.for, defaulted_to_ar = FALSE))
  
}
