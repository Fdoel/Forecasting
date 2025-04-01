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
#' 
#' Inshallah dit werkt
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
    } else {
    nT <- nrow(Z)
    Z_c <- Z[,1:p]
    }
  }
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
  if (is.null(x)){
    x <- "not"
  }
  
  y <- fBasics::vec(y)
  if (length(x) > 1){
    colnum <- NCOL(x)
    
    if (p_C > 0 && p_NC > 0){
      BC1  <- params[1:(2*p_C)]
      BNC1 <- params[((2*p_C)+1):(2*(p_C + p_NC))]
      Bx1  <- params[((2*(p_C+ p_NC) )+ 1):(2*(p_C + p_NC) + 2*colnum)]
      IC1  <- params[(2*(p_C + p_NC) + 2*colnum + 1)]
      sig1 <- params[(2*(p_C + p_NC) + 2*colnum + 2)]
      df1  <- params[(2*(p_C + p_NC) + 2*colnum + 3)]
    } else if (p_NC > 0 && p_C == 0){
      BC1  <- 0
      BNC1 <- params[1:(2*p_NC)]
      Bx1  <- params[((2*p_NC)+1):((2*p_NC) + 2*colnum)]
      IC1  <- params[(2*p_NC + 2*colnum + 1)]
      sig1 <- params[(2*p_NC + 2*colnum + 2)]
      df1  <- params[(2*p_NC + 2*colnum + 3)]
    } else if (p_C > 0 && p_NC == 0){
      BNC1 <- 0
      BC1  <- params[1:(2*p_C)]
      Bx1  <- params[(2*p_C + 1):(2*p_C + 2*colnum)]
      IC1  <- params[(2*p_C + 2*colnum + 1)]
      sig1 <- params[(2*p_C + 2*colnum + 2)]
      df1  <- params[(2*p_C + 2*colnum + 3)]
    } else if (p_C == 0 && p_NC == 0){
      BNC1  <- 0
      BC1   <- 0
      Bx1   <- params[(1:2*colnum)]
      IC1   <- params[(2*colnum + 1)]
      sig1  <- params[(2*colnum + 2)]
      df1   <- params[(2*colnum + 3)]
    }
  } else{
    colnum <- 0
    if (p_C > 0 && p_NC > 0){
      BC1  <- params[1:(2*p_C)]
      BNC1 <- params[((2*p_C)+1):(2*(p_C + p_NC))]
      IC1  <- params[(2*(p_C + p_NC) + 1)]
      sig1 <- params[(2*(p_C + p_NC) + 2)]
      df1  <- params[(2*(p_C + p_NC) + 3)]
    } else if (p_NC > 0 && p_C == 0){
      BC1  <- 0
      BNC1 <- params[1:(2*p_NC)]
      IC1  <- params[(2*p_NC + 1)]
      sig1 <- params[(2*p_NC + 2)]
      df1  <- params[(2*p_NC + 3)]
    } else if (p_C > 0 && p_NC == 0){
      BNC1 <- 0
      BC1  <- params[1:(2*p_C)]
      IC1  <- params[(2*p_C + 1)]
      sig1 <- params[(2*p_C + 2)]
      df1  <- params[(2*p_C + 3)]
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
  if (p_C == 1){
    ZC2 <- fBasics::vec(ZC2)
  }

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
  
  if(colnum > 1){
    for (i in 1:colnum){
      x[,i] <- rev(x[,i])
    }
  } else {
    x <- rev(x)
  }
  
  if (length(x) > 1){
    if (colnum > 1){
      x <- x[(p_NC +1):length(U),]
    } else {
      x <- x[(p_NC + 1):length(U)]
      x <- fBasics::vec(x)
    }
  } else {
    x = "not"
  }
  
  if (p_NC == 1){
    ZNC2 <- fBasics::vec(ZNC2)
  }
  
  if (length(x) > 1){
    
    if (p_NC > 0){
      E <- rev(ZNC1 - (ZNC2 %*% BNC1) - IC1)
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

# DEZE FUNCTIE MOET NOG AANGEPAST WORDEN
MART <- function(y, x, p_C, p_NC, c) {
  #print(match.call())
  nargin <- length(as.list(match.call())) - 2
  
  if (is.null(x)){
    x <- "not"
  }
  
  if (length(x) == 1) {
    numcol <- 0
  } else {
    numcol <- NCOL(x)
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
    numcol <- ncol(x)
    
    if (p_C > 0 && p_NC > 0){
      B_C  <- PARAMS[1:(2*p_C)]
      B_NC <- PARAMS[(2*p_C+1):(2*p_C + 2*p_NC)]
      B_x  <- PARAMS[(2*(p_C + p_NC) + 1):(2*(p_C + p_NC) + 2*numcol)]
      IC   <- PARAMS[(2*(p_C + p_NC) + 2*numcol + 1)]
      sig  <- PARAMS[(2*(p_C + p_NC) + 2*numcol + 2)]
      df   <- PARAMS[(2*(p_C + p_NC) + numcol + 3)]
    }
    else if (p_NC > 0 && p_C == 0){
      B_C  <- 0
      B_NC <- PARAMS[1:(2*p_NC)]
      B_x  <- PARAMS[(2*p_NC + 1):(2*p_NC + 2*numcol)]
      IC   <- PARAMS[(2*p_NC + 2*numcol + 1)]
      sig  <- PARAMS[(2*p_NC + 2*numcol + 2)]
      df   <- PARAMS[(2*p_NC + 2*numcol + 3)]
    }
    else if (p_C > 0 && p_NC == 0){
      B_NC <- 0
      B_C  <- PARAMS[1:(2*p_C)]
      B_x  <- PARAMS[(2*p_C + 1):(2*p_C + 2*numcol)]
      IC   <- PARAMS[(2*p_C + 2*numcol + 1)]
      sig  <- PARAMS[(2*p_C + 2*numcol + 2)]
      df   <- PARAMS[(2*p_C + 2*numcol + 3)]
    }
    else if(p_C == 0 && p_NC == 0){
      B_NC  <- 0
      B_C   <- 0
      B_x   <- PARAMS[(2*p_C + 3):(2*p_C + 2 + 2*numcol)]
      IC    <- PARAMS[(2*p_C + 2*numcol + 3)]
      sig   <- PARAMS[(2*p_C + 2*numcol + 4)]
      df    <- PARAMS[(2*p_C + 2*numcol + 5)]
    }
  } else{
    numcol <- 0
    B_x <- 0
    if (p_C > 0 && p_NC > 0){
      B_C  <- PARAMS[1:(2*p_C)]
      B_NC <- PARAMS[(2*p_C+1):(2*(p_C + p_NC))]
      IC   <- PARAMS[(2*(p_C + p_NC) + 1)]
      sig  <- PARAMS[(2*(p_C + p_NC) + 2)]
      df   <- PARAMS[(2*(p_C + p_NC) + 3)]
    }
    else if (p_NC > 0 && p_C == 0){
      B_C  <- 0
      B_NC <- PARAMS[1:(2*p_NC)]
      IC   <- PARAMS[(2*(p_NC) + 1)]
      sig  <- PARAMS[(2*(p_NC) + 2)]
      df   <- PARAMS[(2*(p_NC) + 3)]
    }
    else if (p_C > 0 && p_NC == 0){
      B_NC <- 0
      B_C  <- PARAMS[1:2*(p_C)]
      IC   <- PARAMS[(2*p_C + 1)]
      sig  <- PARAMS[(2*p_C + 2)]
      df   <- PARAMS[(2*p_C + 3)]
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
  ZC2 <- regressor.matrix_T(y,"not",p_C,c)
  
  if (p_C == 1){
    ZC2 <- fBasics::vec(ZC2)
  }
  
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
  
  if(numcol > 1){
    for (i in 1:numcol){
      x[,i] <- rev(x[,i])
    }
  }else{
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
  } else{
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
  
  return(list(coef.c1 = B_C[1:p_C], coef.c2 = B_C[(p_C+1):(2*p_C)], coef.nc1 = B_NC[1:p_NC], coef.nc2 = B_NC[(p_NC+1):(2*p_NC)], coef.exo1 = B_x[1:(length(B_x)/2)], coef.exo2 = B_x[(length(B_x)/2 +1): length(B_x)], coef.int = IC, scale = sig,df = df,residuals = E, se.dist = se.dist))
}
