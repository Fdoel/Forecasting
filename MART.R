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

regressor.matrix <- function(y,x,p) {
  
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

ll.max <- function(params,y,x,p_C,p_NC){
  
  if (is.null(x)){
    x <- "not"
  }
  
  y <- fBasics::vec(y)
  if (length(x) > 1){
    colnum <- NCOL(x)
    
    if (p_C > 0 && p_NC > 0){
      BC1  <- params[1:p_C]
      BC2 <- params[(p_C+1):2*p_C]
      BNC1 <- params[(2*p_C+1):(2*p_C + p_NC)]
      BNC2 <- params[(2*p_C + p_NC + 1):(2*(p_C + p_NC))]
      Bx1  <- params[(2*(p_C+ p_NC) + 1):(p_C + p_NC + colnum)]
      IC1  <- params[(2*(p_C + p_NC) + colnum + 1)]
      sig1 <- params[(2*(p_C + p_NC) + colnum + 2)]
      df1  <- params[(2*(p_C + p_NC) + colnum + 3)]
    }
    else if (p_NC > 0 && p_C == 0){
      BC1  <- 0
      BC2  <- 0
      BNC1 <- params[1:p_NC]
      BNC2 <- params[(p_NC+1):(2*p_NC)]
      Bx1  <- params[(2*p_NC + 1):(2*p_NC + colnum)]
      IC1  <- params[(2*p_NC + colnum + 1)]
      sig1 <- params[(2*p_NC + colnum + 2)]
      df1  <- params[(2*p_NC + colnum + 3)]
    }
    else if (p_C > 0 && p_NC == 0){
      BNC1 <- 0
      BNC2 <- 0
      BC1  <- params[1:p_C]
      BC2  <- params[(p_C + 1):(2*p_C)]
      Bx1  <- params[(2*p_C + 1):(2*p_C + colnum)]
      IC1  <- params[(2*p_C + colnum + 1)]
      sig1 <- params[(2*p_C + colnum + 2)]
      df1  <- params[(2*p_C + colnum + 3)]
    }
    else if (p_C == 0 && p_NC == 0){
      BNC1  <- 0
      BNC2  <- 0
      BC1   <- 0
      BC2   <- 0
      Bx1   <- params[1:colnum]
      IC1   <- params[(colnum + 1)]
      sig1  <- params[(colnum + 2)]
      df1   <- params[(colnum + 3)]
    }
  }
  
  else{
    colnum <- 0
    
    if (p_C > 0 && p_NC > 0){
      BC1 <- params[1:p_C]
      BC2 <- params[(p_C+1):2*p_C]
      BNC1 <- params[(2*p_C+1):(2*p_C + p_NC)]
      BNC2 <- params[(2*p_C + p_NC + 1):(2*(p_C + p_NC))]
      IC1 <- params[(2*(p_C + p_NC) + 1)]
      sig1 <- params[(2*(p_C + p_NC) + 2)]
      df1 <- params[(2*(p_C + p_NC) + 3)]
    }
    else if (p_NC > 0 && p_C == 0){
      BC1  <- 0
      BC2  <- 0
      BNC1 <- params[1:p_NC]
      BNC2 <- params[(p_NC+1):(2*p_NC)]
      IC1  <- params[(2*p_NC + 1)]
      sig1 <- params[(2*p_NC + 2)]
      df1  <- params[(2*p_NC + 3)]
    }
    else if (p_C > 0 && p_NC == 0){
      BNC1 <- 0
      BNC2 <- 0
      BC1  <- params[1:p_C]
      BC2  <- params[(p_C + 1):(2*p_C)]
      IC1  <- params[(2*p_C + 1)]
      sig1 <- params[(2*p_C + 2)]
      df1  <- params[(2*p_C + 3)]
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
    V1 <- ZC1 - (ZC2 %*% BC1)
    V2 <- ZC1 - (ZC2 %*% BC2)
  }
  else{
    V1 <- ZC1
    V2 <- ZC1
  }
  
  U1 <- rev(V1)
  U2 <- rev(V2)
  U1 <- fBasics::vec(U1)
  U2 <- fBasics::vec(U2)
  
  ZNC1 <- U[(p_NC + 1):length(U)]
  ZNC1 <- fBasics::vec(ZNC1)
  ZNC21 <- regressor.matrix(U1,"not",p_NC)
  ZNC22 <- regressor.matrix(U2,"not",p_NC)
  
  
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
      x <- x[(p_NC +1):length(U1),]
    }
    else{
      x <- x[(p_NC + 1):length(U1)]
      x <- fBasics::vec(x)
    }
  }
  else{
    x = "not"
  }
  
  if (p_NC == 1){
    ZNC21 <- fBasics::vec(ZNC21)
    ZNC22 <- fBasics::vec(ZNC22)
  }
  
  if (length(x) > 1){
    if (p_NC > 0){
      E1 <- rev(ZNC1 - (ZNC2 %*% BNC1) - IC1 - (x %*% Bx1))
      E2 <- rev(ZNC1 - (ZNC2 %*% BNC2) - IC1 - (x %*% Bx1))
    }
    else{
      E1 <- rev(ZNC1 - IC1 - (x %*% Bx1))
      E2 <- rev(ZNC1 - IC1 - (x %*% Bx1))
    }
  }
  else{
    if (p_NC > 0){
      E1 <- rev(ZNC1 - (ZNC21 %*% BNC1) - IC1)
      E2 <- rev(ZNC1 - (ZNC22 %*% BNC2) - IC1)
    }
    else{
      E1 <- rev(ZNC1 - IC1)
      E2 <- rev(ZNC1 - IC1)
    }
  }
  
  n <- length(E)
  
  loglik_eval <- -(n*lgamma((df1+1)/2) - n*log(sqrt(df1*pi*sig1^2)) - n*lgamma(df1/2) - ((df1+1)/2)*log(1+((E1 * vec(E1 >= 0)  + (1 - vec(E1 >= 0)) * E2)/sig1)^2/df1) %*% matlab::ones(n,1))
  
  return(neg.loglikelihood = loglik_eval)
}
MART <- function(y, x, p_C, p_NC) {
  
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
}