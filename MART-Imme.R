
# DEZE FUNCTIE MOET NOG AANGEPAST WORDEN
MART <- function(y, x, p_C, p_NC) {
  
  # Vgm hoeft dit niet aangepast te worden toch?
  
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
  
  # Dit weet ik niet precies hoe dit moet...
  
  if (nargin < 5){
    y    <- fBasics::vec(y)
    z    <- rev(y)
    z    <- fBasics::vec(z)
    BC0_1  <- arx.ls(y,x,p_C)[[2]] # Fit een AR model en pak de phi's
    BC0_2  <- arx.ls(y,x,p_C)[[2]] # ??
    
    # Deze moet niet in de twee regimes geswitcht worden toch?
    Bx0  <- arx.ls(y,x,p_C)[[3]] # Fit een AR model en pak de beta's
    
    BNC0_1 <- arx.ls(z,x.rev,p_NC)[[2]] # Fit een AR model op de omgedraaide volgorde, dus basically de toekomst
    BNC0_2 <- arx.ls(z,x.rev,p_NC)[[2]] # ??
    IC0  <- 0
    df0  <- 20
    sig0 <- 2
    
    
    # Dit heb ik wel aangepast weer, even checken of t klopt :)
    BC0_1 <- fBasics::vec(BC0_1)
    BC0_2 <- fBasics::vec(BC0_2)
    BNC0_1 <- fBasics::vec(BNC0_1)
    BNC0_2 <- fBasics::vec(BNC0_2)
    Bx0 <- fBasics::vec(Bx0)
    
    if (length(x) > 1){
      if (p_C > 0 && p_NC > 0){
        params0 <- rbind(BC0_1,BC0_2,BNC0_1,BNC0_2,Bx0,IC0,sig0,df0)
      }
      else if (p_NC > 0 && p_C == 0){
        params0 <- rbind(BNC0_1,BNC0_2,Bx0,IC0,sig0,df0)
      }
      else if (p_C > 0 && p_NC == 0){
        params0 <- rbind(BC0_1,BC0_2,Bx0,IC0,sig0,df0)
      }
      else if (p_C == 0 && p_NC == 0){
        params0 <- rbind(Bx,IC0,sig0,df0)
      }
    }
    else{
      if (p_C > 0 && p_NC > 0){
        params0 <- rbind(BC0_1,BC0_2,BNC0_1,BNC0_2,IC0,sig0,df0)
      }
      else if (p_NC > 0 && p_C == 0){
        params0 <- rbind(BNC0_1,BNC0_2,IC0,sig0,df0)
      }
      else if (p_C > 0 && p_NC == 0){
        params0 <- rbind(BC0_1,BC0_2,IC0,sig0,df0)
      }
      else if (p_C == 0 && p_NC == 0){
        params0 <- rbind(IC0,sig0,df0)
      }
    }
  }
  
  # Dit blijft het zelfde denk ik?
  optimization_results <- stats::optim(params0,ll.max,gr=NULL,y=fBasics::vec(y),p_C=p_C,p_NC=p_NC,x=x,method="BFGS",hessian=TRUE)
  PARAMS <- optimization_results$par
  
  # Vanaf hier heb ik t aangepast naar MART, maar wel nog even naar kijken hoor twijfel of t klopt :))
  if (length(x) > 1){
    numcol <- NCOL(x)
    
    if (p_C > 0 && p_NC > 0){
      B_C1  <- PARAMS[1:p_C]
      B_C2  <- PARAMS[(p_C+1):2*p_C]
      B_NC1 <- PARAMS[(2*p_C+1):(2*p_C + p_NC)]
      B_NC2 <- PARAMS[(2*p_C + p_NC + 1):(2*(p_C + p_NC))]
      B_x   <- PARAMS[(2*(p_C+ p_NC) + 1):(p_C + p_NC + numcol)]
      IC    <- PARAMS[(2*(p_C + p_NC) + numcol + 1)]
      sig   <- PARAMS[(2*(p_C + p_NC) + numcol + 2)]
      df    <- PARAMS[(2*(p_C + p_NC) + numcol + 3)]
    }
    
    else if (p_NC > 0 && p_C == 0){
      B_C1  <- 0
      B_C2  <- 0
      B_NC1 <- PARAMS[1:p_NC]
      B_NC2 <- PARAMS[(p_NC+1):(2*p_NC)]
      B_x   <- PARAMS[(2*p_NC + 1):(2*p_NC + numcol)]
      IC    <- PARAMS[(2*p_NC + numcol + 1)]
      sig   <- PARAMS[(2*p_NC + numcol + 2)]
      df    <- PARAMS[(2*p_NC + numcol + 3)]
    }
    else if (p_C > 0 && p_NC == 0){
      B_NC1 <- 0
      B_NC2 <- 0
      B_C1  <- PARAMS[1:p_C]
      B_C2  <- PARAMS[(p_C + 1):(2*p_C)]
      B_x  <- PARAMS[(2*p_C + 1):(2*p_C + numcol)]
      IC  <- PARAMS[(2*p_C + numcol + 1)]
      sig <- PARAMS[(2*p_C + numcol + 2)]
      df  <- PARAMS[(2*p_C + numcol + 3)]
    }
    else if (p_C == 0 && p_NC == 0){
      BNC1 <- 0
      BNC2 <- 0
      B_C1  <- 0
      B_C2  <- 0
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
      B_C1  <- PARAMS[1:p_C]
      B_C2  <- PARAMS[(p_C+1):2*p_C]
      B_NC1 <- PARAMS[(2*p_C+1):(2*p_C + p_NC)]
      B_NC2 <- PARAMS[(2*p_C + p_NC + 1):(2*(p_C + p_NC))]
      IC    <- PARAMS[(2*(p_C + p_NC) + 1)]
      sig   <- PARAMS[(2*(p_C + p_NC) + 2)]
      df    <- PARAMS[(2*(p_C + p_NC) + 3)]
    }
    else if (p_NC > 0 && p_C == 0){
      B_C1  <- 0
      B_C2  <- 0
      B_NC1 <- PARAMS[1:p_NC]
      B_NC2 <- PARAMS[(p_NC+1):(2*p_NC)]
      IC    <- PARAMS[(2*p_NC + 1)]
      sig   <- PARAMS[(2*p_NC + 2)]
      df    <- PARAMS[(2*p_NC + 3)]
    }
    else if (p_C > 0 && p_NC == 0){
      BNC1 <- 0
      BNC2 <- 0
      BC1  <- PARAMS[1:p_C]
      BC2  <- PARAMS[(p_C + 1):(2*p_C)]
      IC  <- PARAMS[(2*p_C + 1)]
      sig <- PARAMS[(2*p_C + 2)]
      df  <- PARAMS[(2*p_C + 3)]
    }
    else if (p_C == 0 && p_NC == 0){
      BNC1 <- 0
      BNC2 <- 0
      B_C1  <- 0
      B_C2  <- 0
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
    V1 <- ZC1 - ZC2 %*% B_C1
    V2 <- ZC1 - ZC2 %*% B_C2
    
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
  
  
  if(numcol > 1){
    for (i in 1:numcol){
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
  
  se <- sqrt(diag(solve(optimization_results$hessian)))
  se.dist <- se[(length(se)-1):length(se)]
  se.dist <- rev(se.dist)
  
  return(list(coef.c1 = B_C1, coef.c2 = B_C2, coef.nc1 = B_NC1, coef.nc2 = B_NC2, coef.exo = B_x, coef.int = IC, scale = sig,df = df,residuals = E, se.dist = se.dist))
}