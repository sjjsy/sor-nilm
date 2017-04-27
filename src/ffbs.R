C_Yt_Xt <- function(Y,Z,t,k,theta,sigma) {
  y <- Y[t]
  z <- Z[t,]
  z0 <- z
  z0[k] <- 0
  z1 <- z
  z1[k] <- 1
  N0 <- dnorm(y,sum(z0*theta),sigma)
  N1 <- dnorm(y,sum(z1*theta),sigma)
  N0 <- N0/(N0+N1)
  N1 <- 1-N0
  
  v <- c(N0,N1)
  #If both possibilities seem unlikely, do not change
  if (is.nan(N0) & is.nan(N1)) {
    v <- c(1,0)
    if (z[k] == z1[k]) {
      v <- c(0,1)
    }
  }
  return(v)
}

FFi <- function(Y,Z,k,W,theta,sigma) {
  N <- C_Yt_Xt(Y,Z,1,k,theta,sigma)
  FF <- matrix(NA,length(Y),2)
  FF[1,1] <- N[1]*W[1,1]
  FF[1,2] <- N[2]*W[1,2]
  FF[1,1] <- FF[1,1]/(sum(FF[1,]))
  FF[1,2] <- 1 - FF[1,1]
  
  
  for (t in 2:length(Y)) {
    N <- C_Yt_Xt(Y,Z,t,k,theta,sigma)
    FF[t,1] <- N[1]*(W[1,1]*FF[t-1,1]+W[2,1]*FF[t-1,2])
    FF[t,2] <- N[2]*(W[1,2]*FF[t-1,1]+W[2,2]*FF[t-1,2])
    FF[t,1] <- FF[t,1]/(sum(FF[t,]))
    FF[t,2] <- 1 - FF[t,1]
  }
  return(FF)
}

BSi <- function(Y,Z,k,W,theta,sigma,Kact,mu) {
  #Forward filtering
  FF <- FFi(Y,Z,k,W,theta,sigma)
  # Calculate multiplier mu^*
  mustar <- mu[Kact]
  
  #Check if k corresponds to the last active column
  if (k == Kact) {
    mus <- mu[which(colSums(matrix(Z[,1:k], nrow = length(Y)))>0)]
    if (length(mus) == 1) {
      mustar0 <- 1
    } else {
      mustar0 <- min(min(mus[-which(mus==min(mus))]),mustar)
    }
    mu.coef <- c(mustar0,mustar)
    f.ind <- min(which(Z[,k]==1))
    
    sampled1 <- F
    trans.coef <- c(1,1)
    for (i in length(Y):1) {
      
      if (i < length(Y)) {
        trans.coef <- W[,Z[i+1,k]+1]
      }
      if (i > f.ind || sampled1) {
        props <- FF[i,]*trans.coef
        p1 <- props[2]/sum(props)
        Zi <- rbinom(1,1,p1)
        Z[i,k] <- Zi
      } else {
        props <- 1/mu.coef*FF[i,]*trans.coef
        p1 <- props[2]/sum(props)
        Zi <- rbinom(1,1,p1)
        Z[i,k] <- Zi
      }
      if (Zi == 1) {
        sampled1 <- T
      }
      
    }
  }
  else if (k < Kact){
    trans.coef <- c(1,1)
    for (i in length(Y):1){
      if (i < length(Y)){
        trans.coef <- W[,Z[i+1,k]+1]
      }
      props <- FF[i,]*trans.coef
      p1 <- props[2]/sum(props)
      Z[i,k] <- rbinom(1,1,p1)
    }
  } else {
    sampled1 <- F
    trans.coef <- c(1,1)
    for (i in length(Y):1){
      if (i < length(Y)){
        trans.coef <- W[,Z[i+1,k]+1]
      }
      
      if (sampled1) {
        props <- FF[i,]*trans.coef
        p1 <- props[2]/sum(props)
        Zi <- rbinom(1,1,p1)
        Z[i,k] <- Zi
      } else {
        mu.coef <- c(mustar,mu[k])
        props <- 1/mu.coef*FF[i,]*trans.coef
        p1 <- props[2]/sum(props)
        Zi <- rbinom(1,1,p1)
        Z[i,k] <- Zi
      }
      if (Zi == 1) {
        sampled1 <- T
      }
    }
  }
  return(Z)
}
