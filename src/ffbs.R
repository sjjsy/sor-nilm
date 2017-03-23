C_Yt_Xt <- function(Y,Z,t,k,theta,sigma){
  y <- Y[t]
  z <- Z[t,]
  z0 <- z
  z0[k] <- 0
  z1 <- z
  z1[k] <- 1
  N0 <- dnorm(y,z0%*%t(theta),sigma)
  N1 <- dnorm(y,z1%*%t(theta),sigma)
  N0 <- N0/(N0+N1)
  N1 <- 1-N0
  return(c(N0,N1))
}

FFi <- function(Y,Z,k,W,theta,sigma){
  N <- C_Yt_Xt(Y,Z,1,k,theta,sigma=sigma)
  FF <- matrix(NA,length(Y),2)
  FF[1,1] <- N[1]*W[1,1]
  FF[1,2] <- N[2]*W[1,2]
  FF[1,1] <- FF[1,1]/(sum(FF[1,]))
  FF[1,2] <- 1 - FF[1,1]
  
  
  for (t in 2:length(Y)){
    N <- C_Yt_Xt(Y,Z,t,k,theta,sigma)
    FF[t,1] <- N[1]*(W[1,1]*FF[t-1,1]+W[2,1]*FF[t-1,2])
    FF[t,2] <- N[2]*(W[1,2]*FF[t-1,1]+W[2,2]*FF[t-1,2])
    FF[t,1] <- FF[t,1]/(sum(FF[t,]))
    FF[t,2] <- 1 - FF[t,1]
  }
  
  return(FF)
}

BSi <- function(Y,Z,k,W,theta,sigma,Kact,o.mu){
  #Forward filtering
  FF <- FFi(Y,Z,k,W,theta,sigma=sigma)
  # Calculate multiplier mu^*
  mustar <- o.mu[Kact]
  
  #Check if k  corresponds to the last active column
  if (k == Kact){
    mus <- o.mu[which(colSums(matrix(Z[,1:k], nrow = length(Y)))>0)]
    mustar0 <- min(mus[-which(mus==min(mus))])
    mu.coef <- c(mustar0,mustar)
    f.ind <- min(which(Z[,k]==1))
    
    sampled1 <- F
    trans.coef <- c(1,1)
    for (i in length(Y):1){
      
      if (i < length(Y)){
        trans.coef <- W[,Z[i+1,k]+1]
      }
      if (i > f.ind || sampled1){
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
      if (Zi == 1){
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
      
      if (sampled1){
        props <- FF[i,]*trans.coef
        p1 <- props[2]/sum(props)
        Zi <- rbinom(1,1,p1)
        Z[i,k] <- Zi
      } else {
        mu.coef <- c(mustar,o.mu[k])
        props <- 1/mu.coef*FF[i,]*trans.coef
        p1 <- props[2]/sum(props)
        Zi <- rbinom(1,1,p1)
        Z[i,k] <- Zi
      }
      if (Zi == 1){
        sampled1 <- T
      }
    }
  }
    
  
  return(Z)
}
