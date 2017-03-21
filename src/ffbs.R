C_Yt_Xt <- function(Y,Z,t,k,theta,sigma){
  y <- Y[t]
  z <- Z[t,]
  z0 <- z
  z0[k] <- 0
  z1 <- z
  z1[k] <- 1
  N0 <- dnorm(y,z0%*%theta,sigma)
  N1 <- dnorm(y,z1%*%theta,sigma)
  N0 <- N0/(N0+N1)
  N1 <- 1-N0
  return(c(N0,N1))
}

FFi <- function(Y,Z,k,W,theta,sigma){
  N <- C_Yt_Xt(Y,Z,1,k,theta,sigma)
  FF <- matrix(NA,length(Y),2)
  FF[1,1] <- N[1]*W[1,1]
  FF[1,2] <- N[2]*W[1,2]
  FF[1,1] <- FF[1,1]/(sum(FF[1,]))
  FF[1,2] <- 1 - FF[1,1]
  
  
  for (t in 2:length(Y)){
    N <- C_Yt_Xt(Y,Z,t,k,theta,sigma)
    FF[t,1] <- N[1]*(W[1,1]*FF[t-1,1]+W[2,1]*FF[t-1,2])
    FF[t,2] <- N[2]*(W(1,2)*FF[t-1,1]+W[2,2]*FF[t-1,2])
    FF[t,1] <- FF[t,1]/(sum(FF[t,]))
    FF[t,2] <- 1 - FF[t,1]
  }
  
  return(FF)
}

BSi <- function(Y,Z,k,W,theta,sigma){
  
  FF <- FFi(Y,Z,k,theta,sigma)
  Z[length(Y),k] <- rbinom(1,1,FF[length(Y),2])
  for (i in (length(Y)-1):1){
    props <- FF[i,]*(W[,Z[length(Y),k]+1])
    p1 <- props[2]/sum(props)
    Z[i,k] <- rbinom(1,1,p1)
  }
  
  return(Z)
}
