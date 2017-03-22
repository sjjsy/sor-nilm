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
  mustar <- min(o.mu[which(colSums(Z)>0)])
  
  #Check if k  corresponds to the last active column
  if (k == Kact){
    # Check if z_{1:(end-1),k} are zeros
    if (sum(Z[1:(length(Y)-1),k])==0){
      #New multiplier for p(z_{end,k}==0)
      mus <- o.mu[which(colSums(Z[,1:k])>0)]
      mustar0 <- min(mus[-which(mus==min(mus))])
      p1 <- 1/mustar*FF[length(Y),2]/(1/mustar*FF[length(Y),2]+1/mustar0*FF[length(Y),1])
      Z[length(Y),k] <- rbinom(1,1,p1) }
     else{
      Z[length(Y),k] <- rbinom(1,1,FF[length(Y),2])
    }
    #For next z_{t,k} we need to check if both z_{1:(t-1),k} and z_{(t+1):end,k} are zeros
      for (i in (length(Y)-1):1){
        if (sum(Z[1:(i-1),k])==0 & sum(Z[(i+1):length(Y),k])==0){
          #New multiplier for p(z_{t,k}==0)
          mus <- o.mu[which(colSums(Z[,1:k])>0)]
          mustar0 <- min(mus[-which(mus==min(mus))])
          props <- c(1/mustar0,1/mustar)*FF[i,]*(W[,Z[i+1,k]+1])
          p1 <- props[2]/sum(props)
          Z[i,k] <- rbinom(1,1,p1) }
       else {
           props <- FF[i,]*(W[,Z[i+1,k]+1])
           p1 <- props[2]/sum(props)
           Z[i,k] <- rbinom(1,1,p1) } 
      }
    
  }
  # Not the last active column
  else {

  Z[length(Y),k] <- rbinom(1,1,FF[length(Y),2])
  for (i in (length(Y)-1):1){
    props <- FF[i,]*(W[,Z[i+1,k]+1])
    p1 <- props[2]/sum(props)
    Z[i,k] <- rbinom(1,1,p1)
  }
  }
  
  return(Z)
}
