# Auxiliary functions ----
# Padding log-distribution for mu_k
fmuk <- function(mu, alpha, t, N = 10) {
  sum <- 0
  for (i in 1:N) {
    sum <- sum + 1/i*((1-mu)^i)
  }
  return((alpha*sum + t*log(1-mu) + (alpha-1)*log(mu)))
}

# Derivative of padding log-distribution for mu_k
dfmuk <- function(mu, alpha, t, N = 10) {
  return(alpha*((1-mu)^N-1)/mu - t/(1-mu) + (alpha-1)/mu)
}

# log-functions for new mu's (ARS)
bmuk <- function(mu,ck00,ck01) {
  h <- ck00*log(1-mu)+(ck01-1)*log(mu)
  return(h)
}

dbmuk <- function(mu,ck00,ck01) {
  dh <- -ck00/(1-mu) + (ck01-1)/mu
  return(dh)
}

# Computation of "c" constants
cfun <- function(i, j, k, Z) {
  column <- Z[,k]
  i.inds <- which(column == i)
  c.val <- 0
  for (i.ind in i.inds) {
    if (i.ind < length(column)) {
      if (column[i.ind + 1] == j) {
        c.val <- c.val + 1
      }
    }
  }
  if (i == 0) {
    c.val <- c.val + 1 # Because of 0-th dummy column
  }
  return(c.val)
}

grid_sample_alpha <- function(mu,Kdag,lb=1,ub=5) {
  
  grid <- seq(lb,ub,by=0.001)
  post <- rep(1,length(grid))
  
  for (i in 1:length(grid)) {
    for (k in 1:Kdag) {
      post[i] <- post[i]*dbeta(mu[k],grid[i]/Kdag,1)
    }
  }
  
  post <- post/sum(post)
  cpost <- cumsum(post)
  
  u <- runif(1)
  ind <- min(which(cpost > u))
  #ind <- which.max(post) # Maximum likelihood method
  return(grid[ind])
}

grid_sample_gamma <- function(b,Kdag,delta,lb=0,ub=5) {
  
  grid <- seq(lb,ub,by=0.001)
  post <- rep(1,length(grid))
  
  for (i in 1:length(grid)) {
    for (k in 1:Kdag) {
      post[i] <- post[i]*dbeta(b[k],grid[i],delta)
    }
  }
  
  post <- post/sum(post)
  cpost <- cumsum(post)
  
  u <- runif(1)
  ind <- min(which(cpost > u))
  #ind <- which.max(post) # Maximum likelihood method
  return(grid[ind])
}

grid_sample_delta <- function(b,Kdag,gamma,lb=0,ub=5) {
  
  grid <- seq(lb,ub,by=0.001)
  post <- rep(1,length(grid))
  
  for (i in 1:length(grid)) {
    for (k in 1:Kdag) {
      post[i] <- post[i]*dbeta(b[k],gamma,grid[i])
    }
  }
  
  post <- post/sum(post)
  cpost <- cumsum(post)
  
  u <- runif(1)
  ind <- min(which(cpost > u))
  #ind <- which.max(post) # Maximum likelihood method
  return(grid[ind])
}

