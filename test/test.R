# Load libraries
library(ars)
source("../src/ffbs.R")

# Read file ----
file <- read.csv("test.csv")

# Segment of file
data <- file[1:500,]

shifter <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}

# Observed data
Y <- data$fridge + shifter(data$fridge, 9)*0.65

# Real values of parameters
theta.real <- matrix(c(1000,2000,500,200),1,4)
Z.real <- matrix(as.matrix(data[,c(-1,-6)]), ncol = 4)
for (col in 1:ncol(Z.real)) {
  Z.real[which(Z.real[,col] != 0),col] <- 1
}

plot(Y, type = "l")

# Initialize NFHMM parameters ----
# Hyperparameters
# TODO: Conjugate priors for hyperparameters (hyperpriors)
gamma <- 1 # b, beta distribution
alpha <- 1 # mu, beta distribution, strength parameter of IBP
delta <- 1 # b, beta distribution
mu_theta <- mean(Y, na.rm = T)*2
sigma_epsilon <- sd(Y, na.rm = T)/2
sigma_theta <- sd(Y, na.rm = T)/4

# Parameters
Kdag  <- 2 # Number of active appliances + 1
Z <- matrix(1, length(Y), (Kdag - 1)) # State matrix
Z[sample(1:length(Y), round(0.3*length(Y))),1] <- 0 # Add some zeros
Z <- cbind(Z, matrix(0, nrow(Z), 1)) # Add empty column

o.mu <- matrix(rbeta(Kdag, alpha/Kdag, 1), 1, Kdag) # State transition probability
o.mu <- o.mu[order(o.mu, decreasing = T)] # Order state transition probabilities
b <- matrix(rbeta(Kdag, gamma, delta), 1, Kdag) # State transition probability
theta <- matrix(rnorm(Kdag, mu_theta, sigma_theta), 1, Kdag) # State levels

# Auxiliary functions ----
# Indicator function (inclusive)
indicator <- function(x, upb, lob = 0) {
  if ((x >= lob) && (x <= upb)) {
    return(1)
  } else {
    return(0)
  }
}

# Padding log-distribution for mu_k
fmuk <- function(x, alpha, t, N=10) {
  sum <- 0
  for (i in 1:N) {
    sum <- sum + 1/i*((1-x)^i)
  }
  return((alpha*sum + t*log(1-x) + (alpha-1)*log(x)))    
}

# Derivative of padding log-distribution for mu_k
dfmuk <- function(x, alpha, t, N=10) {
  return(alpha*((1-x)^N-1)/x - t/(1-x) + (alpha-1)/x)
}

# Computation of c's
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

# Iterative sampling for NFHMM ----
# Number of iterations
IterNum <- 4
while (IterNum > 0) {
  
  # K-active: The number of active appliances
  Kact <- max(which(colSums(Z) > 0))
  
  # K-dagger: The index of the highest non-active appliance
  Kdag <- max(which(colSums(Z) == 0))
  
  # Put mu in order
  o.mu <- o.mu[order(o.mu, decreasing = T)]
  
  # Sample auxiliary variable s
  s <- runif(1, min = 0, max = o.mu[Kact])
  
  # K-star: The highest index for which o.mu[K-star] > s
  Kstar <- max(which(o.mu > s))
  
  # Expand representations
  while (Kstar >= Kdag) {
    # Sample mu_k with ARS
    # Starting points for ARS
    x <- runif(3, 0, o.mu[length(o.mu)])
    x <- x[order(x)]
    
    # Sample mu_k
    mu_k <- ars(1,fmuk,dfmuk,x=x,m=3,lb=T,xlb=0,ub=T,
                xub=o.mu[length(o.mu)],alpha=alpha,t=nrow(Z))
    
    # Pad mu_k to ordered mu
    o.mu <- c(o.mu, mu_k)
    o.mu <- o.mu[order(o.mu, decreasing = T)]
    
    # Add empty column to Z
    Z <- cbind(Z, matrix(0, nrow(Z), 1))
    
    # Update K-star
    Kstar <- max(which(o.mu > s))
    
    # Update K-dagger
    Kdag <- max(which(colSums(Z) == 0))
    
    # Sample b_k and theta_k from their priors
    b <- c(b, rbeta(1, gamma, delta))
    theta <- c(theta, rnorm(1, mu_theta, sigma_theta))
    
  }
  
  # Sample Z with blocked Gibbs and run FFBS on each column of Z
  # Note: Z is updated only up to column k <= K*
  # Note: After sampling Z, for any columns k <= K* that are non-active, delete those columns
  for (k in 1:Kstar) {
    W <- matrix(c(1-o.mu[k],o.mu[k],1-b[k],b[k]),byrow = T,nrow=2)
    
    Z <- BSi(Y,Z,k,W,theta,sigma_epsilon,Kact,o.mu)
    Kact <- max(which(colSums(Z) > 0))
  }
  
  # Delete non-active columns for k <= K*
  nonact <- which(colSums(matrix(Z[,1:Kstar], nrow = length(Y))) == 0)
  if (length(nonact) > 0) {
    Z <- Z[,-nonact]
    o.mu <- o.mu[-nonact]
    b <- b[-nonact]
    theta <- theta[-nonact]
  }
  
  # Update K-dagger
  Kdag <- max(which(colSums(Z) == 0))
  
  
  # Sample theta, mu, b, from their conditionals
  # Sample theta
  for (k in 1:Kdag) {
    sigma_theta_p2 <- (1/sigma_theta^2 + 1/sigma_epsilon^2*sum(Z[,k]))^-1
    sum1 <- Z[,k]%*%Y
    thetadk <- theta[-k]
    if (length(thetadk)==1){
      sum2 <- Z[,-k] * thetadk
    } else {
      sum2 <- Z[,-k] %*% thetadk
    }
    
    sum2 <- sum(sum2*Z[,k])
    mu_theta_g <- sigma_theta_p2*((mu_theta/sigma_theta^2) + (sum1-sum2)/sigma_epsilon^2)
    sigma_theta_g <- sqrt(sigma_theta_p2)
    theta[k] <- rnorm(1, mu_theta_g, sigma_theta_g)
  }
  
  # Sample mu_k
  for (k in 1:(Kdag-1)) {
    # Find lower and upper bounds 
    if (k == 1) {
      ubmu <- 1
    } else {
      ubmu <- o.mu[k-1]
    }
    if (is.na(o.mu[k+1])) {
      lbmu <- 0
    } else {
      lbmu <- o.mu[k+1]
    }
    
    # Compute c's
    ck00 <- cfun(0,0,k,Z)
    ck01 <- cfun(0,1,k,Z)
    
    # Draw mu_k using inverse grid sampling
    bound <- seq(lbmu, ubmu, by = 0.0001)
    dbound <- dbeta(bound,ck01,ck00+1)
    dbound <- dbound/sum(dbound)
    cbound <- cumsum(dbound)
    u <- runif(1)
    if (u <= cbound[1]) {
      ind <- 1
    } else if (u >= cbound[length(cbound)]) {
      ind <- length(cbound)
    } else {
      ind <- min(which(cbound > u))
    }
    o.mu[k] <- bound[ind]
  }
  
  # Sample mu_k for k = K-dagger (ARS)
  # Starting points for ARS
  x <- runif(3, 0, o.mu[Kdag-1])
  x <- x[order(x)]
  
  # Sample mu_k
  o.mu[Kdag] <- ars(1,fmuk,dfmuk,x=x,m=3,lb=T,xlb=0,ub=T,
                    xub=o.mu[Kdag-1],alpha=alpha,t=nrow(Z))
  
  
  # Sample b_k
  for (k in 1:Kdag) {
    # Compute c's
    ck11 <- cfun(1,1,k,Z)
    ck10 <- cfun(1,0,k,Z)
    
    # Draw b_k
    b[k] <- rbeta(1, ck11+gamma,ck10+delta+1)
  }
  
  
  # TODO: Sample alpha, gamma, delta, mu_theta, sigma_theta, sigma_epsilon from their posteriors (conjugacy)
  
  IterNum <- IterNum - 1
  
}
# TODO: Why does theta need transpose in some cases?
try(ts.plot(Z%*%t(theta)))
try(ts.plot(Z%*%theta))
ncol(Z)

