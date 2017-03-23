# Load libraries
library(stats)
library(ars)
source("ffbs.R")

# Read file ----
file <- read.csv("../in/mdr_2016-03_export_99.csv", sep = ";")
# Data preprocessing ----
# Check number of different DeviceIDs
if (length(unique(file$DeviceId)) > 1) {
  warning("Multiple DeviceIDs in input file.")
}

# Store each phase to a separate data frame
f1 <- file[which(file$TypeId == 1),]
f2 <- file[which(file$TypeId == 2),]
f3 <- file[which(file$TypeId == 3),]

# Create data frame for data
data <- data.frame(matrix(NA, length(unique(file$ValueTimeStamp)), 6))
colnames(data) <- c("DeviceId", "ValueTimeStamp", "P1", "P2", "P3", "S")
data[["DeviceId"]] <- file$DeviceId[1:nrow(data)]
data[["ValueTimeStamp"]] <- unique(file$ValueTimeStamp)

# Match each phase to its time
# This is done only once so duplicates are ignored
data[match(f1$ValueTimeStamp, data$ValueTimeStamp),"P1"] <- f1$DataValue
data[match(f2$ValueTimeStamp, data$ValueTimeStamp),"P2"] <- f2$DataValue
data[match(f3$ValueTimeStamp, data$ValueTimeStamp),"P3"] <- f3$DataValue
data[,"S"] <- data[,"P1"] + data[,"P2"] + data[,"P3"]

# Delete empty data rows
data <- data[-which(is.na(data[,"S"])),]

# Initialize NFHMM parameters ----
# Hyperparameters
# TODO: Conjugate priors for hyperparameters (hyperpriors)
alpha <- 1.5 # mu, beta distribution, strength parameter of IBP
gamma <- 1.5 # b, beta distribution
delta <- 1 # b, beta distribution
mu_theta <- mean(data$S, na.rm = T)
sigma_theta <- sd(data$S, na.rm = T)
sigma_epsilon <- sd(data$S, na.rm = T)

# Parameters
Kdag  <- 2 # Number of active appliances + 1
Z <- matrix(1, nrow(data), (Kdag - 1)) # State matrix
Z[sample(1:nrow(data), round(0.3*nrow(data))),1] <- 0 # Add some zeros
Z <- cbind(Z, matrix(0, nrow(Z), 1)) # Add empty column
Y <- data$S # Observed data

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

# log-functions for new mu's (ARS)
bmuk <- function(x,ck00,ck01){
  h <- ck00*log(1-x)+(ck01-1)*log(x)
  return(h)
}
dbmuk <- function(x,ck00,ck01){
  dh <- -ck00/(1-x)+(ck01-1)/x
  return(dh)
}

# Iterative sampling for NFHMM ----
# Number of iterations
IterNum <- 1
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
  
  # TODO: Sample Z with blocked Gibbs and run FFBS on each column of Z
  # Note: Z is updated only up to column k <= K*
  # Note: After sampling Z, for any columns k <= K* that are non-active, delete those columns
  for (k in 1:Kstar) {
    W <- matrix(c(1-o.mu[k],o.mu[k],1-b[k],b[k]),byrow = T,nrow=2)

    Z <- BSi(Y,Z,k,W,theta,sigma=sigma_epsilon,Kact,o.mu)
    Kact <- max(which(colSums(Z) > 0))
    }

  
  # Sample theta, mu, b, from their conditionals
    
  for (k in 1:Kdag) {
    # Sample theta
    sigma_theta_p2 <- (1/sigma_theta^2 + 1/sigma_epsilon^2*sum(Z[,k]))^-1
    sum1 <- Z[,k]*Y
    sum2 <- matrix(Z[,-k]) %*% matrix(theta, nrow = 1)[,-k]
    sum2 <- sum2*Z[,k]
    mu_theta_g <- sigma_theta_p2*((mu_theta/sigma_theta^2) + (sum1-sum2)/sigma_epsilon^2)
    sigma_theta_g <- sqrt(sigma_theta_p2)
    theta[k] <- rnorm(1, mu_theta_g, sigma_theta_g)
  }
    
  for (k in 1:(Kdag)) {
    # Sample mu_k
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
    
  for (k in 1:Kdag) {
    # Sample b_k
    # Compute c's
    ck11 <- cfun(1,1,k,Z)
    ck10 <- cfun(1,0,k,Z)
    
    # Draw b_k
    b[k] <- rbeta(1, ck11+gamma,ck10+1)
  }
    
  
  # TODO: Sample alpha, gamma, delta, mu_theta, sigma_theta, sigma_epsilon from their posteriors (conjugacy)
  
  IterNum <- IterNum - 1

}
