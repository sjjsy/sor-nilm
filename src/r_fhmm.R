# Load libraries
library(stats)
library(ars)

# Read file ----
file <- read.csv("../in/mdr_2016-03_export_99.csv", sep = ";")

# Data preprocessing ----
# Check number of different DeviceIDs
if (length(unique(file$DeviceId)) > 1) {
  warning("Multiple DeviceIDs in given file.")
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

# Initialize NFHMM parameters ----
# Hyperparameters
# TODO: Conjugate priors for hyperparameters (hyperpriors)
alpha <- 1.5 # mu, beta distribution, strength parameter of IBP
gamma <- 1.5 # b, beta distribution
delta <- 1 # b, beta distribution
mu_theta <- mean(data$S, na.rm = T)
sigma_theta <- sd(data$S, na.rm = T)*2
sigma_epsilon <- sd(data$S, na.rm = T)

# Parameters
Kdag  <- 2 # Number of active appliances + 1
Z <- matrix(1, nrow(data), (Kdag - 1)) # State matrix
Z[sample(1:nrow(data), round(0.3*nrow(data))),1] <- 0
Z <- cbind(Z, matrix(0, nrow(Z), 1))

mu <- matrix(rbeta(Kdag, alpha/Kdag, 1), 1, Kdag) # State transition probability
o.mu <- mu[order(mu, decreasing = T)]
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

# Iterative sampling for NFHMM ----
# Number of iterations
IterNum <- 100
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
  
  # TODO: Sample theta, mu, b, from their conditionals
  
  # TODO: Sample alpha, gamma, delta, mu_theta, sigma_theta, sigma_epsilon from their posteriors (conjugacy)
  
  
  IterNum <- IterNum - 1
}
