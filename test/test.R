# Load libraries
library(ars)
source("src/ffbs.R")

# Color palette by Paul Tol
tol5 <- c("#1B9E77", "#66A61E", "#7570B3", "#D95F02", "#E6AB02", "#E7298A")



# Read file ----
file <- read.csv("test/test.csv")

# Segment of file

data <- file[6000:8300,]
# Observed data
Y <- data$aggregate

# Add noise
Y <- Y + rnorm(length(Y), 0, 20)

# Real values of parameters
theta.real <- c(1000,2000,500,200)
Z.real <- matrix(as.matrix(data[,c(-1,-6)]), ncol = 4)
for (col in 1:ncol(Z.real)) {
  Z.real[which(Z.real[,col] != 0),col] <- 1
}

# Plot observed signal and its components
plot(Y, type = "l", ylab = "Power (W)", xlab = "Time (s)", col = tol5[1])
lines(Z.real[,1]*theta.real[1], type = "l", lty = 1, col = tol5[2], lwd = 1.2)
lines(Z.real[,2]*theta.real[2], type = "l", lty = 1, col = tol5[3], lwd = 1.2)
lines(Z.real[,3]*theta.real[3], type = "l", lty = 1, col = tol5[4], lwd = 1.2)
lines(Z.real[,4]*theta.real[4], type = "l", lty = 1, col = tol5[5], lwd = 1.2)
legend("topleft", c("Aggregate","Oven","Heater","Fridge","Light"),
       fill = c(tol5[1],tol5[2],tol5[3],tol5[4],tol5[5]))

# Initialize NFHMM parameters ----
# Hyperparameters
# TODO: Conjugate priors for hyperparameters (hyperpriors)
alpha <- 4
gamma <- 1
delta <- 1
mu_theta <- mean(Y, na.rm = T)
sigma_epsilon <- 80#sd(Y, na.rm = T)/8
sigma_theta <- sd(Y, na.rm = T)

# Parameters
Kdag <- 2 # Number of active appliances + 1
Z <- matrix(1, length(Y), (Kdag - 1)) # State matrix
Z[sample(1:length(Y), round(0.3*length(Y))),1] <- 0 # Add some zeros
Z <- cbind(Z, matrix(0, nrow(Z), 1)) # Add empty column

mu <- rbeta(Kdag, alpha/Kdag, 1) # State transition probability
mu <- mu[order(mu, decreasing = T)] # Order state trabsition probabilities
b <- rbeta(Kdag, gamma, delta) # State transition probability
theta <- rnorm(Kdag, mu_theta, sigma_theta) # State levels

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

# Iterative sampling for NFHMM ----
# Number of iterations
IterNum.total <- 10
IterNum <- IterNum.total

# Should each iteration be plotted?
iter.plot <- T
while (IterNum > 0) {
  # K-active: The number of active appliances
  # Active appliances
  ActApp <- which(colSums(Z) > 0)
  if (length(ActApp) == 0) {
    Kact <- 0
  } else {
    Kact <- max(ActApp)
  }
  
  # Sample auxiliary variable s
  # Define upper bounds (mu-star)
  if (Kact == 0) {
    ub <- 1
  } else {
    ub <- mu[Kact]
  }
  # Draw auxiliary variable s
  s <- runif(1, 0, ub)
  
  # K-dagger: The index of the highest non-active appliance
  Kdag <- max(which(colSums(Z) == 0))
  
  # Put parameters in correct order
  order.mu <- order(mu, decreasing = T)
  mu <- mu[order.mu]
  b <- b[order.mu]
  theta <- theta[order.mu]
  
  # K-star: The highest index for which mu[K-star] > s
  Kstar <- max(which(mu > s))
  
  # Expand representations
  while (Kstar >= Kdag) {
    # Sample mu_k with ARS
    # Starting points and bounds for ARS
    ub <- mu[length(mu)]
    sp <- runif(4, 0, ub)
    sp <- sp[order(sp)]
    # Draw mu_k with ARS
    mu_k <- ars(1, fmuk, dfmuk, x = sp, m = 4, lb = T, xlb = 0,
                ub = T, xub = ub, alpha = alpha, t = nrow(Z))
    
    # Sample b_k and theta_k from their priors
    b_k <- rbeta(1, gamma, delta)
    theta_k <- rnorm(1, mu_theta, sigma_theta) # TODO: Disallow negatives?
    
    # Pad sampled parameters to their vectors
    mu <- c(mu, mu_k)
    b <- c(b, b_k)
    theta <- c(theta, theta_k)
    
    # Order parameter vectors
    order.mu <- order(mu, decreasing = T)
    mu <- mu[order.mu]
    b <- b[order.mu]
    theta <- theta[order.mu]
    
    # Add empty column to Z
    Z <- cbind(Z, matrix(0, nrow(Z), 1))
    
    # Update K-dagger
    Kdag <- max(which(colSums(Z) == 0))
    
    # Update K-star
    K.b <- Kstar
    Kstar <- max(which(mu > s)) # TODO: Is this correct?
    K.a <- Kstar
    # if (K.a != K.b) {
    #   print("ERRORRRRRRRRRRRRRRRR")
    # }
  }
  
  # Sample Z with blocked Gibbs and run FFBS on each column of Z
  # Note: Z is updated only up to column k <= K*
  # Note: After sampling Z, for any columns k <= K* that are inactive, delete those columns
  for (k in 1:Kstar) {
    # Transition matrix for appliance k
    W <- matrix(c(1-mu[k],mu[k],1-b[k],b[k]), byrow = T, nrow = 2)
    # FFBS for column k
    Z <- BSi(Y,Z,k,W,theta,sigma_epsilon,Kact,mu)
    # Update number of active appliances
    Kact <- max(which(colSums(Z) > 0))
  }
  # Delete inactive columns for k <= K*
  inact <- which(colSums(matrix(Z[,1:Kstar], nrow = length(Y))) == 0)
  if (length(inact) > 0) {
    Z <- Z[,-inact]
    mu <- mu[-inact]
    b <- b[-inact]
    theta <- theta[-inact]
  }
  
  # Update K-dagger
  Kdag <- max(which(colSums(Z) == 0))
  
  # Sample theta from its conditional
  for (k in 1:Kdag) {
    sigma_theta_p2 <- (1/sigma_theta^2 + 1/sigma_epsilon^2*sum(Z[,k]))^-1
    sum1 <- Z[,k]%*%Y
    thetadk <- theta[-k]
    if (length(thetadk)==1) {
      sum2 <- Z[,-k]*thetadk
    } else {
      sum2 <- Z[,-k]%*%thetadk
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
      ub <- 1
    } else {
      ub <- mu[k-1]
    }
    if (is.na(mu[k+1])) {
      lb <- 0
    } else {
      lb <- mu[k+1]
    }
    
    # Compute "c" constants
    ck00 <- cfun(0,0,k,Z)
    ck01 <- cfun(0,1,k,Z)
    
    # Draw mu_k using inverse grid sampling
    # Starting points for ARS
    sp <- seq(lb, ub, length.out = 15)
    sp <- sp[-c(1,15)]
    mu[k] <- ars(n = 1, bmuk, dbmuk, x = sp, m = 13, lb = T, xlb = lb,
                 ub = T, xub = ub, ck00 = ck00, ck01 = ck01)
  }
  
  # Sample mu_K-dagger with ARS
  # Starting point for ARS and upper bound
  ub <- mu[Kdag-1]
  sp <- runif(4, 0, ub)
  sp <- sp[order(sp)]
  
  # Draw mu_K-dagger
  mu[Kdag] <- ars(1, fmuk, dfmuk, x = sp, m = 4, lb = T, xlb = 0,
                  ub = T, xub = ub, alpha = alpha, t = nrow(Z))
  
  # Sample b_k
  for (k in 1:Kdag) {
    # Compute "c" constants
    ck11 <- cfun(1,1,k,Z)
    ck10 <- cfun(1,0,k,Z)
    # Draw b_k
    b[k] <- rbeta(1, ck11+gamma, ck10+delta+1)
  }
  
  # TODO: Sample hyperparameters from their posteriors (conjugacy)
  
  IterNum <- IterNum - 1
  
  # Plot observed signal and its components for this iteration
  if (iter.plot == T) {
    Z.i <- Z[,-Kdag]
    theta.i <- theta[-Kdag]
    
    try(computed.aggregate <- Z.i*theta.i)
    try(computed.aggregate <- Z.i%*%theta.i)
    plot(computed.aggregate, type = "l", ylab = "Power (W)", xlab = "Time (s)",
         ylim = c(min(min(computed.aggregate), min(theta.i)), max(max(computed.aggregate), max(theta.i))),
         col = tol5[1], main = paste(ncol(Z.i), "devices extracted")) 
    if (is.null(ncol(Z.i)) == F) {
      for (k in 1:ncol(Z.i)) {
        lines(Z.i[,k]*theta.i[k], type = "l", lty = 1, col = tol5[(k+1)%%6], lwd = 1.2)
      }
    } 
  }
  
}

# Return represenations to MIBP by removing inactive appliances
Z <- Z[,-Kdag]
theta <- theta[-Kdag]

print(paste("Retrieved", ncol(Z), "appliances."))

# Plot observed signal and its components
try(computed.aggregate <- Z*theta)
try(computed.aggregate <- Z%*%theta)
plot(computed.aggregate, type = "l", ylab = "Power (W)", xlab = "Time (s)",
     ylim = c(min(min(computed.aggregate), min(theta)), max(max(computed.aggregate), max(theta))),
     col = tol5[1], main = paste(ncol(Z), "devices extracted")) 
if (is.null(ncol(Z)) == F) {
  for (k in 1:ncol(Z)) {
    lines(Z[,k]*theta[k], type = "l", lty = 1, col = tol5[(k+1)%%6], lwd = 1.2)
  }
}

