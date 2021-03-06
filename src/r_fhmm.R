#!/usr/bin/env Rscript
# Load libraries
library(ars)
source("src/auxiliary.R")
source("src/ffbs.R")

## Built-in help string ----

"r_fhmm.R -- an R implementation of a non-parametric factorial hidden markov
model for non-intrusive energy disaggregation research.

The script can be used interactively eg. with RStudio or via the CLI.

The CLI interface works as follows:

  Usage: r_fhmm.R [OPTS]

    -h --help               show this help
    -b --batch              operate in batch mode (eg. skips plotting)
    -i --in CSV_FILE        specify input file path
    -s --segrange RANGE     specify input data segment range
    -o --out CSV_FILE       specify output file path
    -p --sehp VAL           specify sigma_epsilon heuristic parameter value [40]
    -n --niter VAL          specify the number of iterations to run [100]
    -a --ina VAL            specify the initial number of appliances [4]
    -q --quiet              print less info
    -v --verbose            print more info

  Example:
    Rscript src/r_fhmm.R -b -i test/test.csv -o out/test.csv -s 6000:6199 -v

Formats:

  The input file must have at least two columns first of which is time and the
  last of which is the aggregate consumption data.

  The output file will have the same time column followed by a column for each
  disaggregated appliance.

Happy disaggregating!
" -> help_str

## Simple logging function ----

ts.start <- as.numeric(Sys.time())
pinfo <- function(verbosity, message) {
  if (verbosity <= mode.verbosity) {
    cat(sprintf("I%s: %s\n", formatC((as.numeric(Sys.time()) - ts.start)*1000,
        format="d", width=7, flag="0"), message))
  }
}

## Define default argument values ----

# Default run mode is interactive instead of batch
mode.batch <- F

# Default verbosity level
mode.verbosity <- 1

# Default input data file path and segment
in.path <- "test/test.csv"
in.segment <- 6000:7400

# Default output data file path
out.path <- "out/test.csv"

# Should each iteration be plotted?
iter.plot <- T

# Should negative appliances be allowed?
negative.thetas <- F

# Parameter for the sigma_epsilon heuristic
# It might be worth experimenting with this parameter
# if you get bad fit or experience other weird behaviour.
# Noisier data benefits from a lower heuristic parameter
# whereas data with low noise benefits from a higher one.
# Warning: Setting this too low might cause errors.
hp <- 10

# The initial value for the number of active appliances in the data set.
# This can be tuned if the user has some knowledge about the possible
# number of appliances. This is only the initial value, the algorithm
# will still try to extract the most probable number of appliances.
ina <- 4

# Number of iterations
IterNum.total <- 100

# Color palette by Paul Tol
tol5 <- c("#1B9E77", "#66A61E", "#7570B3", "#D95F02", "#E6AB02", "#E7298A")

## CLI argument parsing ----

args <- commandArgs(trailingOnly=TRUE)
i <- 1
while (i <= length(args)) {
  arg = args[i]
  if ((arg == '-h') || (arg == '--help')) {
    cat(help_str)
    quit(save="no", status=0, runLast=FALSE)
  } else if ((arg == '-b') || (arg == '--batch')) {
    mode.batch <- T
    iter.plot <- F
  } else if ((arg == '-i') || (arg == '--in')) {
    i <- i + 1
    in.path <- args[i]
    in.segment <- 1:1
  } else if ((arg == '-s') || (arg == '--segrange')) {
    i <- i + 1
    in.segment <- unlist(strsplit(args[i], ":"))
    in.segment <- strtoi(in.segment[1]):strtoi(in.segment[2])
  } else if ((arg == '-o') || (arg == '--out')) {
    i <- i + 1
    out.path <- args[i]
  } else if ((arg == '-p') || (arg == '--sehp')) {
    i <- i + 1
    hp <- strtoi(args[i])
  } else if ((arg == '-n') || (arg == '--niter')) {
    i <- i + 1
    IterNum.total <- strtoi(args[i])
  } else if ((arg == '-a') || (arg == '--ina')) {
    i <- i + 1
    ina <- strtoi(args[i])
  } else if ((arg == '-q') || (arg == '--quiet')) {
    mode.verbosity <- mode.verbosity - 1
  } else if ((arg == '-v') || (arg == '--verbose')) {
    mode.verbosity <- mode.verbosity + 1
  } else {
    # Abort if any arguments remain unrecognized
    cat(2, sprintf("Error: Invalid argument \"%s\"!", arg))
    quit(save="no", status=1, runLast=FALSE)
  }
  i <- i + 1
}

## Data input ----

# Read the specified segment or all of the input file into `data`
pinfo(2, sprintf("Reading input data from %s...", in.path))
in.file <- read.csv(in.path)
if (length(in.segment) > 1) {
  data <- in.file[in.segment,]
} else {
  data <- in.file
}
# Print info about the read data
pinfo(2, sprintf("- Read a %dx%d matrix!", nrow(data), ncol(data)))
if (mode.verbosity > 2) { summary(data) }

# Observed data
Y <- data[,ncol(data)]

# Add noise
#epsilon <- rnorm(length(Y), 0, 10)
#Y <- Y + abs(epsilon)

# Real values of parameters
theta.real <- c(1000,2000,500,200)
Z.real <- matrix(as.matrix(data[,c(-1,-6)]), ncol = 4)
for (col in 1:ncol(Z.real)) {
  Z.real[which(Z.real[,col] != 0),col] <- 1
}

if (mode.batch == F) {
  # Plot observed signal and its components
  pinfo(2, sprintf("Plotting ground truth..."))
  plot(Y, type = "l", ylab = "Power (W)", xlab = "Time (s)", col = tol5[1])
  lines(Z.real[,1]*theta.real[1], type = "l", lty = 1, col = tol5[2], lwd = 1.2)
  lines(Z.real[,2]*theta.real[2], type = "l", lty = 1, col = tol5[3], lwd = 1.2)
  lines(Z.real[,3]*theta.real[3], type = "l", lty = 1, col = tol5[4], lwd = 1.2)
  lines(Z.real[,4]*theta.real[4], type = "l", lty = 1, col = tol5[5], lwd = 1.2)
  legend("topleft", c("Aggregate","Oven","Heater","Fridge","Light"),
         fill = c(tol5[1],tol5[2],tol5[3],tol5[4],tol5[5]))
}

## Initialize NFHMM parameters ----

# Hyperparameters
alpha <- 1
gamma <- 1
delta <- 1
mu_theta <- mean(Y, na.rm = T)
sigma_epsilon <- median(abs(diff(Y)))*hp
sigma_theta <- sd(Y, na.rm = T)

# Need to have at least some noise
if (sigma_epsilon == 0) {
  sigma_epsilon <- 0.01*mean(Y) # "Noise" variance
}
epsilon <- rep(0, length(Y))

# Parameters
Kdag <- ina + 1 # Number of active appliances + 1
Z <- matrix(1, length(Y), (Kdag - 1)) # State matrix
Z[sample(1:length(Y), round(0.3*length(Y))),1] <- 0 # Add some zeros
Z <- cbind(Z, matrix(0, nrow(Z), 1)) # Add empty column

mu <- rbeta(Kdag, alpha/Kdag, 1) # State transition probability
mu <- mu[order(mu, decreasing = T)] # Order state trabsition probabilities
b <- rbeta(Kdag, gamma, delta) # State transition probability
theta <- rnorm(Kdag, mu_theta, sigma_theta) # State levels

## Iterative sampling for NFHMM ----

pinfo(2, sprintf("Running %d sampling iterations over %d data points...",
    IterNum.total, length(Y)))
ts.iter <- ts.start
IterNum <- IterNum.total

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
    Kstar <- max(which(mu > s)) # TODO: Is this correct?
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
    theta_cand <- rnorm(1, mu_theta_g, sigma_theta_g)
    if (negative.thetas == F) {
      while (theta_cand < 0.05*mean(Y)) {
        theta_cand <- rnorm(1, mu_theta_g, sigma_theta)
      }
    }
    theta[k] <- theta_cand
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

    # Draw mu_k using ARS
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

  # TODO: Sample hyperparameters sigma_theta, mu_theta
  # from their posteriors (conjugacy)
  alpha <- grid_sample_alpha(mu,Kdag)
  gamma <- grid_sample_gamma(b,Kdag,delta)
  delta <- grid_sample_delta(b,Kdag,gamma)

  # This is a temporary heuristic workaround for conjugacy
  m <- median(abs(diff(Y-Z%*%theta)))

  # Need to have at least some noise
  if (m == 0) {
    m = 0.01*mean(Y)
  }

  sigma_epsilon <- abs(rnorm(1, m, 0.5*m))*hp/Kact

  # Decrease iteration counter
  IterNum <- IterNum - 1

  # Plot observed signal and its components for this iteration
  if (iter.plot == T) {
    Z.i <- Z[,-Kdag]
    theta.i <- theta[-Kdag]

    try(computed.aggregate <- Z.i*theta.i, silent = T)
    try(computed.aggregate <- Z.i%*%theta.i, silent = T)
    plot(computed.aggregate, type = "l", ylab = "Power (W)", xlab = "Time (s)",
         ylim = c(min(min(computed.aggregate), min(theta.i)), max(max(computed.aggregate), max(theta.i))),
         col = tol5[1], main = paste(ncol(Z.i), "devices extracted"))
    if (is.null(ncol(Z.i)) == F) {
      for (k in 1:ncol(Z.i)) {
        lines(Z.i[,k]*theta.i[k], type = "l", lty = 1, col = tol5[(k+1)%%6], lwd = 1.2)
      }
    }
  }
  if (mode.verbosity > 1) {
    # Print a progress message at regular intervals
    ts.now <- as.numeric(Sys.time())
    if (ts.now > ts.iter + 10.0 / mode.verbosity) {
      ts.iter = ts.now
      pinfo(1, sprintf("- Iteration %d/%d finished...", IterNum.total - IterNum, IterNum.total))
    }
  }
}

# Return represenations to MIBP by removing inactive appliances
Z <- Z[,-Kdag]
theta <- theta[-Kdag]

## Print and/or plot some stats ----

# Number of extracted appliances
pinfo(1, sprintf("Iterative sampling finished: Retrieved %d appliances!", ncol(Z)))
# Mean absolute error
pinfo(2, sprintf("- MAE:   %.3f", mean(abs(Y-epsilon-Z%*%theta))))
# Mean absolute percentage error
pinfo(2, sprintf("- MAPE:  %.3f", mean(abs((Y-epsilon-Z%*%theta)/(Y-epsilon)))))

if (mode.batch == F) {
  # Plot observed signal and its components ----
  pinfo(2, sprintf("Plotting the results..."))
  try(computed.aggregate <- Z*theta, silent = T)
  try(computed.aggregate <- Z%*%theta, silent = T)
  plot(computed.aggregate, type = "l", ylab = "Power (W)", xlab = "Time (s)",
       ylim = c(min(min(computed.aggregate), min(theta)), max(max(computed.aggregate), max(theta))),
       col = tol5[1], main = paste(ncol(Z), "devices extracted"))
  if (is.null(ncol(Z)) == F) {
    for (k in 1:ncol(Z)) {
      lines(Z[,k]*theta[k], type = "l", lty = 1, col = tol5[(k+1)%%6], lwd = 1.2)
    }
  }
}

## Store results in the output file ----

pinfo(2, sprintf("Writing results to %s...", out.path))
out.table <- data.frame(matrix(NA, nrow(Z), ncol(Z)+1))
colnames(out.table) <- c("Timestamp", paste("App", seq(1:ncol(Z))))
out.table[,1] <- data[,1] # Index
for (column in 2:ncol(out.table)) {
  out.table[,column] <- Z[,column-1]*theta[column-1]
}
write.table(out.table, out.path, sep = ",", row.names = F)

# Exit
quit(save="no", status=0, runLast=FALSE)