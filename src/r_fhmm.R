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
K <- 3 # Number of appliances
Z <- matrix(0, nrow(data), K) # State matrix
mu <- matrix(0, 1, K) # State transition probability
b <- matrix(0, 1, K) # State transition probability
theta <- matrix(0, 1, K) # State levels

# Iterative sampling for NFHMM ----
# Number of iterations
IterNum <- 1
while (IterNum > 0) {
  IterNum <- IterNum - 1
}
