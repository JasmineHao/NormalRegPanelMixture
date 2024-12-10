# Get the directory of the current file
# Install this.path if not already installed
if (!requireNamespace("this.path", quietly = TRUE)) {
  install.packages("this.path")
}

# Get the current file path
current_file_path <- this.path::this.path()

current_file_dir <- dirname(current_file_path)
source(file.path(current_file_dir, "functions.R"))

#Generate Data
M <- 2 #Number of Type
r.test <- 2 # test the null hypothesis of 2
n.grid.X <- 2 
n.grid <- 3 # partition each t into 2 intervals
p <- 0 #Number of Z
q <- 1 #Number of X
nrep <- 100
cl <- makeCluster(15)


N <- 200 
T <- 3

alpha <- c(0.5,0.5)
sigma <- c(0.8,1.2)
mu <- c(-1,1)
beta <- c(1,1)

phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
               beta = beta, N = N, T = T, M = M, p = p, q = q, X=NULL)
phi.data.pair <- GenerateSample(phi,nrep)
Data = phi.data.pair$Data
phi = phi.data.pair$phi
an <- anFormula(phi,M,N,T,q=1) #The an function according the the empirical         

T.triplet.list <- triplets_combinations(1:T)

# For each monte carlo simulation
ii <- 1
data <- Data[,ii]


data$Y 
data$X 

X.matrix <- matrix(data$X, nrow=T, ncol=N)


indicator_list.X <- lapply(1:T, function(t) {
  # Calculate the quantiles
  quantiles <- quantile(X.matrix[t,], probs = seq(0, 1, length.out = n.grid.X + 1))
  quantiles[1] <- -Inf
  quantiles[n.grid.X + 1] <- Inf
  
  # Use vectorized operation to create the indicator matrix
  cut_indices <- cut(X.matrix[t,], breaks = quantiles, labels = FALSE, include.lowest = TRUE)
  indicator_matrix <- matrix(0, nrow = N, ncol = n.grid.X)
  indicator_matrix[cbind(1:N, cut_indices)] <- 1
  
  return(indicator_matrix)
})  

# Generate a list based on combinations of indicator_list.X[[1]] and indicator_list.X[[3]]
X_combinations_list <- lapply(1:(n.grid.X * n.grid.X), function(combo) {
  # Decode the combination index into the respective bins for X[[1]] and X[[3]]
  bin1 <- (combo - 1) %/% n.grid.X + 1  # Bin for X[[1]]
  bin3 <- (combo - 1) %% n.grid.X + 1  # Bin for X[[3]]
  
  # Find individuals belonging to this combination
  indices <- which(indicator_list.X[[1]][, bin1] == 1 &
                     indicator_list.X[[3]][, bin3] == 1)
  
  return(indices)
})



k <- 1
T.triplet <- T.triplet.list[[k]]
result_matrix <- t(sapply(1:N, function(n) Reduce(kronecker, lapply(T.triplet[-1], function(t) indicator_list.Y[[t]][n, ]))))
P_k <- t(indicator_list.Y.ngrid[[triplet_k[1]]]) %*% result_matrix


