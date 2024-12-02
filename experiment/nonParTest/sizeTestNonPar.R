library(NormalRegPanelMixture)
library(foreach)
library(Matrix)
library(expm)
library(MASS)

#Generate Data
M <- 2 #Number of Type
r.test <- 2 # test the null hypothesis of 2
n.grid <- 3 # partition each t into 2 intervals
p <- 0 #Number of Z
q <- 0 #Number of X
nrep <- 500
cl <- makeCluster(15)

GenerateSample <- function(phi,nrep){
  p = phi$p
  q = phi$q
  N = phi$N
  T = phi$T
  M = phi$M
  alpha = phi$alpha
  sigma = phi$sigma
  mu = phi$mu
  gamma = phi$gamma
  beta = phi$beta
  
  Data <- replicate(nrep,generateData(alpha,mu,sigma,gamma,beta,N,T,M,p,q))
  return(list(phi=phi,Data=Data))
}

matrix_sqrt_svd <- function(mat) {
  if (!is.matrix(mat)) {
    stop("Input must be a matrix.")
  }
  
  # Check if the matrix is square (number of rows equals the number of columns)
  if (nrow(mat) != ncol(mat)) {
    stop("Input must be a square matrix.")
  }
  
  # Compute the SVD of the matrix
  svd_decomp <- svd(mat)
  
  # Compute the square root of the singular values
  sqrt_singular_values <- sqrt(svd_decomp$d)
  
  # Construct the square root matrix using the SVD components
  sqrt_mat <- svd_decomp$u %*% diag(sqrt_singular_values,nrow=nrow(mat),ncol=ncol(mat)) %*% t(svd_decomp$v)
  
  return(sqrt_mat)
}

# Function to return all pairwise combinations of elements in a list or array
pairwise_combinations <- function(input_array) {
  
  # Get the length of the input
  n <- length(input_array)
  
  # Initialize an empty list to store combinations
  combinations <- list()
  
  # Counter for indexing the combinations
  counter <- 1
  
  # Loop through all pairwise combinations
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # Create a 2-element array for each pair
      combinations[[counter]] <- c(input_array[[i]], input_array[[j]])
      counter <- counter + 1
    }
  }
  
  return(combinations)
}


# Function to find all triplets in a given list
triplets_combinations <- function(input_array) {
  # Check if the input has at least 3 elements
  if (length(input_array) < 3) {
    stop("The input list must have at least 3 elements to form triplets.")
  }
  
  # Generate all combinations of 3 elements
  triplets_unordered <- combn(input_array, 3, simplify = FALSE)
  triplets <- list()
  for (i in seq_along(triplets_unordered)){
    triplet <- triplets_unordered[[i]]
    triplets[[(i-1)*3 + 1]] <- c(triplet[1], triplet[2], triplet[3])  # First element is 1
    triplets[[(i-1)*3 + 2]] <- c(triplet[2], triplet[1], triplet[3])  # First element is 2
    triplets[[(i-1)*3 + 3]] <- c(triplet[3], triplet[1], triplet[2])  # First element is 3    
  }
  
  # Return the list of triplets
  return(triplets)
}


calculate_P_matrix_t_pair <- function(data_c, T.pair.list, weights, n.grid=3){
  # T <- nrow(data_c)
  # N <- ncol(data_c)
  # Create a list to store the indicator matrices
  
  # Optimized code
  indicator_list <- lapply(1:T, function(t) {
    # Calculate the quantiles
    quantiles <- quantile(data_c[t,], probs = seq(0, 1, length.out = n.grid + 1))
    quantiles[1] <- -Inf
    quantiles[n.grid + 1] <- Inf
    
    # Use vectorized operation to create the indicator matrix
    cut_indices <- cut(data_c[t,], breaks = quantiles, labels = FALSE, include.lowest = TRUE)
    indicator_matrix <- matrix(0, nrow = N, ncol = n.grid)
    indicator_matrix[cbind(1:N, cut_indices)] <- 1
    
    return(indicator_matrix)
  })  
  # Initialize the result matrix
  P_k_list <- list() 
  Sigma_P_k_list <- list()
  # Iterate over the rows and compute the Kronecker product for each row
  for (k in 1:length(T.pair.list)){
    T.pair <- T.pair.list[[k]]
    P_k <- t(weights * indicator_list[[T.pair[1]]]) %*% indicator_list[[T.pair[2]]]    
    P_k_list[[k]] <- P_k 

    P_k_vec <- as.vector(P_k)
    W_P_s <- (diag(P_k_vec) - P_k_vec %*% t(P_k_vec)) 
    Sigma_P_k_list[[k]] <- W_P_s
  }
  return(list(P_k_list=P_k_list,Sigma_P_k_list = Sigma_P_k_list))
}


calculate_P_matrix_t_triplet <- function(data_c, T.triplets.list, weights, n.grid=3){
  T <- nrow(data_c)
  N <- ncol(data_c)
  # Create a list to store the indicator matrices
  # Optimized code
  indicator_list <- lapply(1:T, function(t) {
    # Calculate the quantiles
    quantiles <- quantile(data_c[t,], probs = seq(0, 1, length.out = n.grid + 1))
    quantiles[1] <- -Inf
    quantiles[n.grid + 1] <- Inf
    
    # Use vectorized operation to create the indicator matrix
    cut_indices <- cut(data_c[t,], breaks = quantiles, labels = FALSE, include.lowest = TRUE)
    indicator_matrix <- matrix(0, nrow = N, ncol = n.grid)
    indicator_matrix[cbind(1:N, cut_indices)] <- 1
    
    return(indicator_matrix)
  })  
  # Initialize the result matrix
  P_k_list <- list() 
  Sigma_P_k_list <- list()
  # Iterate over the rows and compute the Kronecker product for each row
  for (k in 1:length(T.triplets.list)){
    T.triplets <- T.triplets.list[[k]]
    # Initialize the result matrix
    result_matrix <- do.call(rbind, lapply(1:N, function(n) {
      # Compute the Kronecker product for each n
      Reduce(kronecker, lapply(T.triplets[-1], function(t) indicator_list[[t]][n, ]))
    }))

    # Compute P_k without the loop
    P_k <- t(indicator_list[[1]]) %*% result_matrix / N
    P_k_list[[k]] <- P_k 

    P_k_vec <- as.vector(P_k)
    W_P_s <- (diag(P_k_vec) - P_k_vec %*% t(P_k_vec)) 
    Sigma_P_k_list[[k]] <- W_P_s
  }
  return(list(P_k_list=P_k_list,Sigma_P_k_list = Sigma_P_k_list))
}

calculate_P_matrix_modified <- function(data_c, n.grid=3){
  
  T <- nrow(data_c)
  N <- ncol(data_c)
  # Create a list to store the indicator matrices
  indicator_list <- lapply(1:T, function(t) {
    # Calculate the quantiles
    quantiles <- quantile(data_c[t,], probs = seq(0, 1, length.out = n.grid + 1))
    quantiles[1] <- -Inf
    quantiles[n.grid + 1] <- Inf
    
    # Use vectorized operation to create the indicator matrix
    cut_indices <- cut(data_c[t,], breaks = quantiles, labels = FALSE, include.lowest = TRUE)
    indicator_matrix <- matrix(0, nrow = N, ncol = n.grid)
    indicator_matrix[cbind(1:N, cut_indices)] <- 1
    
    return(indicator_matrix)
  })
   
  # Initialize the result matrix
  k <- 1
  P_dim = length(Reduce(kronecker, lapply( (1:T)[-k], function(t) indicator_list[[t]][1, ])))
  
  # Iterate over the rows and compute the Kronecker product for each row
  P_k_list <- list()
  Q_k_list <- list()
  Sigma_P_k_list <- list()
  Sigma_Q_k_list <- list()
  for (k in 1:T){
    result_matrix <- matrix(0, nrow = N, ncol = P_dim)
    for (n in 1:N) {
      kronecker_t <- Reduce(kronecker, lapply((1:T)[-k], function(t) indicator_list[[t]][n, ]))
      # Add the Kronecker product to the result matrix
      result_matrix[n, ] <- kronecker_t
    }
    P_k <- t(indicator_list[[k]]) %*% result_matrix / N
    P_k_list[[k]] <- P_k
    Q_k_list[[k]] <- P_k %*% t(P_k) 
    P_k_vec <- as.vector(P_k)
    W_P_s <- (diag(P_k_vec) - P_k_vec %*% t(P_k_vec)) 
    
    J_P_k <- kronecker( P_k ,diag(1, nrow = n.grid)) 
    W_Q_s <- J_P_k %*% W_P_s %*% t(J_P_k) #use delta method to obtain 
    Sigma_P_k_list[[k]] <- W_P_s
    Sigma_Q_k_list[[k]] <- W_Q_s
    
  }
  list(P_k_list = P_k_list, Q_k_list = Q_k_list, Sigma_P_k_list=Sigma_P_k_list, Sigma_Q_k_list=Sigma_Q_k_list)
}

invert_matrix <- function(mat, epsilon = 1e-8) {
  # Check if the matrix is square
  if (!is.matrix(mat) || nrow(mat) != ncol(mat)) {
    stop("The input must be a square matrix.")
  }
  
  # Check if the matrix is singular (determinant close to zero)
  det_val <- det(mat)
  if (abs(det_val) < epsilon) {
    # cat("Matrix is singular or nearly singular. Regularizing by adding epsilon to the diagonal.\n")
    # Add epsilon to the diagonal elements
    mat <- mat + diag(epsilon, nrow(mat))
  }
  
  # Return the inverse of the (possibly regularized) matrix
  return(solve(mat))
}

# Define the function
matrix_svd_decomposition <- function(P, r.test) {
  # Perform Singular Value Decomposition (SVD) on matrix P
  P_svd <- svd(P)
  
  # Extract the singular values and matrices
  D <- P_svd$d
  U <- P_svd$u
  V <- P_svd$v
  
  # Submatrices of U and V
  U_12 <- U[1:r.test, (r.test + 1):ncol(U)]
  V_12 <- V[1:r.test, (r.test + 1):ncol(V)]
  
  U_22 <- U[(r.test + 1):ncol(U), (r.test + 1):ncol(U)]
  V_22 <- V[(r.test + 1):ncol(V), (r.test + 1):ncol(V)]
  
  # Construct the A_q_o and B_q_o matrices
  A_q_o <- t(sqrtm(U_22 %*% t(U_22)) %*% invert_matrix(t(U_22)) %*% cbind(t(U_12), t(U_22)))
  B_q_o <- sqrtm(V_22 %*% t(V_22)) %*% invert_matrix(t(V_22)) %*% cbind(t(V_12), t(V_22))
    
  # Compute the Kronecker product of B_q_o and t(A_q_o)
  kron_BA_o <- kronecker(B_q_o, t(A_q_o))
  
  # Return all results as a list
  return(list(
    D = D,
    U = U,
    V = V,
    U_12 = U_12,
    V_12 = V_12,
    U_22 = U_22,
    V_22 = V_22,
    A_q_o = A_q_o,
    B_q_o = B_q_o,
    kron_BA_o = kron_BA_o
  ))
}

construct_stat_KP <- function(P, Sigma_P, r.test, n_size, lambda_c=0, transform="P"){
  allowed_strings <- c("P","Q")
  # Check if the input string is in the allowed list
  if (!(transform %in% allowed_strings)) {
    stop(paste("Invalid input! The input_string must be one of:", 
               paste(allowed_strings, collapse = ", ")))
  }
  if ( (transform == "P") & (dim(P)[1] != dim(P)[2] ) ){
    stop("If we construct the P stats, the P matrix should be square matrix.")
  }

  if (transform == "P"){
    # Perform SVD decomposition of the matrix
    P_svd <- matrix_svd_decomposition(P, r.test)
    lambda_q <- t(P_svd$A_q_o) %*% P %*% t(P_svd$B_q_o) - lambda_c
    Omega_q <-  P_svd$kron_BA_o %*% Sigma_P %*% t(P_svd$kron_BA_o)

    # if (qr(Omega_q)$rank == nrow(Omega_q)) {
    r <- nrow(Omega_q)
    rk_c <- n_size * sum(as.vector(lambda_q) * invert_matrix(Omega_q) %*% as.vector(lambda_q))
  } else {
    Q <- P %*% t(P)
    J_P <- kronecker(P ,diag(1, nrow = n.grid)) 
    Sigma_Q <- J_P %*% Sigma_P %*% t(J_P)

    Q_svd <- matrix_svd_decomposition(Q, r.test)
    lambda_q <- t(Q_svd$A_q_o) %*% Q %*% t(Q_svd$B_q_o) - lambda_c
    Omega_q <- Q_svd$kron_BA_o %*% Sigma_Q %*% t(Q_svd$kron_BA_o)

    r <- nrow(Omega_q)
    rk_c <- n_size * sum(as.vector(lambda_q) * invert_matrix(Omega_q) %*% as.vector(lambda_q))
  }

  AIC_c = rk_c - 2*r
  BIC_c = rk_c - log(n_size)*r
  HQ_c  = rk_c - 2*log(log(n_size))*r

  return(list(rk_c = rk_c, lambda_c=lambda_q, Omega_q = Omega_q))
}

construct_stat_KP_smoothed_nonpar_bootstrap  <- function(data, T.pair.list, N, BB, r.test, lambda_c,  n.grid = 3, transform="P") {
  # Generate the ru matrix (random weights normalized by row sum)
  ru <- matrix(rexp(N * BB, rate = 1), nrow = BB, ncol = N)
  ru <- apply(ru, 1, function(row) row / sum(row))
  
  # Initialize result matrices
  rk_b <- matrix(0, nrow = BB, ncol = length(T.pair.list))
  lambda_b <- matrix(0, nrow = BB, ncol = length(T.pair.list))
  omega_b <- matrix(0, nrow = BB, ncol = length(T.pair.list))
  
  # Bootstrap loop
  for (i in 1:BB) {
    # Calculate bootstrapped P and Q matrices
    data_P_W_b <- calculate_P_matrix_t_pair(data$Y, T.pair.list, ru[, i], n.grid = n.grid)
    
    # Loop over T.pair.list
    for (k in 1:length(T.pair.list)) {
      P_k <- data_P_W_b$P_k_list[[k]]
      Sigma_P_k <- data_P_W_b$Sigma_P_k_list[[k]]
      
      # Compute KP statistics for the k-th pair
      stat_KP <- construct_stat_KP(P_k, Sigma_P_k, r.test, N, lambda_c = lambda_c[k], transform = transform)
      
      # Update result matrices
      rk_b[i, k] <- stat_KP$rk_c
      lambda_b[i, k] <- stat_KP$lambda_c
      omega_b[i, k] <- stat_KP$Omega_q
    }
  }
  
  # Return results as a list
  return(list(
    rk_b = rk_b,
    lambda_b = lambda_b,
    omega_b = omega_b
  ))
}


construct_stat_KP_P_bootstrap  <- function(P_k_list, Sigma_P_list, T.pair.list, N, BB, lambda_c, n.grid = 3, transform="P") {
  # Initialize result matrices
  rk_b <- matrix(0, nrow = BB, ncol = length(T.pair.list))
  lambda_b <- matrix(0, nrow = BB, ncol = length(T.pair.list))
  omega_b <- matrix(0, nrow = BB, ncol = length(T.pair.list))
  
  # Bootstrap loop
  for (k in 1:length(T.pair.list)) {
    if (transform == "Q"){ 
      vec_k_0  <- P_k_list[[k]] %*% t(P_k_list[[k]])
      J_P <- kronecker(P_k_list[[k]] ,diag(1, nrow = n.grid))
      Sigma_k_0 <- J_P %*% Sigma_P_list[[k]] %*% t(J_P)
      
    } else{
      Sigma_k_0 <- Sigma_P_list[[k]]
      vec_k_0 <- P_k_list[[k]]
    }
    vec_BB <- mvrnorm(n = BB, mu = as.vector(vec_k_0), Sigma = Sigma_k_0 /  N )
    
    # mvrnorm(n = 1, mu = as.vector(P_k_0), Sigma = Sigma_P_k / 1e10
    # )
    for (i in 1:BB) {
      # Loop over T.pair.list
      vec_b <- matrix(vec_BB[i,], nrow = nrow(vec_k_0), ncol = ncol(vec_k_0)) 
      # Compute KP statistics for the k-th pair
      stat_KP <- construct_stat_KP(vec_b, Sigma_k_0, r.test, N, lambda_c = lambda_c[[k]], transform = "P")
      
      # Update result matrices
      rk_b[i, k] <- stat_KP$rk_c
      lambda_b[i, k] <- stat_KP$lambda_c
      omega_b[i, k] <- stat_KP$Omega_q
    }
  }
  
  # Return results as a list
  return(list(
    rk_b = rk_b,
    lambda_b = lambda_b,
    omega_b = omega_b
  ))
}



set.seed(123)  # Set the random seed
T <- 2
alpha <- c(0.5,0.5)
sigma <- c(0.8,1.2)
mu <- c(-1,1)
N <- 200
transform="P"
# for (N in c(200, 400, 800)){
for (N in c(200)){
  # for (mu in list(c(-1,1), c(-2,2))){
  for (mu in list(c(-1,1))){
    print(mu)
    print(N)
    start_time <- Sys.time()
    # Begin simulation
    phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = NULL, beta = NULL, N = N, T = T, M = M, p = p, q = q)
    phi.data.pair <- GenerateSample(phi,nrep)
    Data = phi.data.pair$Data
    phi = phi.data.pair$phi
    
    ii <- 1
    BB <- 199
    T.pair.list <- pairwise_combinations(1:T)
    registerDoParallel(cl)

    # for (ii in 1:nrep){
    results <- foreach (ii = 1:nrep, .packages = c("expm", "Matrix", "NormalRegPanelMixture", "MASS"))%dopar% {
      # Record the start time
      data <- Data[, ii]
      weights <- rep(1/N, N)
      data_P_W <- calculate_P_matrix_t_pair(data$Y, T.pair.list, weights, n.grid=3) 
      # data_P_W_modified <- calculate_P_matrix_modified(data$Y, n.grid = 3)

      rk <- numeric(length(T.pair.list))
      lambda_c <- numeric(length(T.pair.list))
      omega_c <- numeric(length(T.pair.list))
      
      rk_Q <- numeric(length(T.pair.list))
      lambda_c_Q <- numeric(length(T.pair.list))
      Sigma_P_list <- list()
      P_k_list <- list()
      for (k in 1:length(T.pair.list)){
        # For P stats
        P_k <- data_P_W$P_k_list[[k]]
        Sigma_P_k <- data_P_W$Sigma_P_k_list[[k]]
        stat_KP <- construct_stat_KP(P_k, Sigma_P_k, r.test, N, transform=transform)
        # stat_KP <-  
        rk[k] <- stat_KP$rk_c
        lambda_c[k] <- stat_KP$lambda_c
        omega_c[k] <- stat_KP$Omega_q
        Sigma_P_list[[k]] <- Sigma_P_k
        P_k_list[[k]] <- P_k
      }
      # stats_KP_boot <- construct_stat_KP_smoothed_nonpar_bootstrap(data, T.pair.list, N, BB, r.test, lambda_c,  n.grid = 3, transform=transform)
      stats_KP_boot <-  construct_stat_KP_P_bootstrap(P_k_list, Sigma_P_list, T.pair.list, N, BB, lambda_c, n.grid = n.grid, transform=transform)
       
      quantile(stats_KP_boot$rk_b, 0.95)

      # Check, these should be of the same scale
      # tmp <- sweep(lambda_b, 2, sqrt(omega_c / (N * diag(cov(lambda_b)))), "*")
      # print(N * diag(cov(tmp)))
      # print(N * diag(cov(lambda_b)))
      # omega_c
      # colMeans(omega_b)
      # Record the end time
      # list(rk = rk, rk_b = sweep(
        # rk_b, 2, omega_c / (N * diag(cov(lambda_b))), "*"
      # ), omega_c = omega_c)
      list(rk = rk, rk_b = stats_KP_boot$rk_b, omega_c = omega_c, rk_Q=rk_Q)
    }

    # quantile(rk_b * omega_c / ( N * diag(cov(lambda_b)) ) ,0.95)
    # qchisq(0.95, 1)

    # Calculate the duration
    end_time <- Sys.time()
    duration <- end_time - start_time
    print(duration)  # Prints the time difference

    rk_1  <- matrix(0, nrow = nrep, ncol = 1)
    rk_1_Q  <- matrix(0, nrow = nrep, ncol = 1)

    rk_max  <- matrix(0, nrow = nrep, ncol = 1)
    rk_mean <- matrix(0, nrow = nrep, ncol = 1)

    rk_1.crit  <- matrix(0, nrow = nrep, ncol = 1)
    rk_max.crit  <- matrix(0, nrow = nrep, ncol = 1)
    rk_mean.crit <- matrix(0, nrow = nrep, ncol = 1)

    for (ii in 1:nrep) {
      rk_1[ii] <- results[[ii]]$rk[1]
      rk_1_Q[ii] <- results[[ii]]$rk_Q[1]
      
      rk_max[ii] <- max(results[[ii]]$rk) 
      rk_max.crit[ii] <- quantile( apply(results[[ii]]$rk_b, 1, max), 0.95)
      rk_1.crit[ii] <- quantile(results[[ii]]$rk_b[,1], 0.95)
      
      rk_mean[ii] <- mean(results[[ii]]$rk) 
      rk_mean.crit[ii] <- quantile(rowMeans(results[[ii]]$rk_b), 0.95)
    }

    print(mean(rk_1 > qchisq(0.95,1)))
    # print(mean(rk_1_Q > qchisq(0.95,1)))
    print(mean(rk_1 > rk_1.crit))
    print(mean(rk_mean > rk_mean.crit))
    print(mean(rk_max > rk_max.crit))
  }
}


# Use triplets Q stats
set.seed(123)  # Set the random seed
T <- 3
alpha <- c(0.5,0.5)
sigma <- c(0.8,1.2)
mu <- c(-1,1)
N <- 200
transform="Q"
# for (N in c(200, 400, 800)){
for (N in c(200)){
  # for (mu in list(c(-1,1), c(-2,2))){
  for (mu in list(c(-1,1))){
    print(mu)
    print(N)
    start_time <- Sys.time()
    # Begin simulation
    phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = NULL, beta = NULL, N = N, T = T, M = M, p = p, q = q)
    phi.data.pair <- GenerateSample(phi,nrep)
    Data = phi.data.pair$Data
    phi = phi.data.pair$phi
    
    ii <- 1
    BB <- 199
    # T.triplets.list <- triplets_combinations(1:T)
    T.triplets.list <- list(c(1,2,3))
    registerDoParallel(cl)

    # for (ii in 1:nrep){
    results <- foreach (ii = 1:nrep, .packages = c("expm", "Matrix", "NormalRegPanelMixture", "MASS"))%dopar% {
      # Record the start time
      data <- Data[, ii]
      weights <- rep(1/N, N)
      data_P_W <- calculate_P_matrix_t_triplet(data$Y, T.triplets.list, weights, n.grid=3)
      
      # data_P_W_modified <- calculate_P_matrix_modified(data$Y, n.grid = 3)
      rk <- numeric(length(T.triplets.list))
      lambda_c <- numeric(length(T.triplets.list))
      omega_c <- numeric(length(T.triplets.list))
      
      rk_Q <- numeric(length(T.triplets.list))
      lambda_c_Q <- numeric(length(T.triplets.list))
      Sigma_P_list <- list()
      P_k_list <- list()
      for (k in 1:length(T.triplets.list)){
        # For P stats
        P_k <- data_P_W$P_k_list[[k]]
        
        Sigma_P_k <- data_P_W$Sigma_P_k_list[[k]]
        stat_KP <- construct_stat_KP(P_k, Sigma_P_k, r.test, N, lambda_c = 0, transform=transform)
        
        # stat_KP <-  
        rk[k] <- stat_KP$rk_c
        lambda_c[k] <- stat_KP$lambda_c
        omega_c[k] <- stat_KP$Omega_q
        Sigma_P_list[[k]] <- Sigma_P_k
        P_k_list[[k]] <- P_k
      }
      # stats_KP_boot <- construct_stat_KP_smoothed_nonpar_bootstrap(data, T.pair.list, N, BB, r.test, lambda_c,  n.grid = 3, transform=transform)
      stats_KP_boot <- construct_stat_KP_P_bootstrap(P_k_list, Sigma_P_list, T.triplets.list, N, BB, lambda_c, n.grid = n.grid, transform=transform)
       
      # quantile(stats_KP_boot$rk_b[,1], 0.95)

      list(rk = rk, rk_b = stats_KP_boot$rk_b, omega_c = omega_c, rk_Q=rk_Q)
    }

    # quantile(rk_b * omega_c / ( N * diag(cov(lambda_b)) ) ,0.95)
    # qchisq(0.95, 1)

    # Calculate the duration
    end_time <- Sys.time()
    duration <- end_time - start_time
    print(duration)  # Prints the time difference

    rk_1  <- matrix(0, nrow = nrep, ncol = 1)
    rk_1_Q  <- matrix(0, nrow = nrep, ncol = 1)

    rk_max  <- matrix(0, nrow = nrep, ncol = 1)
    rk_mean <- matrix(0, nrow = nrep, ncol = 1)

    rk_1.crit  <- matrix(0, nrow = nrep, ncol = 1)
    rk_max.crit  <- matrix(0, nrow = nrep, ncol = 1)
    rk_mean.crit <- matrix(0, nrow = nrep, ncol = 1)

    for (ii in 1:nrep) {
      rk_1[ii] <- results[[ii]]$rk[1]
      rk_1_Q[ii] <- results[[ii]]$rk_Q[1]
      
      rk_max[ii] <- max(results[[ii]]$rk) 
      rk_max.crit[ii] <- quantile( apply(results[[ii]]$rk_b, 1, max), 0.95)
      rk_1.crit[ii] <- quantile(results[[ii]]$rk_b[,1], 0.95)
      
      rk_mean[ii] <- mean(results[[ii]]$rk) 
      rk_mean.crit[ii] <- quantile(rowMeans(results[[ii]]$rk_b), 0.95)
    }

    print(mean(rk_1 > qchisq(0.95,1)))
    # print(mean(rk_1_Q > qchisq(0.95,1)))
    print(mean(rk_1/4 > rk_1.crit ))
    print(mean(rk_mean > rk_mean.crit))
    print(mean(rk_max > rk_max.crit))
  }
}
