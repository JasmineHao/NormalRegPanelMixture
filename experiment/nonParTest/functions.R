library(NormalRegPanelMixture)
library(foreach)
library(Matrix)
library(expm)
library(MASS)


GenerateSample <- function(phi,nrep){
  p = phi$p
  q = phi$q
  N = phi$N
  T = phi$T
  M = phi$M
  alpha = phi$alpha
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
  P_svd <- svd(P, nu=nrow(P),nv=ncol(P))
  
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
  # if ( (transform == "P") & (dim(P)[1] != dim(P)[2] ) ){
  #   stop("If we construct the P stats, the P matrix should be square matrix.")
  # }

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

# Define the function
compute_rk_statistics_pairwise_T <- function(data, T.pair.list, N, r.test,   n.grid = 3) {
  # Initialize weights
  weights <- rep(1 / N, N)
  
  # Compute P and Sigma matrices based on T.pair.list
  data_P_W <- calculate_P_matrix_t_pair(data$Y, T.pair.list, weights, n.grid = n.grid)
  
  # Initialize vectors and lists to store results
  rk <- numeric(length(T.pair.list))
  lambda_c <- numeric(length(T.pair.list))
  omega_c <- numeric(length(T.pair.list))
  
  Sigma_P_list <- list()
  P_k_list <- list()
  
  # Loop through T.pair.list to compute statistics
  for (k in 1:length(T.pair.list)) {
    # Extract P_k and Sigma_P_k from the data_P_W object
    P_k <- data_P_W$P_k_list[[k]]
    Sigma_P_k <- data_P_W$Sigma_P_k_list[[k]]
    
    # Compute KP statistics for the k-th pair
    stat_KP <- construct_stat_KP(P_k, Sigma_P_k, r.test, N, transform = "P")
    
    # Store results
    rk[k] <- stat_KP$rk_c
    lambda_c[k] <- stat_KP$lambda_c
    omega_c[k] <- stat_KP$Omega_q
    Sigma_P_list[[k]] <- Sigma_P_k
    P_k_list[[k]] <- P_k
  }
  
  # Return results as a list
  return(list(
    rk = rk,
    lambda_c = lambda_c,
    omega_c = omega_c,
    Sigma_P_list = Sigma_P_list,
    P_k_list = P_k_list
  ))
}

# for test purpose: compare triplet P and pairwise P
# -------------------------------------------------------
calculate_P_matrix_t_triplet <- function(data_c, T.triplet.list, weights, n.grid = n.grid){
  # Optimized code
  T <- nrow(data_c)
  N <- ncol(data_c)
  
  indicator_list.Y <- lapply(1:T, function(t) {
    # Calculate the quantiles
    quantiles <- quantile(data_c[t,], probs = seq(0, 1, length.out = 3))
    quantiles[1] <- -Inf
    quantiles[3] <- Inf
    
    # Use vectorized operation to create the indicator matrix
    cut_indices <- cut(data_c[t,], breaks = quantiles, labels = FALSE, include.lowest = TRUE)
    indicator_matrix <- matrix(0, nrow = N, ncol = 2)
    indicator_matrix[cbind(1:N, cut_indices)] <- 1
    
    return(indicator_matrix)
  })  
  
  indicator_list.Y.ngrid <- lapply(1:T, function(t) {
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
  for (k in 1:length(T.triplet.list)){
    T.triplet <- T.triplet.list[[k]]
    result_matrix <- t(sapply(1:N, function(n) Reduce(kronecker, lapply(T.triplet[-1], function(t) indicator_list.Y[[t]][n, ]))))
    P_k <- t( weights * indicator_list.Y.ngrid[[T.triplet[1]]] ) %*% result_matrix
    P_k_list[[k]] <- P_k
    
    P_k_vec <- as.vector(P_k)
    W_P_s <- (diag(P_k_vec) - P_k_vec %*% t(P_k_vec)) 
    Sigma_P_k_list[[k]] <- W_P_s
  }
  
  return(list(P_k_list=P_k_list,Sigma_P_k_list = Sigma_P_k_list))
}


# Define the function
compute_rk_statistics_triplet_T <- function(data, T.triplet.list, N, r.test, n.grid = 3) {
  # Initialize weights
  weights <- rep(1 / N, N)
  
  # Compute P and Sigma matrices based on T.pair.list
  data_P_W <- calculate_P_matrix_t_triplet(data$Y, T.triplet.list, weights, n.grid = n.grid)
  
  # Initialize vectors and lists to store results
  rk <- numeric(length(T.triplet.list))
  lambda_c <- list()
  omega_c  <- list()
  Sigma_P_list <- list()
  P_k_list <- list()
  
  # Loop through T.pair.list to compute statistics
  for (k in 1:length(T.triplet.list)) {
    # Extract P_k and Sigma_P_k from the data_P_W object
    P_k <- data_P_W$P_k_list[[k]]
    Sigma_P_k <- data_P_W$Sigma_P_k_list[[k]]
    
    # Compute KP statistics for the k-th pair
    stat_KP <- construct_stat_KP(P_k, Sigma_P_k, r.test, N, transform = "P")
    
    # Store results
    rk[k] <- stat_KP$rk_c
    lambda_c[[k]] <- stat_KP$lambda_c
    omega_c[[k]] <- stat_KP$Omega_q
    Sigma_P_list[[k]] <- Sigma_P_k
    P_k_list[[k]] <- P_k
  }
  
  # Return results as a list
  return(list(
    rk = rk,
    lambda_c = lambda_c,
    omega_c = omega_c,
    Sigma_P_list = Sigma_P_list,
    P_k_list = P_k_list
  ))
}



construct_stat_KP_P_triplet_bootstrap  <- function(P_k_list, Sigma_P_list, T.triplet.list, N, BB, lambda_c, n.grid = 3, transform="P") {
  # Initialize result matrices
  rk_b <- matrix(0, nrow = BB, ncol = length(T.triplet.list))
  lambda_b <- matrix(0, nrow = BB, ncol = length(T.triplet.list))
  omega_b <- matrix(0, nrow = BB, ncol = length(T.triplet.list))
  
  # Bootstrap loop
  for (k in 1:length(T.triplet.list)) {
    
    Sigma_k_0 <- Sigma_P_list[[k]]
    vec_k_0 <- P_k_list[[k]]
  
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
    }
  }
  
  # Return results as a list
  return(list(
    rk_b = rk_b
  ))
}




construct_stat_KP_smoothed_nonpar_triplet_bootstrap  <- function(data, T.triplet.list, N, BB, r.test, lambda_c,  n.grid = 3, transform="P") {
  # Generate the ru matrix (random weights normalized by row sum)
  ru <- matrix(rexp(N * BB, rate = 1), nrow = BB, ncol = N)
  ru <- apply(ru, 1, function(row) row / sum(row))
  
  # Initialize result matrices
  rk_b <- matrix(0, nrow = BB, ncol = length(T.triplet.list))
  lambda_b <- matrix(0, nrow = BB, ncol = length(T.triplet.list))
  omega_b <- matrix(0, nrow = BB, ncol = length(T.triplet.list))
  
  # Bootstrap loop
  for (i in 1:BB) {
    # Calculate bootstrapped P and Q matrices
    data_P_W_b <- calculate_P_matrix_t_triplet(data$Y, T.triplet.list, ru[, i], n.grid = n.grid)
    
    # Loop over T.triplet.list
    for (k in 1:length(T.triplet.list)) {
      P_k <- data_P_W_b$P_k_list[[k]]
      Sigma_P_k <- data_P_W_b$Sigma_P_k_list[[k]]
      
      # Compute KP statistics for the k-th pair
      stat_KP <- construct_stat_KP(P_k, Sigma_P_k, r.test, N, lambda_c = lambda_c[[k]], transform = "P")
      
      # Update result matrices
      rk_b[i, k] <- stat_KP$rk_c
    }
  }
  
  # Return results as a list
  return(list(
    rk_b = rk_b  ))
}
