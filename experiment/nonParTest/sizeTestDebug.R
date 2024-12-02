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

N <- 200 
T <- 4

alpha <- c(0.5,0.5)
sigma <- c(0.8,1.2)
mu <- c(-1,1)

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




calculate_P_matrix_modified <- function(data_c, n.grid=3){
  
  T <- nrow(data_c)
  N <- ncol(data_c)
  # Create a list to store the indicator matrices
  indicator_list <- list()
  
  for (t in 1:T){
    # Calculate the quantiles
    quantiles <- quantile(data_c[t,], probs = seq(0, 1, length.out = n.grid + 1))
    quantiles[1] <- -Inf
    quantiles[n.grid+1] <- Inf
    # Create the indicator matrix
    indicator_matrix <- matrix(0, nrow = N, ncol = n.grid)
    for (n in 1:length(data_c[t,])) {
      for (m in 1:n.grid) {
        if (data_c[t,n] > quantiles[m] & data_c[t,n] <= quantiles[m + 1]) {
          indicator_matrix[n, m] <- 1
        }
      }
    }
    # Add the indicator matrix to the list
    indicator_list[[t]] <- indicator_matrix
  }
  
  # Initialize the result matrix
  k <- 1
  P_dim = length(Reduce(kronecker, lapply( (1:T)[-k], function(t) indicator_list[[t]][1, ])))
  
  # Iterate over the rows and compute the Kronecker product for each row
  P_k_list <- list()
  for (k in 1:T){
    result_matrix <- matrix(0, nrow = N, ncol = P_dim)
    for (n in 1:N) {
      kronecker_t <- Reduce(kronecker, lapply((1:T)[-k], function(t) indicator_list[[t]][n, ]))
      # Add the Kronecker product to the result matrix
      result_matrix[n, ] <- kronecker_t
    }
    P_k <- t(indicator_list[[k]]) %*% result_matrix
    P_k_list[[k]] <- P_k / N
  }
  P_k_list
}



calculate_W_P_modified <- function(data, n.grid=3){
  data_c <- data$Y
  n_size <- ncol(data_c)
  
  P_k_list <- calculate_P_matrix_modified(data_c, n.grid=n.grid)
  
  P_k <- P_k_list[[1]]
  Q_k <- P_k %*% t(P_k)
  n_element <- length(as.vector(P_k))
  n_element_q <- length(as.vector(Q_k))
  
  # init the matrix for estimating Omega / Sigma
  # vec_P_b <- matrix(0, nrow = BB, ncol = n_element * T)
  # vec_Q_b <- matrix(0, nrow = BB, ncol = n_element_q * T)
  # vec_P_b_list <- list()
  # vec_Q_b_list <- list()
  # for (k in 1:T){
  #   vec_P_b_list[[k]] <- matrix(0, nrow = BB, ncol = n_element)
  #   vec_Q_b_list[[k]] <- matrix(0, nrow = BB, ncol = n_element_q)
  # }
  
  
  # Calculate mean_vec_P_b
  Sigma_P_k_list <- list()
  Sigma_Q_k_list <- list()
  mean_vec_P <- matrix(0, nrow = 1, ncol = n_element * T)
  mean_vec_Q <- matrix(0, nrow = 1, ncol = n_element_q * T)

    
  for (k in 1:T){
    P_k <- P_k_list[[k]]
    Q_k <-P_k %*% t(P_k)
    P_k_vec <- as.vector(P_k)
    W_P_s <- (diag(P_k_vec) - P_k_vec %*% t(P_k_vec)) 
    J_P_k <- kronecker( P_k ,diag(1, nrow = n.grid)) + kronecker(diag(1, nrow = n.grid),P_k)
    W_Q_s <- J_P_k %*% W_P_s %*% t(J_P_k) #use delta method to obtain 
    Sigma_P_k_list[[k]] <- W_P_s 
    Sigma_Q_k_list[[k]] <- W_Q_s
  }
  
  return(list(P_k_list = P_k_list, Sigma_P_k_list = Sigma_P_k_list, Sigma_Q_k_list= Sigma_Q_k_list))
  
}

# data_P_W <- calculate_W_P_modified(data, n.grid=3)

construct_stat_Q_KP <- function(Q, Sigma_Q, r.test, n_size){
  # Q is the matrix for performing the test
  Q_svd <- svd(Q)
  D <- Q_svd$d
  U <- Q_svd$u
  U_12 <- U[1:(r.test),(r.test+1):ncol(U)]
  
  if (r.test == 1){
    U_12 <- t(as.matrix(U_12))
  }
  if ((dim(Q)[1] - r.test)==1){
    U_12 <- as.matrix(U_12) 
  }
  
  U_22 <- U[(r.test+1):ncol(U),(r.test+1):ncol(U)]
  svd_U22_square <- svd( U_22 %*% t(U_22))
  svd_U22 <- svd( U_22 )
  sqrtm_u_22_square <- svd_U22_square$u %*% diag(svd_U22_square$d, nrow=(dim(Q)[1] - r.test), ncol=(dim(Q)[1] - r.test)) %*% t(svd_U22_square$u)
  inv_u_22 <- svd_U22$u %*% diag(1/svd_U22$d, nrow=(dim(Q)[1] - r.test), ncol=(dim(Q)[1] - r.test)) %*% t(svd_U22$v)
  A_rk <- rbind(U_12, U_22) %*% inv_u_22 %*% sqrtm_u_22_square
  A_rk_kron <- kronecker(A_rk, A_rk)
  lambda_M <- t(A_rk) %*% Q %*% A_rk
  Omega_M <- t(A_rk_kron) %*% Sigma_Q %*% A_rk_kron

  Omega_M_pinv <- ginv(Omega_M)
  return(n_size *  t(lambda_M) %*% Omega_M_pinv %*% lambda_M)
}

construct_stat_KP <- function(P, Sigma_P, r.test, n_size){
  # Perform SVD decomposition of the matrix
  P_svd <- svd(P)
  tol_s <- 1e-10
  # Extract the U, D, and V components of the SVD decomposition
  D <- P_svd$d
  U <- P_svd$u
  V <- P_svd$v

  U_12 <- U[1:(r.test),(r.test+1):ncol(U)]
  V_12 <- V[1:(r.test),(r.test+1):ncol(U)]

  U_22 <- U[(r.test+1):ncol(U),(r.test+1):ncol(U)]
  V_22 <- V[(r.test+1):ncol(V),(r.test+1):ncol(V)]

  # Construct the A_q_o and B_q_o matrix.
  A_q_o <- t(sqrtm(U_22 %*% t(U_22)) %*% solve(t(U_22)) %*% cbind(t(U_12), t(U_22)))
  B_q_o <- sqrtm(V_22 %*% t(V_22)) %*% solve(t(V_22)) %*% cbind(t(V_12), t(V_22))
  lambda_q <- t(A_q_o) %*% P %*% t(B_q_o)
  kron_BA_o <- kronecker(B_q_o, t(A_q_o))
  Omega_q <-  kron_BA_o %*% Sigma_P %*% t(kron_BA_o)

  if (qr(Omega_q)$rank == nrow(Omega_q)) {
    r <- nrow(Omega_q)
    rk_c <- n_size * sum(as.vector(lambda_q) * solve(Omega_q) %*% as.vector(lambda_q))
  } else {
    svd_result <- svd(Omega_q)
    s <- svd_result$d
    r <- sum(s > tol_s * max(s))
    inv_Omega_q <- svd_result$v[, 1:r] %*% diag(1 / s[1:r]) %*% t(svd_result$u[, 1:r])
    rk_c <- n_size * sum( t(as.vector(lambda_q)) %*% inv_Omega_q %*% as.vector(lambda_q))
  }

  AIC_c = rk_c - 2*r
  BIC_c = rk_c - log(n_size)*r
  HQ_c  = rk_c - 2*log(log(n_size))*r
  return(rk_c)
}
  

construct_stat_KP_modified <- function(data, data_P_W, r.test){
  
  P_k_list <- data_P_W$P_k_list
  Sigma_P_k_list <- data_P_W$Sigma_P_k_list
  Sigma_Q_k_list <- data_P_W$Sigma_Q_k_list
  # Sigma_P <- data_P_W$Sigma_P
  # Sigma_Q <- data_P_W$Sigma_Q
  n_size <- ncol(data$Y)
  
  q_size <- dim(Sigma_Q_k_list[[1]])[1]
  lambda_M_k_list <- list()
  
  A_rk_kron_list <- list()
  for (k in 1:1){  
    P_k <- P_k_list[[k]] 
    Sigma_Q_k <- Sigma_Q_k_list[[k]]
    Q_k <- P_k %*% t(P_k)
    Q_svd <- svd(Q_k)
    D <- Q_svd$d
    U <- Q_svd$u
    U_12 <- U[1:(r.test),(r.test+1):ncol(U)]
    
    if (r.test == 1){
      U_12 <- t(as.matrix(U_12))
    }
    if ((dim(Q_k)[1] - r.test)==1){
      U_12 <- as.matrix(U_12) 
    }
    
    U_22 <- U[(r.test+1):ncol(U),(r.test+1):ncol(U)]
    svd_U22_square <- svd( U_22 %*% t(U_22))
    svd_U22 <- svd( U_22 )
    sqrtm_u_22_square <- svd_U22_square$u %*% diag(svd_U22_square$d, nrow=(dim(Q_k)[1] - r.test), ncol=(dim(Q_k)[1] - r.test)) %*% t(svd_U22_square$u)
    inv_u_22 <- svd_U22$u %*% diag(1/svd_U22$d, nrow=(dim(Q_k)[1] - r.test), ncol=(dim(Q_k)[1] - r.test)) %*% t(svd_U22$v)
    A_rk <- rbind(U_12, U_22) %*% inv_u_22 %*% sqrtm_u_22_square
    A_rk_kron <- kronecker(A_rk, A_rk)
    A_rk_kron_list[[k]] <- A_rk_kron
    # Omega_M_k <- t(A_rk_kron) %*% Sigma_Q_k %*% A_rk_kron
    # Omega_M_k_sqrt <- matrix_sqrt_svd(Omega_M_k)
    # Omega_M_k_sqrt_list[[k]] <- Omega_M_k_sqrt
    Lambda_M_k <- t(A_rk) %*% Q_k %*% A_rk
    lambda_M_k_list[[k]] <- as.vector(Lambda_M_k)
    Sigma_M_k <- t(A_rk_kron) %*% Sigma_Q_k  %*%  A_rk_kron
    # rk_list[[k]] <- n_size * t(lambda_M_k) %*% solve(Omega_M_k)  %*% lambda_M_k
  }
  
  # lambda_size <- dim(A_rk_kron_list[[1]])[2]
  # Omega_M <- matrix(0, nrow = T * lambda_size, ncol = T * lambda_size )
  # for (k in 1:T){
  #   for (j in 1:T){
  #     Omega_M[((k-1)*lambda_size +1):(k*lambda_size), ((j-1)*lambda_size +1):(j*lambda_size) ] <- t(A_rk_kron_list[[k]]) %*% Sigma_Q[((k-1)*q_size +1):(k*q_size), ((j-1)*q_size +1):(j*q_size)] %*% A_rk_kron_list[[j]]
  #   }
  # }
  # lambda_M <- as.matrix(unlist(lambda_M_k_list))
  # # Omega_M_sqrt_vec <- do.call(rbind, Omega_M_k_sqrt_list)
  # 
  # # Omega_M <- Omega_M_sqrt_vec %*% t(Omega_M_sqrt_vec)
  # Omega_M_pinv <- ginv(Omega_M)
  # return(n_size *  t(lambda_M) %*% Omega_M_pinv %*% lambda_M)
  return(n_size * t(Lambda_M_k) %*%  ginv(Sigma_M_k) %*% Lambda_M_k
  )
}


# Begin simulation
phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = NULL, beta = NULL, N = N, T = T, M = M, p = p, q = q)
phi.data.pair <- GenerateSample(phi,nrep)
Data = phi.data.pair$Data
phi = phi.data.pair$phi
ii <- 1
registerDoParallel(cl)
results <- foreach (ii = 1:nrep, .packages = c("expm", "Matrix", "NormalRegPanelMixture", "MASS"))%dopar% {
  #for (k in 1:nrep) { 
  data <- Data[,ii]  
  data_P_W <- calculate_W_P_modified(data, n.grid=3)  
  P_1 <- data_P_W$P_k_list[[1]]
  Q_1 <- P_1 %*% t(P_1)
  Sigma_Q_1 <- data_P_W$Sigma_Q_k_list[[1]]
  Sigma_P_1 <- data_P_W$Sigma_P_k_list[[1]]
  rk_1_q <- construct_stat_Q_KP(Q_1, Sigma_Q_1, r.test, N)
  rk_1_p <- construct_stat_KP(P_1, Sigma_P_1, r.test, N)
  
  
  # rk_c <- construct_stat_KP_modified(data, data_P_W, r.test = 2)
  list(rk_1_q = rk_1_q, rk_1_p = rk_1_p)
}

rk_p_all <- matrix(0, nrow = nrep, ncol = 1)
rk_q_all <- matrix(0, nrow = nrep, ncol = 1)
for (k in 1:nrep) {
  rk_p_all[k] <- results[[k]]$rk_1_p
  rk_q_all[k] <- results[[k]]$rk_1_q
}

quantile(rk_p_all,0.95)

qchisq(0.95,1)
