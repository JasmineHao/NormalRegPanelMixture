library(NormalRegPanelMixture)
library(foreach)
library(Matrix)
library(expm)
library(MASS)

#Generate Data
M <- 1 #Number of Type
r.test <- 1 # test the null hypothesis of 2
n.grid <- 2 # partition each t into 2 intervals
p <- 0 #Number of Z
q <- 0 #Number of X
nrep <- 500
cl <- makeCluster(15)

set.seed(123456)
Nset <- c(200,400)
Tset <- c(5)

alphaset <- list(c(0.5,0.5),c(0.2,0.8))
alphaset <- list(c(1))
muset <- list(c(-1,1),c(-0.5,0.5))
muset <- list(c(0))
sigmaset <- list(c(0.8,1.2))



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


NTset <- expand.grid(Nset,Tset)
Parset <- expand.grid(muset,alphaset)
nNT <- dim(NTset)[1]
nPar <- dim(Parset)[1]

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




calculate_W_P_modified <- function(data, n.grid=3, BB=199){
  data_c <- data$Y
  n_size <- ncol(data_c)
  
  P_k_list <- calculate_P_matrix_modified(data_c, n.grid=n.grid)
  
  P_k <- P_k_list[[1]]
  Q_k <- P_k %*% t(P_k)
  n_element <- length(as.vector(P_k))
  n_element_q <- length(as.vector(Q_k))
  
  # init the matrix for estimating Omega / Sigma
  vec_P_b <- matrix(0, nrow = BB, ncol = n_element * T)
  vec_Q_b <- matrix(0, nrow = BB, ncol = n_element_q * T)
  vec_P_b_list <- list()
  vec_Q_b_list <- list()
  for (k in 1:T){
    vec_P_b_list[[k]] <- matrix(0, nrow = BB, ncol = n_element)
    vec_Q_b_list[[k]] <- matrix(0, nrow = BB, ncol = n_element_q)
  }
  
  # draw the bootstrap sample
  ru <- matrix(runif(n_size * BB), nrow = n_size, ncol = BB)
  # record the bootstrapped P and Q
  for (i in 1:BB) {
    index <- ceiling(ru[, i] * n_size)
    data_b <- data$Y[,index]
    P_b_list <- calculate_P_matrix_modified(data_b, n.grid=n.grid)
    # Loop to generate vec_P_b
    for (k in 1:T){
      P_b_k <- P_b_list[[k]]
      Q_b_k <- P_b_k%*% t(P_b_k)
      vec_P_b_list[[k]][i, ] <- as.vector(P_b_k)
      vec_Q_b_list[[k]][i, ] <- as.vector(Q_b_k)
      vec_P_b[i, ((k-1)*n_element+1):(k*n_element)] <- as.vector(P_b_k)
      vec_Q_b[i, ((k-1)*n_element_q+1):(k*n_element_q)] <- as.vector(Q_b_k)
    }
  }
  
  # Calculate mean_vec_P_b
  Sigma_P_k_list <- list()
  Sigma_Q_k_list <- list()
  mean_vec_P <- matrix(0, nrow = 1, ncol = n_element * T)
  mean_vec_Q <- matrix(0, nrow = 1, ncol = n_element_q * T)
  for (k in 1:T){
    P_k <- P_k_list[[k]]
    Q_k <-P_k %*% t(P_k)
    mean_vec_P_b <- as.vector(P_k) # colMeans(vec_P_b)
    mean_vec_Q_b <- as.vector(Q_k) # colMeans(vec_P_b)
    mean_vec_P[1, ((k-1)*n_element+1):(k*n_element)] <- mean_vec_P_b
    mean_vec_Q[1, ((k-1)*n_element_q+1):(k*n_element_q)] <- mean_vec_Q_b
    
    diff_q <- sweep(vec_Q_b_list[[k]], 2, mean_vec_Q_b, FUN = "-")
    W_Q_b <- n_size * ( t(diff_q) %*% diff_q ) / (BB - 1)
    
    diff_p <- sweep(vec_P_b_list[[k]], 2, mean_vec_P_b, FUN = "-")
    W_P_b <- n_size * ( t(diff_p) %*% diff_p ) / (BB - 1)
    
    Sigma_P_k_list[[k]] <- W_P_b
    Sigma_Q_k_list[[k]] <- W_Q_b
    
  }
  
  diff_q <- sweep(vec_Q_b, 2, mean_vec_Q, FUN = "-")
  diff_p <- sweep(vec_P_b, 2, mean_vec_P, FUN = "-")
  Sigma_Q <- n_size * ( t(diff_q) %*% diff_q ) / (BB - 1)
  Sigma_P <- n_size * ( t(diff_p) %*% diff_p ) / (BB - 1)
    
  return(list(P_k_list = P_k_list, Sigma_P_k_list = Sigma_P_k_list, Sigma_Q_k_list= Sigma_Q_k_list, Sigma_P=Sigma_P, Sigma_Q=Sigma_Q))
    
}
  
# data_P_W <- calculate_W_P_modified(data, n.grid=3, BB=199)

construct_stat_KP_modified <- function(data, data_P_W, r.test){
  
  P_k_list <- data_P_W$P_k_list
  Sigma_P_k_list <- data_P_W$Sigma_P_k_list
  Sigma_Q_k_list <- data_P_W$Sigma_Q_k_list
  Sigma_P <- data_P_W$Sigma_P
  Sigma_Q <- data_P_W$Sigma_Q
  n_size <- ncol(data$Y)

  q_size <- dim(Sigma_Q_k_list[[1]])[1]
  lambda_M_k_list <- list()
  #Omega_M_k_sqrt_list <- list()
  A_rk_kron_list <- list()
  for (k in 1:T){  
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
    # rk_list[[k]] <- n_size * t(lambda_M_k) %*% solve(Omega_M_k)  %*% lambda_M_k
  }
  
  lambda_size <- dim(A_rk_kron_list[[1]])[2]
  Omega_M <- matrix(0, nrow = T * lambda_size, ncol = T * lambda_size )
  for (k in 1:T){
    for (j in 1:T){
      Omega_M[((k-1)*lambda_size +1):(k*lambda_size), ((j-1)*lambda_size +1):(j*lambda_size) ] <- t(A_rk_kron_list[[k]]) %*% Sigma_Q[((k-1)*q_size +1):(k*q_size), ((j-1)*q_size +1):(j*q_size)] %*% A_rk_kron_list[[j]]
    }
  }
  lambda_M <- as.matrix(unlist(lambda_M_k_list))
  # Omega_M_sqrt_vec <- do.call(rbind, Omega_M_k_sqrt_list)
  
  # Omega_M <- Omega_M_sqrt_vec %*% t(Omega_M_sqrt_vec)
  Omega_M_pinv <- ginv(Omega_M)
  return(n_size *  t(lambda_M) %*% Omega_M_pinv %*% lambda_M)
}



NTset <- expand.grid(Nset,Tset)
Parset <- expand.grid(muset,alphaset)
nNT <- dim(NTset)[1]
nPar <- dim(Parset)[1]
result <- matrix(0,nr=(nNT),nc=nPar)

for (r in 1:nNT){
  N <-  NTset[r,1]
  T <-  NTset[r,2]
  
  # Create a sequence from 1 to T
  
  # Partition the sequence into even and odd numbers
  # T.even <- T.sequence[T.sequence %% 2 == 0]
  # T.odd <- T.sequence[T.sequence %% 2 == 1]
  for (count in 1:nPar){
    t <- Sys.time()
    mu <- Parset[count,1][[1]]
    alpha <- Parset[count,2][[1]]
    sigma <- sigmaset[[1]]
    
    print(N)
    print(T)
    print(mu)
    print(alpha)
    
    
    phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = NULL, beta = NULL, N = N, T = T, M = M, p = p, q = q)
    phi.data.pair <- GenerateSample(phi,nrep)
    
    Data = phi.data.pair$Data
    phi = phi.data.pair$phi
    
    
    registerDoParallel(cl)
    results <- foreach (ii = 1:nrep, .packages = c("expm", "Matrix", "NormalRegPanelMixture", "MASS"))%dopar% {
      #for (k in 1:nrep) { 
      data <- Data[,ii]
      data_P_W <- calculate_W_P_modified(data, n.grid = 3)
      rk_c <- construct_stat_KP_modified(data, data_P_W, r.test = 2)
      
      list(rk_c = rk_c)
    }
    
    rk_c_all <- matrix(0, nrow = nrep, ncol = 1)
    # P_c_list <- vector("list", nrep)  # List to store all P_c matrices
    
    # Extract results from the list
    for (k in 1:nrep) {
      rk_c_all[k] <- results[[k]]$rk_c
      #   P_c_list[[k]] <- results[[k]]$P_c
    }
    
    # Stop the parallel backend
    quantile_95 <- qchisq(0.95, df = T)
    
    
    # s_1 <- dim(P_c_list[[1]])[1]
    # s_2 <- dim(P_c_list[[2]])[2]
    # df <- (s_1 - r.test) * (s_2 - r.test)
    # crit.0.95 <- qchisq(0.95, df)
    
    result[r, count] <- mean(rk_c_all > quantile_95)
    print(Sys.time() - t)
    
  }
}
# stopCluster(cl)

rownames(result) <- apply(NTset,1,paste,collapse = ",")
colnames(result) <- apply(Parset,1,paste,collapse = ",")

write.csv(rbind(result), file="results/sizeTestM2_nonpar_polynomial.csv")
