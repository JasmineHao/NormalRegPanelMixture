library(NormalRegPanelMixture)
library(foreach)
library(Matrix)
library(expm)


#Generate Data
M <- 2 #Number of Type
r.test <- 2 # test the null hypothesis of 2
n.grid <- 2 # partition each t into 2 intervals
p <- 0 #Number of Z
q <- 0 #Number of X
nrep <- 500
cl <- makeCluster(15)

set.seed(123456)
Nset <- c(200,400)
Tset <- c(5)

alphaset <- list(c(0.5,0.5),c(0.2,0.8))
muset <- list(c(-1,1),c(-0.5,0.5))
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
  sqrt_mat <- svd_decomp$u %*% diag(sqrt_singular_values) %*% t(svd_decomp$v)
  
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
    P_k_list[[k]] <- P_k
  }
  P_k_list
}




calculate_W_P_modified <- function(data, n.grid=3, BB=199){
  data_c <- data$Y
  n_size <- ncol(data_c)
  
  P_k_list <- calculate_P_matrix_modified(data_c, n.grid=3)
  Sigma_k_list <- list()
  for (k in 1:T){
    P_k <- P_k_list[[k]]
    Q_k <- P_k %*% t(P_k)
    ru <- matrix(runif(n_size * BB), nrow = n_size, ncol = BB)
    n_element <- length(as.vector(P_k))
    
    # Initialize the vec_P_b matrix
    vec_P_b <- matrix(0, nrow = BB, ncol = n_element)
    # Loop to generate vec_P_b
    for (i in 1:BB) {
      index <- ceiling(ru[, i] * n_size)
      data_b <- data$Y[,index]
      P_b_list <- calculate_P_matrix_modified(data_b, n.grid=n.grid)
      vec_P_b[i, ] <- as.vector(P_b_list[[k]])
    }
    
    # Calculate mean_vec_P_b
    mean_vec_P_b <- as.vector(P_k) # colMeans(vec_P_b)
    
    # Initialize the W_b matrix
    W_b <- matrix(0, nrow =  n_element, ncol = n_element)
    
    # Fill the W_b matrix
    for (i in 1:n_element) {
      for (j in i:n_element) {
        if (i == j) {
          W_b[i, i] <- (1 / (BB - 1)) * sum((vec_P_b[, i] - mean_vec_P_b[i])^2)
        } else {
          W_b[i, j] <- (1 / (BB - 1)) * sum((vec_P_b[, i] - mean_vec_P_b[i]) * (vec_P_b[, j] - mean_vec_P_b[j]))
          W_b[j, i] <- W_b[i, j]
        }
      }
    }
    W_c = n_size*W_b 
    Sigma_k_list[[k]] <- W_c
  }
  
  return(list(P_k_list = P_k_list, Sigma_k_list = Sigma_k_list))
}



construct_stat_KP_modified <- function(data, data_P_W, r.test){
  
  P_k_list <- data_P_W$P_k_list
  Sigma_k_list <- data_P_W$Sigma_k_list
  
  for (k in 1:T){  
    P_k <- P_k_list[[k]]
    Sigma_k <- Sigma_k_list[[k]]
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
    Lambda_km <- as.vector( t(A_rk) %*% Q_k %*% A_rk)
  }
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
  T.sequence <- 1:T
  
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
    results <- foreach (k = 1:nrep, .packages = c("expm", "Matrix", "NormalRegPanelMixture"))%dopar% {
      #for (k in 1:nrep) { 
      data <- Data[,k]
      data_P_W <- calculate_W_P_modified(data, n.grid = 3)
      
      rk_c <- construct_stat_KP_modified(data, data_P_W, r.test)
      
      list(rk_c = rk_c, P_c = P_c)
    }
    
    rk_c_all <- matrix(0, nrow = nrep, ncol = 1)
    P_c_list <- vector("list", nrep)  # List to store all P_c matrices
    
    # Extract results from the list
    for (k in 1:nrep) {
      rk_c_all[k] <- results[[k]]$rk_c
      P_c_list[[k]] <- results[[k]]$P_c
    }
    
    # Stop the parallel backend
    
    
    s_1 <- dim(P_c_list[[1]])[1]
    s_2 <- dim(P_c_list[[2]])[2]
    df <- (s_1 - r.test) * (s_2 - r.test)
    crit.0.95 <- qchisq(0.95, df)
    
    result[r, count] <- mean(rk_c_all > crit.0.95)
    print(Sys.time() - t)
    
  }
}
# stopCluster(cl)

rownames(result) <- apply(NTset,1,paste,collapse = ",")
colnames(result) <- apply(Parset,1,paste,collapse = ",")

write.csv(rbind(result), file="results/sizeTestM2_nonpar_polynomial.csv")
