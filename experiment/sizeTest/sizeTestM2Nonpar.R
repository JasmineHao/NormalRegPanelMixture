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
cl <- makeCluster(10)

set.seed(123456)
Nset <- c(200,400)
Tset <- c(4,6)


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

NonParametricNumber <- function(data_c, T.even, T.odd, n.grid=2){
    
    T <- nrow(data_c)
    # Initialize an empty list to store the empirical breakpoints for each dimension
    empirical_breakpoints <- vector("list", T)
    
    
    # Calculate the empirical breakpoints for each dimension
    for (t in 1:T) {
      empirical_breakpoints[[t]] <- quantile(data_c[t,], probs = seq(0, 1, length.out = n.grid+1))
      empirical_breakpoints[[t]][1] <- -Inf
      empirical_breakpoints[[t]][n.grid+1] <- Inf
    }

        
    # Partition the data into total_partitions groups
    partition_indices.even <- matrix(0, nrow = ncol(data_c), ncol = length(T.even))
    partition_indices.odd <- matrix(0, nrow = ncol(data_c), ncol = length(T.odd))
    t.even <- 1
    t.odd <- 1
    
    for (t in 1:T) {
      if (t %in% T.even){
        for (j in 1:ncol(data_c)) {
          
        partition_indices.even[j, t.even] <- findInterval(data_c[t,j], empirical_breakpoints[[t]])
        }
        t.even <- t.even + 1
      }
      else{
        for (j in 1:ncol(data_c)) {
          partition_indices.odd[j, t.odd] <- findInterval(data_c[t,j], empirical_breakpoints[[t]])
        }
        t.odd <- t.odd + 1
      }
    }
    
    partition_decimal_indices.even <- apply(partition_indices.even, 1, function(x) sum((x - 1) * n.grid^(seq_along(x) - 1))) + 1
    partition_decimal_indices.odd <- apply(partition_indices.odd, 1, function(x) sum((x - 1) * n.grid^(seq_along(x) - 1))) + 1
    
    P_raw <- table(partition_decimal_indices.even, partition_decimal_indices.odd)
    
    nrow.P_c <- n.grid^(length(T.even))
    ncol.P_c <- n.grid^(length(T.odd))
    
    
    if ( nrow(P_raw) < nrow.P_c | ncol(P_raw) < ncol.P_c ){
      P_c <- matrix(0,nrow=nrow.P_c,ncol=nrow.P_c)

      for (index in rownames(P_raw)) {
        for (index.col in colnames(P_raw)){
          P_c[as.numeric(index), as.numeric(index.col)] <- P_raw[index, index.col]
        }
      }
    }
    else{
      P_c <- P_raw
    }
    P_c
    # P_raw
}
  

calculate_W_P <- function(data,T.even, T.odd, n.grid=2, BB=199){
  data_c <- data$Y
  n_size <- ncol(data_c)
  
  P_c <- NonParametricNumber(data_c, T.even, T.odd, n.grid)
  
  for (kk in 1:100){
  ru <- matrix(runif(n_size * BB), nrow = n_size, ncol = BB)
  n_element <- length(as.vector(P_c))
  # Initialize the vec_P_b matrix
  vec_P_b <- matrix(0, nrow = BB, ncol = n_element)
  
  # Loop to generate vec_P_b
  for (i in 1:BB) {
    index <- ceiling(ru[, i] * n_size)
    data_b <- data$Y[,index]
    P_b <- NonParametricNumber(data_b, T.even, T.odd, n.grid)
    vec_P_b[i, ] <- as.vector(P_b)
  }
  }
  
  # Calculate mean_vec_P_b
  mean_vec_P_b <- colMeans(vec_P_b)
  
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
  
  return(list(P_c = P_c, W_c = W_c))
}

construct_stat_KP <- function(P_c, W_c, r.test, n_size){
  # Perform SVD decomposition of the matrix
  P_svd <- svd(P_c)
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
  lambda_q <- t(A_q_o) %*% P_c %*% t(B_q_o)
  kron_BA_o <- kronecker(B_q_o, t(A_q_o))
  Omega_q <- kron_BA_o %*% W_c %*% t(kron_BA_o)
  
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
  T.even <- T.sequence[T.sequence %% 2 == 0]
  T.odd <- T.sequence[T.sequence %% 2 == 1]
  
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
      data_P_W <- calculate_W_P(data,T.even, T.odd, n.grid=n.grid, BB=199)
      P_c <- data_P_W$P_c 
      W_c <- data_P_W$W_c
      rk_c <- construct_stat_KP(P_c, W_c, r.test, N)
      #rk_c_all[k] <- rk_c
      #P_c_list[[k]] <- results[[k]]$P_c
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

write.csv(rbind(result), file="results/sizeTestM2_nonpar.csv")
