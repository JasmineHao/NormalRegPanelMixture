# Get the current file path
current_file_path <- this.path::this.path()

current_file_dir <- dirname(current_file_path)
source(file.path(current_file_dir, "functions.R"))

#Generate Data
r.test <- 2 # test the null hypothesis of 2
n.grid <- 3 # partition each t into 2 intervals
M <- 2 #Number of Type
p <- 0 #Number of Z
q <- 0 #Number of X


Nset <- c(200,400)
T <- 3
alphaset <- list(c(0.5,0.5))
muset <- list(c(-1,1),c(-0.5,0.5))
sigmaset <- list(c(0.8,1.2))
alpha <- c(0.5,0.5)
sigma <- c(0.8,1.2)


nrep <- 500
BB <- 199

cl <- makeCluster(15)

count <- 0
result_matrix <- matrix(0, nrow = 4, ncol = 4)

# N <- Nset[[1]]
# mu <- muset[[1]]
for (N in Nset){
  # for (mu in list(c(-1,1), c(-2,2))){
  for (mu in muset){
    print(mu)
    print(N)
    count <- count + 1
    start_time <- Sys.time()
    # Begin simulation
    phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = NULL, beta = NULL, N = N, T = T, M = M, p = p, q = q)
    phi.data.pair <- GenerateSample(phi,nrep)
    Data = phi.data.pair$Data
    phi = phi.data.pair$phi
    an <- anFormula(phi,M,N,T) 
    
    ii <- 1
    T.pair.list <- pairwise_combinations(1:T)
    T.triplet.list <- triplets_combinations(1:T)
    registerDoParallel(cl)
    
    # for (ii in 1:nrep){
    results <- foreach (ii = 1:nrep, .packages = c("expm", "Matrix", "NormalRegPanelMixture", "MASS"))%dopar% {
      # Record the start time
      data <- Data[, ii]
      result_rk <- compute_rk_statistics_pairwise_T(data, T.pair.list, N, M,   n.grid = 3)
      # stats_KP_boot <-  construct_stat_KP_P_bootstrap(result_rk$P_k_list, result_rk$Sigma_P_list, T.pair.list, N, BB, result_rk$lambda_c, n.grid = n.grid, transform="P")
      stats_KP_boot <- construct_stat_KP_smoothed_nonpar_bootstrap(data, T.triplet.list, N, BB, r.test, result_rk$lambda_c,  n.grid = 3)
      
      result_rk_triplet <- compute_rk_statistics_triplet_T(data, T.triplet.list, N, M,  n.grid = 3)
      # stats_KP_triplet_boot <- construct_stat_KP_P_triplet_bootstrap(result_rk_triplet$P_k_list, result_rk_triplet$Sigma_P_list, T.triplet.list, N, BB, lambda_c, n.grid = 3)
      stats_KP_triplet_boot <- construct_stat_KP_smoothed_nonpar_triplet_bootstrap(data, T.triplet.list, N, BB, r.test, result_rk_triplet$lambda_c,  n.grid = 3)
      
      list(rk = result_rk$rk, rk_b = stats_KP_boot$rk_b, omega_c = result_rk$omega_c, rk_triplet = result_rk_triplet$rk, rk_triplet_b = stats_KP_triplet_boot$rk_b)
    }
    
    # quantile(rk_b * omega_c / ( N * diag(cov(lambda_b)) ) ,0.95)
    # qchisq(0.95, 1)
    
    # Calculate the duration
    end_time <- Sys.time()
    duration <- end_time - start_time
    print(duration)  # Prints the time difference
    
    rk_max  <- matrix(0, nrow = nrep, ncol = 1)
    rk_mean <- matrix(0, nrow = nrep, ncol = 1)
    rk_max.crit  <- matrix(0, nrow = nrep, ncol = 1)
    rk_mean.crit <- matrix(0, nrow = nrep, ncol = 1)
    
    rk_max_triplet  <- matrix(0, nrow = nrep, ncol = 1)
    rk_mean_triplet <- matrix(0, nrow = nrep, ncol = 1)
    rk_max_triplet.crit  <- matrix(0, nrow = nrep, ncol = 1)
    rk_mean_triplet.crit <- matrix(0, nrow = nrep, ncol = 1)
    
    for (ii in 1:nrep) {
      
      rk_max[ii] <- max(results[[ii]]$rk) 
      rk_max.crit[ii] <- quantile( apply(results[[ii]]$rk_b, 1, max), 0.95)
      
      rk_mean[ii] <- mean(results[[ii]]$rk) 
      rk_mean.crit[ii] <- quantile(rowMeans(results[[ii]]$rk_b), 0.95)
     
      
      rk_max_triplet[ii] <- max(results[[ii]]$rk_triplet) 
      rk_max_triplet.crit[ii] <- quantile( apply(results[[ii]]$rk_triplet_b, 1, max), 0.95)
      
      rk_mean_triplet[ii] <- mean(results[[ii]]$rk_triplet) 
      rk_mean_triplet.crit[ii] <- quantile(rowMeans(results[[ii]]$rk_triplet_b), 0.95)
      
    }
    
    
    # print(mean(rk_1_Q > qchisq(0.95,1)))
    result_matrix[count, ] <- c(mean(rk_mean > rk_mean.crit),  mean(rk_max > rk_max.crit), mean(rk_mean_triplet > rk_mean_triplet.crit),  mean(rk_max_triplet > rk_max_triplet.crit))
  }
}

colnames(result_matrix) <- c('rk mean', 'rk max', 'rk mean P triplet', 'rk max P triplet')
rownames(result_matrix) <- c('200,-1,1', '200,-0.5,0.5', '400,-1,1', '400,-0.5,0.5')

write.csv(100 * result_matrix, "size_test_M2.csv", row.names = TRUE)
