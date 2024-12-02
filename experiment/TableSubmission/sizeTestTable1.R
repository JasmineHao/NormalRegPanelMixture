# Get the directory of the current file
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

# for (N in c(200, 400, 800)){
count <- 0
result_matrix <- matrix(0, nrow = 4, ncol = 5)

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
    registerDoParallel(cl)

    # for (ii in 1:nrep){
    results <- foreach (ii = 1:nrep, .packages = c("expm", "Matrix", "NormalRegPanelMixture", "MASS"))%dopar% {
      # Record the start time
      data <- Data[, ii]
      result_rk <- compute_rk_statistics_pairwise_T(data, T.pair.list, N, M,   n.grid = 3)
      stats_KP_boot <-  construct_stat_KP_P_bootstrap(result_rk$P_k_list, result_rk$Sigma_P_list, T.pair.list, N, BB, result_rk$lambda_c, n.grid = n.grid, transform="P")
      
      out.h0 <- normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none")
      out.h1 <- normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=(M+1),vcov.method = "none")
      
      aic_rej <- out.h1$aic < out.h0$aic
      bic_rej <- out.h1$bic < out.h0$bic
      out.h1.m <- normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an=(an),parallel = FALSE)
      
      lr.m <- 2 *  max(out.h1.m$penloglik - out.h0$loglik)
      lr.crit.m <- regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, an=(an), z = data$Z, parallel = FALSE, nbtsp = 199, ninits = 10)$crit

      list(rk = result_rk$rk, rk_b = stats_KP_boot$rk_b, omega_c = result_rk$omega_c, lr.m = lr.m, lr.crit.m = lr.crit.m, aic_rej=aic_rej, bic_rej=bic_rej)
    }

    # quantile(rk_b * omega_c / ( N * diag(cov(lambda_b)) ) ,0.95)
    # qchisq(0.95, 1)

    # Calculate the duration
    end_time <- Sys.time()
    duration <- end_time - start_time
    print(duration)  # Prints the time difference

    rk_max  <- matrix(0, nrow = nrep, ncol = 1)
    rk_mean <- matrix(0, nrow = nrep, ncol = 1)
    lr_stat <- matrix(0, nrow = nrep, ncol = 1)
    rk_max.crit  <- matrix(0, nrow = nrep, ncol = 1)
    rk_mean.crit <- matrix(0, nrow = nrep, ncol = 1)
    lr_stat.crit <- matrix(0, nrow = nrep, ncol = 1)
    aic_rej.stat <- matrix(0, nrow = nrep, ncol = 1)
    bic_rej.stat <- matrix(0, nrow = nrep, ncol = 1)
    
    for (ii in 1:nrep) {
      
      rk_max[ii] <- max(results[[ii]]$rk) 
      rk_max.crit[ii] <- quantile( apply(results[[ii]]$rk_b, 1, max), 0.95)
      lr_stat[ii] <- results[[ii]]$lr.m

      rk_mean[ii] <- mean(results[[ii]]$rk) 
      rk_mean.crit[ii] <- quantile(rowMeans(results[[ii]]$rk_b), 0.95)
      lr_stat.crit[ii] <- results[[ii]]$lr.crit.m[1]
      
      aic_rej.stat[ii] <- results[[ii]]$aic_rej
      bic_rej.stat[ii] <- results[[ii]]$bic_rej
    }


    # print(mean(rk_1_Q > qchisq(0.95,1)))
    result_matrix[count, ] <- c(mean(rk_mean > rk_mean.crit),  mean(rk_max > rk_max.crit), mean(lr_stat > lr_stat.crit), mean(aic_rej.stat), mean(bic_rej.stat))
  }
}

