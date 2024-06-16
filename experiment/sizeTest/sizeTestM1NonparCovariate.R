library(NormalRegPanelMixture)
library(foreach)
library(Matrix)
library(expm)


#Generate Data
M <- 1 #Number of Type
r.test <- 1 # test the null hypothesis of 2
n.grid <- 2 # partition each t into 2 intervals
p <- 0 #Number of Z
q <- 1 #Number of X
nrep <- 500
cl <- makeCluster(15)

set.seed(123456)
Nset <- c(200,400)
Tset <- c(4,6)


alphaset <- list(c(1))
muset <- list(c(-1))
betaset <- list(c(1))
sigmaset <- list(c(0.8), c(1.2))

anFormula.alt <- function(parlist, m, n, t){
  omega <- omega.12(parlist)
  omega <- pmin(pmax(omega, 1e-16), 0.5-1e-16)  # an becomes NaN if omega[j]=0 or 1
  omega.term <- log(omega /(1-omega))
  b <-   c(-0.8112790,  -0.2882271,   4.6374028,  -0.1012959,  -0.1973225)
  x <- (  b[1] + b[2]/t + b[3]/n + b[5] * omega.term ) / b[4]   # maxa=1
  an <- 1 / (1 + exp(x))
  an
}

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



NTset <- expand.grid(Nset,Tset)
Parset <- expand.grid(muset,betaset)
nNT <- dim(NTset)[1]
nPar <- dim(Parset)[1]
result <- matrix(0,nr=(nNT),nc=nPar)
nBB <- 199
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
    beta <- Parset[count,2][[1]]
    sigma <- sigmaset[[1]]
    
    print(N)
    print(T)
    print(mu)
    alpha <- 1
    
    phi = list(alpha = alpha,mu = mu, sigma = sigma, gamma = NULL, beta = beta, N = N, T = T, M = M, p = p, q = q)
    phi.data.pair <- GenerateSample(phi,nrep)
    
    Data = phi.data.pair$Data
    phi = phi.data.pair$phi
    
    
    registerDoParallel(cl)
    results <- foreach (k = 1:nrep, .packages = c("expm", "Matrix", "NormalRegPanelMixture"))%dopar% {
      #for (k in 1:nrep) { 
      data <- Data[,k]
      lm.data <- lm(as.vector(data$Y) ~ data$X)
      
      data.res <- list(Y = matrix(lm.data$residuals,nrow=T))
      
      data_P_W <- calculate_W_P(data.res, T.even, T.odd, n.grid=2, BB=199, type="indicator")
      P_c <- data_P_W$P_c 
      W_c <- data_P_W$W_c
      rk_c <- construct_stat_KP(P_c, W_c, r.test, N)
      #rk_c_all[k] <- rk_c
      #P_c_list[[k]] <- results[[k]]$P_c
      
      resample_index <- matrix(sample(N, nBB*N, replace=TRUE), nrow=nBB)
      rk_c_boot <- matrix(0,nrow=1,ncol=nBB)
      for (b in 1:nBB){
        data.res.b <- list(Y = data.res$Y[,resample_index[b,]] )
        data_P_W <- calculate_W_P(data.res.b, T.even, T.odd, n.grid=2, BB=199, type="indicator")
        P_c <- data_P_W$P_c 
        W_c <- data_P_W$W_c
        rk_c_b <- construct_stat_KP(P_c, W_c, r.test, N)
        rk_c_boot[b] <- rk_c_b
      }
      list(rk_c = rk_c, P_c = P_c, rk_c_boot=rk_c_boot)
    }
    
    rk_c_all <- matrix(0, nrow = nrep, ncol = 1)
    rk_c_crit_all <- matrix(0, nrow = nrep, ncol = 1)
    P_c_list <- vector("list", nrep)  # List to store all P_c matrices
    
    # Extract results from the list
    for (k in 1:nrep) {
      rk_c_all[k] <- results[[k]]$rk_c
      P_c_list[[k]] <- results[[k]]$P_c
      rk_c_crit_all[[k]] <- results[[k]]$crit.0.95
    }
    
    # Stop the parallel backend
    
    s_1 <- dim(P_c_list[[1]])[1]
    s_2 <- dim(P_c_list[[2]])[2]
    df <- (s_1 - r.test) * (s_2 - r.test)
    crit.0.95 <- qchisq(0.95, df)
    
    result[r, count] <- mean(rk_c_all > crit.0.95)
    result.boot[r, count] <- mean(rk_c_all > rk_c_crit_all) 
    
    print(Sys.time() - t)
    
  }
}
# stopCluster(cl)

rownames(result) <- apply(NTset,1,paste,collapse = ",")
colnames(result) <- apply(Parset,1,paste,collapse = ",")

rownames(result.boot) <- apply(NTset,1,paste,collapse = ",")
colnames(result.boot) <- apply(Parset,1,paste,collapse = ",")


write.csv(rbind(result, result.boot), file="results/sizeTestM1_nonpar_covariates.csv")
