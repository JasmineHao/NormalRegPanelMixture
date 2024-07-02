library(NormalRegPanelMixture)
library(foreach)
library(Matrix)
library(expm)


#Generate Data
M <- 2 #Number of Type
r.test <- 1 # test the null hypothesis of 1 component
n.grid <- 2 # partition each t into 2 intervals
p <- 0 #Number of Z
q <- 3 #Number of X
p.0 <- round(( p - 1 ) / 2)
q.0 <- round(( q - 1 ) / 2) 


nrep <- 500
cl <- makeCluster(12)

set.seed(123456)
Nset <- c(200,400)
Tset <- c(4,6)


alphaset <- list(c(0.5,0.5))

rho <- c(0.5,0.5)
beta.r <- c(-1,1)
mu.r <- c(-1,1)

muset <- list(mu.r * (1-rho))
mu0set <- list(mu.r / sqrt(1- rho**2))


betaset <- list( t(matrix(rbind(rho, beta.r, -rho * beta.r),nrow = q)) )
beta0set <- list(beta.r / sqrt(1 - rho**2))


sigma <- c(0.75,0.75)
sigma0 <- sqrt(sigma**2 / (1 - rho**2))


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
  p.0 = phi$p.0
  q.0 = phi$q.0
  
  N = phi$N
  T = phi$T
  M = phi$M
  alpha = phi$alpha
  sigma = phi$sigma
  sigma0 = phi$sigma0
  mu = phi$mu
  mu0 = phi$mu0
  gamma = phi$gamma
  beta = phi$beta
  

  Data <- replicate(nrep,generateDataAR1(
    alpha,mu,sigma,gamma,beta,mu0,sigma0,gamma0,
    beta0,N, T ,M,p,q,p.0,q.0))
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
result.boot <- matrix(0,nr=(nNT),nc=nPar)

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
    mu <- Parset[count,1][[count]]
    beta <- Parset[count,2][[count]]
    sigma <- sigmaset[[count]]
    
    mu0 <- mu0set[[count]]
    beta0 <- beta0set[[count]]
    sigma0 <- sigma0set[[count]]
    print(N)
    print(T)
    print(mu)
    alpha <- alphaset[[1]]
    
    phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = NULL, beta = beta, 
               mu0 = mu0, sigma0 = sigma0 , N = N, T = T, M = M, p = p, q = q, p.0 = p.0, q.0 = q.0)
    phi.0 <- list(alpha = alpha, gamma = NULL, beta = NULL, 
                  mu = mu0, sigma = sigma0)
    
    phi.data.pair <- GenerateSample(phi,nrep)
    
    Data = phi.data.pair$Data
    phi = phi.data.pair$phi
    
    
    registerDoParallel(cl)
    results <- foreach (k = 1:nrep, .packages = c("expm", "Matrix", "NormalRegPanelMixture"))%dopar% {
      #for (k in 1:nrep) { 
      data <- Data[,k]
      lm.data <- lm(as.vector(data$Y) ~ data$X )
      
      data.res <- list(Y = matrix(lm.data$residuals,nrow=T))
      
      data_P_W <- calculate_W_P(data.res, T.even, T.odd, n.grid=2, BB=199, type="indicator")
      P_c <- data_P_W$P_c 
      W_c <- data_P_W$W_c
      rk_c <- construct_stat_KP(P_c, W_c, r.test, N)
      #rk_c_all[k] <- rk_c
      #P_c_list[[k]] <- results[[k]]$P_c
      
      resample_index <- matrix(sample(N, nBB*N, replace=TRUE), nrow=nBB)
      rk_c_boot <- matrix(0,nrow=1,ncol=nBB)
      # for (b in 1:nBB){
      #   data.res.b <- list(Y = data.res$Y[,resample_index[b,]] )
      #   data_P_W <- calculate_W_P(data.res.b, T.even, T.odd, n.grid=2, BB=199, type="indicator")
      #   P_c <- data_P_W$P_c
      #   W_c <- data_P_W$W_c
      #   rk_c_b <- construct_stat_KP(P_c, W_c, r.test, N)
      #   rk_c_boot[b] <- rk_c_b
      # }
      list(rk_c = rk_c, P_c = P_c, crit.0.95=quantile(rk_c_boot, 0.95))
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
    # result.boot[r, count] <- mean(rk_c_all > rk_c_crit_all) 
    print(Sys.time() - t)
    
  }
}
# stopCluster(cl)


rownames(result) <- apply(NTset,1,paste,collapse = ",")
colnames(result) <- apply(Parset,1,paste,collapse = ",")

rownames(result.boot) <- apply(NTset,1,paste,collapse = ",")
colnames(result.boot) <- apply(Parset,1,paste,collapse = ",")

write.csv(rbind(result, result.boot), file="results/powerTestM1_nonpar_AR1.csv")
