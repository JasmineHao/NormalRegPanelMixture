library(NormalRegPanelMixture)
library(foreach)
library(Matrix)
library(expm)


#Generate Data
M <- 2 #Number of Type
r.test <- 2 # test the null hypothesis of 2
n.grid <- 3 # partition each t into 2 intervals
p <- 0 #Number of Z
q <- 0 #Number of X
nrep <- 500
cl <- makeCluster(8)
M_max <- 5  

set.seed(123456)
# Nset <- c(200,400)
# Tset <- c(4,6)

Nset <- c(200)
Tset <- c(6)

# alphaset <- list(c(0.5,0.5),c(0.2,0.8))
# muset <- list(c(-1,1),c(-0.5,0.5))
alphaset <- list(c(0.5,0.5))
muset <- list(c(-0.5,0.5))
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



count_freq <- function(stats,M_max=5){
  stats_freq_tab <- rep(0,M_max)
  for(m in 1:M_max){
    stats_freq_tab[m] = sum(stats == m)
  }
  stats_freq_tab <- stats_freq_tab/length(stats)
  return(stats_freq_tab)
}

NTset <- expand.grid(Nset,Tset)
Parset <- expand.grid(muset,alphaset)
nNT <- dim(NTset)[1]
nPar <- dim(Parset)[1]


aic_table <- matrix(0,nr=(nNT),nc=nPar)
rownames(aic_table) <- apply(NTset,1,paste,collapse = ",")
colnames(aic_table) <- apply(Parset,1,paste,collapse = ",")

bic_table <- matrix(0,nr=(nNT),nc=nPar)
rownames(bic_table) <- apply(NTset,1,paste,collapse = ",")
colnames(bic_table) <- apply(Parset,1,paste,collapse = ",")

mem_table <- matrix(0,nr=(nNT),nc=nPar)
rownames(mem_table) <- apply(NTset,1,paste,collapse = ",")
colnames(mem_table) <- apply(Parset,1,paste,collapse = ",")

mem_table_nonpar_h <- matrix(0,nr=(nNT),nc=nPar)
rownames(mem_table_nonpar_h) <- apply(NTset,1,paste,collapse = ",")
colnames(mem_table_nonpar_h) <- apply(Parset,1,paste,collapse = ",")

mem_table_nonpar_i <- matrix(0,nr=(nNT),nc=nPar)
rownames(mem_table_nonpar_i) <- apply(NTset,1,paste,collapse = ",")
colnames(mem_table_nonpar_i) <- apply(Parset,1,paste,collapse = ",")

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
      
      data <- Data[,k]
      aic <- rep(0,M_max)
      bic <- rep(0,M_max)
      lr.estim <- rep(0,M_max)
      test <- 0
      mem_result <- 0
      test.nonpar.h <- 1
      test.nonpar.i <- 1
      
      data_P_W_h <- calculate_W_P(data,T.even, T.odd, n.grid=n.grid, BB=199, type="polynomial")
      P_c_h <- data_P_W_h$P_c 
      W_c_h <- data_P_W_h$W_c
      s_1 <- dim(P_c_h)[1]
      s_2 <- dim(P_c_h)[2]
      
      data_P_W_i <- calculate_W_P(data,T.even, T.odd, n.grid=n.grid, BB=199, type="indicator")
      P_c_i <- data_P_W_i$P_c
      W_c_i <- data_P_W_i$W_c
      
      for(m in 1:M_max){
        out.h0 <- NormalRegPanelMixture::normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=m,vcov.method = "none")
        aic[m] <- out.h0$aic
        bic[m] <- out.h0$bic
        
        an <- anFormula(out.h0$parlist, m, N, T)
        out.h1 <- NormalRegPanelMixture::normalpanelmixMaxPhi(y=data$Y, parlist=out.h0$parlist,an=(an),parallel = FALSE)
        lr.estim[m] <- 2 * max(out.h1$penloglik - out.h0$loglik)
        if (test == 1){
          mem_result <- m
          suppressWarnings(
            crit <- try(NormalRegPanelMixture::regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, parallel = FALSE, nrep=1000)$crit)
          )
          if (class(crit) == "try-error"){
            crit <- NormalRegPanelMixture::regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, parallel = FALSE)$crit
          }
          if (2 * max(out.h1$penloglik - out.h0$loglik) < crit[2]){
            test = 0
          }
        }
        
        
        rk_c_h <- construct_stat_KP(P_c_h, W_c_h, m, N)
        rk_c_i <- construct_stat_KP(P_c_i, W_c_i, m, N)
        
        df <- (s_1 - m) * (s_2 - m)
        crit.0.95 <- qchisq(0.95, df)
        if (test.nonpar.h == 1){
          mem_result_nonpar_h <- m
          if (rk_c_h < crit.0.95){
            test.nonpar.h = 0
          }
        }
        if (test.nonpar.i == 1){
          mem_result_nonpar_i <- m
          if (rk_c_i < crit.0.95){
            test.nonpar.i = 0
          }
        }
      }
      
      list(lr.estim=lr.estim, aic=aic, bic=bic, mem_result=mem_result, mem_result_nonpar_h=mem_result_nonpar_h, mem_result_nonpar_i=mem_result_nonpar_i)    
    }
    
    aic <- matrix(0,nrow=nrep,ncol=M_max)
    bic <- matrix(0,nrow=nrep,ncol=M_max)
    mem_seq <-matrix(0,nrow=nrep,ncol=1)
    mem_seq_nonpar_h <-matrix(0,nrow=nrep,ncol=1)
    mem_seq_nonpar_i <-matrix(0,nrow=nrep,ncol=1)
    
    lr.estim_table <- matrix(0,nrow=nrep,ncol=M_max)
    
    for (k in 1:nrep) {
      lr.estim_table[k,] <- results[[k]]$lr.estim
      aic[k,] <- results[[k]]$aic
      bic[k,] <- results[[k]]$bic
      mem_seq[k,] <- results[[k]]$mem_result
      mem_seq_nonpar_h[k,] <- results[[k]]$mem_result_nonpar_h
      mem_seq_nonpar_i[k,] <- results[[k]]$mem_result_nonpar_i
    }
    
    
    aic_freq <- apply(aic,1,which.min)
    bic_freq <- apply(bic,1,which.min)
    
    aic_table[r, count] <-paste(count_freq(aic_freq),collapse=",")
    bic_table[r, count] <- paste(count_freq(bic_freq),collapse=",")
    mem_table[r, count] <- paste(count_freq(mem_seq),collapse=",")
    mem_table_nonpar_h[r, count] <- paste(count_freq(mem_seq_nonpar_h),collapse=",")
    mem_table_nonpar_i[r, count] <- paste(count_freq(mem_seq_nonpar_i),collapse=",")
    print(Sys.time() - t)
    
  }
}
# stopCluster(cl)
# Add a new column to each data frame to distinguish where each row comes from
aic_table$source <- "aic_table"
bic_table$source <- "bic_table"
mem_table$source <- "mem_table"
mem_table_nonpar_h$source <- "polynomial"
mem_table_nonpar_i$source <- "indicator"

write.csv(rbind(aic_table,bic_table,mem_table,mem_table_nonpar_h, mem_table_nonpar_i), file="results/seqTestM2_nonpar_polynomial.csv")
