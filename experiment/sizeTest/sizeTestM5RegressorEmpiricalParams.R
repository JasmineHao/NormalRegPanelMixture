library(NormalRegPanelMixture)
library(foreach)

M_max <- 7
#Generate Data
M <- 5 #Number of Type
p <- 0 #Number of Z
q <- 1 #Number of X
nrep <- 100
cl <- makeCluster(16)

# Based on textile industry in Chile normed log revenue  share of intermediate input against ln K, (N,T)
# 5-component model
set.seed(123456)
Nset <- c(196)
Tset <- c(3)

Parset = list(list(c(0.16076522, 0.32454077, 0.09025875, 0.35478905, 0.06964622),c(-1.241241, -0.33803875,  0.4480291,  0.52379553, 1.4139465),c(0.451833, -0.05988709, -0.2453261, -0.03106076, 0.2053708),c(0.9933480, 0.4585760, 0.9954302, 0.4116855, 0.1863346) ))


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
  X = phi$X
  Data <- replicate(nrep,generateData(alpha,mu,sigma,gamma,beta,N,T,M,p,q))
  return(list(phi=phi,Data=Data))
}



getResult <- function(Data,nrep,an,cl,M, parlist){
  registerDoParallel(cl)
  # results <- foreach (k = 1:nrep)%dopar% {
  aic_table <- matrix(0,nrow=nrep,ncol=M_max)
  bic_table <- matrix(0,nrow=nrep,ncol=M_max)
  mem_seq <-matrix(0,nrow=nrep,ncol=1)
  lr.estim_table <- matrix(0,nrow=nrep,ncol=M_max)
  
  for (k in 1:nrep){
    # library(NormalRegPanelMixture)
    t <- Sys.time()
    data <- Data[,k]
    aic <- rep(0,M_max)
    bic <- rep(0,M_max)
    lr.estim <- rep(0,M_max)
    test <- 1
    for(m in 1:M_max){
      out.h0 <- NormalRegPanelMixture::regpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=m,vcov.method = "none")
      aic[m] <- out.h0$aic
      bic[m] <- out.h0$bic
      
      out.h1 <- NormalRegPanelMixture::regpanelmixMaxPhi(y=data$Y,x=data$X, parlist=out.h0$parlist,an=(an), an_0 = 0, parallel = TRUE, cl = cl)
      lr.estim[m] <- 2 * max(out.h1$penloglik - out.h0$loglik)
      if (test == 1){
        mem_result <- m
        crit <- NormalRegPanelMixture::regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, parallel = TRUE, cl = cl)$crit
        # crit <- try(NormalRegPanelMixture::regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, parallel = TRUE, nrep=1000, cl = cl)$crit)
        if (class(crit) == "try-error"){
          crit <- NormalRegPanelMixture::regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, parallel = TRUE, cl = cl)$crit
        }
        if (2 * max(out.h1$penloglik - out.h0$loglik) < crit[2]){
          test = 0
        }
      }
    }
    lr.estim_table[k,] <- lr.estim  
    aic_table[k,] <- aic
    bic_table[k,] <- bic
    mem_seq[k,] <- mem_result
    print(paste("simulation", k, "out of", nrep))
    print(Sys.time() - t) 
  }

    aic_freq <- apply(aic_table,1,which.min)
  bic_freq <- apply(bic_table,1,which.min)
  
  return(list(aic=aic_freq,bic=bic_freq,mem_seq=mem_seq))
}



count_freq <- function(stats){
  stats_freq_tab <- rep(0,7)
  for(m in 1:7){
    stats_freq_tab[m] = sum(stats == m)
  }
  stats_freq_tab <- stats_freq_tab/length(stats)
  return(stats_freq_tab)
}

#GeneratePhiDataPairs
count <- 0

phi.data <- list()

nset <- length(Nset) * length(Tset) * length(Parset)
NTset <- expand.grid(Nset,Tset)
# Parset <- expand.grid(muset,alphaset,betaset,sigmaset)
nNT <- dim(NTset)[1]
nPar <-  length(Parset)



aic_table <- matrix(0,nr=(nNT),nc=nPar)
rownames(aic_table) <- apply(NTset,1,paste,collapse = ",")
# colnames(aic_table) <- apply(Parset,1,paste,collapse = ",")

bic_table <- matrix(0,nr=(nNT),nc=nPar)
rownames(bic_table) <- apply(NTset,1,paste,collapse = ",")
# colnames(bic_table) <- apply(Parset,1,paste,collapse = ",")

mem_table <- matrix(0,nr=(nNT),nc=nPar)
rownames(mem_table) <- apply(NTset,1,paste,collapse = ",")
# colnames(mem_table) <- apply(Parset,1,paste,collapse = ",")



for (r in 1:nNT){
  N <-  NTset[r,1]
  T <-  NTset[r,2]
  for (count in 1:nPar){
    t <- Sys.time()
    alpha = Parset[[count]][[1]]
    mu = Parset[[count]][[2]]
    beta = Parset[[count]][[3]]
    sigma = Parset[[count]][[4]]
    
    t <- Sys.time()
    phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
               beta = beta, N = N, T = T, M = M, p = p, q = q, X=NULL)
    
    phi.data.pair <- GenerateSample(phi,nrep)
    Data = phi.data.pair$Data
    phi = phi.data.pair$phi
    
    an <- anFormula(phi,M,N,T, q = 1)  #The an function according the the empirical regression
    # an <- 0.03
    print(N)
    print(T)
    print(mu)
    print(alpha)
    print(an)
    parlist = list(alpha = alpha, mubeta = mu, sigma=sigma, gam=NULL)
    result <- getResult(Data,nrep,an,cl,M, parlist)
    
    aic_freq <- result$aic
    bic_freq <- result$bic
    mem_seq  <- result$mem_seq
    aic_table[r, count] <-paste(count_freq(aic_freq),collapse=",")
    bic_table[r, count] <- paste(count_freq(bic_freq),collapse=",")
    mem_table[r, count] <- paste(count_freq(mem_seq),collapse=",")
    print(aic_table[r, count])
    
    print(Sys.time() - t)
  }
}




# write.csv(mem_table, file = "/home/haoyu/SizeTest/results/sizeTestM2_mem_table_regressor.csv")
# write.csv(aic_table, file = "/home/haoyu/SizeTest/results/sizeTestM2_aic_table_regressor.csv")
# write.csv(bic_table, file = "/home/haoyu/SizeTest/results/sizeTestM2_bic_table_regressor.csv")

write.csv(mem_table, file = "results/sizeTestM5_mem_table_regressor_empirical_param_boot.csv")
write.csv(aic_table, file = "results/sizeTestM5_aic_table_regressor_empirical_param_boot.csv")
write.csv(bic_table, file = "results/sizeTestM5_bic_table_regressor_empirical_param_boot.csv")
