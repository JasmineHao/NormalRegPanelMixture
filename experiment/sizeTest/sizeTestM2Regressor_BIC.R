library(NormalRegPanelMixture)
library(foreach)

#Generate Data
M <- 2 #Number of Type
p <- 0 #Number of Z
q <- 1 #Number of X
nrep <- 500
cl <- makeCluster(16)

set.seed(123456)
Nset <- c(200,500)
Tset <- c(5)

alphaset <- list(c(0.2,0.8))
muset <- list(c(-1,1),c(-0.5,0.5))
sigmaset <- list(c(0.1,0.1))
betaset <- list(c(1,1),c(-1,1))

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
  results <- foreach (k = 1:nrep)%dopar% {
    library(NormalRegPanelMixture)
    data <- Data[,k]
    aic <- rep(0,5)
    bic <- rep(0,5)
    lr.estim <- rep(0,5)
    test <- 1
    for(m in 1:5){
      out.h0 <- NormalRegPanelMixture::regpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=m,vcov.method = "none")
      aic[m] <- out.h0$aic
      bic[m] <- out.h0$bic
      
      out.h1 <- NormalRegPanelMixture::regpanelmixMaxPhi(y=data$Y,x=data$X, parlist=out.h0$parlist,an=(an),parallel = FALSE)
      lr.estim[m] <- 2 * max(out.h1$penloglik - out.h0$loglik)
      if (test == 1){
        mem_result <- m
        crit <- try(NormalRegPanelMixture::regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, parallel = FALSE, nrep=1000)$crit)
        if (class(crit) == "try-error"){
          crit <- NormalRegPanelMixture::regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, parallel = FALSE)$crit
        }
        if (2 * max(out.h1$penloglik - out.h0$loglik) < crit[2]){
          test = 0
        }
      }
    }
    c(lr.estim,aic,bic,mem_result)
  }
  lr.estimate <- t(t(sapply(results, function(x) x[1:5])))
  # lr.crit <- t(sapply(results, function(x) x[2:5]))
  aic <- t(sapply(results, function(x) x[6:10]))
  bic <- t(sapply(results, function(x) x[11:15]))
  mem_seq <- t(sapply(results, function(x) x[16]))
  aic_freq <- apply(aic,1,which.min)
  bic_freq <- apply(bic,1,which.min)
  return(list(aic=aic_freq,bic=bic_freq,mem_seq=mem_seq))
}



count_freq <- function(stats){
  stats_freq_tab <- rep(0,5)
  for(m in 1:5){
    stats_freq_tab[m] = sum(stats == m)
  }
  stats_freq_tab <- stats_freq_tab/length(stats)
  return(stats_freq_tab)
}

#GeneratePhiDataPairs
count <- 0

phi.data <- list()
nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset)


NTset <- expand.grid(Nset,Tset)
Parset <- expand.grid(muset,alphaset,betaset,sigmaset)
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



for (r in 1:nNT){
  N <-  NTset[r,1]
  T <-  NTset[r,2]
  for (count in 1:nPar){
    mu <- Parset[count,1][[1]]
    alpha <- Parset[count,2][[1]]
    beta <- Parset[count,3][[1]]
    sigma <- Parset[count,4][[1]]
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
    aic_table[r, count] <- paste(count_freq(aic_freq),collapse=",")
    bic_table[r, count] <- paste(count_freq(bic_freq),collapse=",")
    mem_table[r, count] <- paste(count_freq(mem_seq),collapse=",")
    print(aic_table[r, count])
    
    print(Sys.time() - t)
}
}




# write.csv(mem_table, file = "/home/haoyu/SizeTest/results/sizeTestM2_mem_table_regressor.csv")
# write.csv(aic_table, file = "/home/haoyu/SizeTest/results/sizeTestM2_aic_table_regressor.csv")
# write.csv(bic_table, file = "/home/haoyu/SizeTest/results/sizeTestM2_bic_table_regressor.csv")

write.csv(mem_table, file = "results/sizeTestM2_mem_table_regressor.csv")
write.csv(aic_table, file = "results/sizeTestM2_aic_table_regressor.csv")
write.csv(bic_table, file = "results/sizeTestM2_bic_table_regressor.csv")
