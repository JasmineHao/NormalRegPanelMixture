library(NormalRegPanelMixture)
library(foreach)

#Generate Data
M <- 2 #Number of Type
p <- 0 #Number of Z
q <- 0 #Number of X
nrep <- 2000
cl <- makeCluster(64)

set.seed(123456)
Nset <- c(200,400)
Tset <- c(3, 5, 8)

alphaset <- list(c(0.5,0.5),c(0.2,0.8))
muset <- list(c(-1,1),c(-0.5,0.5))
sigmaset <- list(c(0.8,1.2))

anFormula.alt <- function(parlist, m, n, t){
  omega <- omega.12(parlist)
  omega <- pmin(pmax(omega, 1e-16), 0.5-1e-16)  # an becomes NaN if omega[j]=0 or 1
  omega.term <- log(omega /(0.5-omega))
  b <- c(-4.0171723345,  0.6410800569, 41.3622294442, -0.0647872609,  0.0005675749) 

  # b <-  c(-3.663692089, -0.530933177, 13.918584239, -0.135488896, -0.002936244)
  # c(-3.75007650, -0.43387385, 18.54754765, -0.13457721, -0.05850536)
  # No Penalty
  # c(-2.8171633, 0.1052204, 16.0619336, -0.1066912, -0.1597512)
  # Updated An set
  # c(-3.663692089, -0.530933177, 13.918584239, -0.135488896, -0.002936244)
  
  x <- exp(  ( log(5/95)  - b[1] - b[2]/t - b[3]/n - b[5] * omega.term ) / b[4] )  # maxa=1
  an <- 0.5 * x / (1 + x)
  an
}

GenerateSample <- function(phi,nrep){
  p = phi$p
  q = phi$q
  N = phi$N
  T = phi$T
  M = phi$M
  alpha = phi$alpha
  mu = phi$mu
  gamma = phi$gamma
  beta = phi$beta

  Data <- replicate(nrep,generateData(alpha,mu,sigma,gamma,beta,N,T,M,p,q))
  return(list(phi=phi,Data=Data))
}


PerformEMtest <- function (data, an, m = M, z = NULL, parallel) {
  library(doParallel) # workers might need information
  library(NormalRegPanelMixture)# workers might need information

  out.h0 <- normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=m,vcov.method = "none")
  out.h1 <- normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an=an,update.alpha = 1,parallel = FALSE)
  return(2 * max(out.h1$penloglik - out.h0$loglik))
}

PerformCritBoot <- function (data, an, m = M, z = NULL, parallel) {
  library(doParallel) # workers might need information
  library(NormalRegPanelMixture)# workers might need information

  crit <- regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, nbtsp = 199 ,parallel = TRUE,cl=cl)$crit
  return(crit)
}


getEstimateDiffAn <- function(Data,nrep,an,cl,M, parlist){
  lr.crit <- matrix(0.0,nr=nrep,ncol=3)
  lr.estimate.l <- matrix(0.0,nr=nrep,ncol=1)
  lr.estimate.m <- matrix(0.0,nr=nrep,ncol=1)
  lr.estimate.h <- matrix(0.0,nr=nrep,ncol=1)
  lr.size.l <- matrix(0.0,nr=nrep,ncol=1) #Nomimal size
  lr.size.m <- matrix(0.0,nr=nrep,ncol=1) #Nomimal size
  lr.size.h <- matrix(0.0,nr=nrep,ncol=1) #Nomimal size
  registerDoParallel(cl)
  results <- foreach (k = 1:nrep)%dopar% {
    data <- Data[,k]
    out.h0 <- NormalRegPanelMixture::normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none")
    out.h1.l <- NormalRegPanelMixture::normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an=(0.1 * an),update.alpha = 1,parallel = FALSE)
    2 * max(out.h1.l$penloglik - out.h0$loglik)

  }
  lr.estimate.l <- t(t(sapply(results, function(x) x[1])))

  results <- foreach (k = 1:nrep)%dopar% {
    data <- Data[,k]
    out.h0 <- NormalRegPanelMixture::normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none")
    out.h1.m <- NormalRegPanelMixture::normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an=(an),update.alpha = 1,parallel = FALSE)
    2 * max(out.h1.m$penloglik - out.h0$loglik)


  }
  lr.estimate.m <- t(t(sapply(results, function(x) x[1])))

  results <- foreach (k = 1:nrep)%dopar% {
    data <- Data[,k]
    out.h0 <- NormalRegPanelMixture::normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none")
    out.h1.h <- NormalRegPanelMixture::normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an=(10 * an) ,update.alpha = 1,parallel = FALSE)
    2 * max(out.h1.h$penloglik - out.h0$loglik)
  }
  lr.estimate.h <- t(t(sapply(results, function(x) x[1])))

  data = Data[,1]

  crit <- regpanelmixCrit(y=data$Y, x=data$X, parlist=parlist, z = data$Z,cl=cl, parallel = TRUE,nrep=1000)$crit
  for ( k in 1:nrep){
    lr.crit[k,] <- crit
  }

  for ( k in 1:nrep){
    lr.size.l[k,] <- 1 * (lr.estimate.l[k,] > lr.crit[k,2])
    lr.size.m[k,] <- 1 * (lr.estimate.m[k,] > lr.crit[k,2])
    lr.size.h[k,] <- 1 * (lr.estimate.h[k,] > lr.crit[k,2])
  }

  return(list(crit = lr.crit,nominal.size.l = apply(lr.size.l,2,mean),nominal.size.m = apply(lr.size.m,2,mean) , nominal.size.h = apply(lr.size.h,2,mean) ))
}



#GeneratePhiDataPairs
count <- 0

phi.data <- list()
nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset)


NTset <- expand.grid(Nset,Tset)
Parset <- expand.grid(muset,alphaset)
nNT <- dim(NTset)[1]
nPar <- dim(Parset)[1]
result.l <- matrix(0,nr=(nNT),nc=nPar)
rownames(result.l) <- apply(NTset,1,paste,collapse = ",")
colnames(result.l) <- apply(Parset,1,paste,collapse = ",")

result.m <- matrix(0,nr=(nNT),nc=nPar)
rownames(result.m) <- apply(NTset,1,paste,collapse = ",")
colnames(result.m) <- apply(Parset,1,paste,collapse = ",")

result.h <- matrix(0,nr=(nNT),nc=nPar)
rownames(result.h) <- apply(NTset,1,paste,collapse = ",")
colnames(result.h) <- apply(Parset,1,paste,collapse = ",")



for (r in 1:nNT){
  N <-  NTset[r,1]
  T <-  NTset[r,2]
  count <- 0
  for (mu in muset){
    for (alpha in alphaset){
      for (sigma in sigmaset){

      t <- Sys.time()
      phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = NULL, beta = NULL, N = N, T = T, M = M, p = p, q = q)

      phi.data.pair <- GenerateSample(phi,nrep)
      count <- count + 1
      Data = phi.data.pair$Data
      phi = phi.data.pair$phi
      
      if (T > 5) {
        T_an <- T 
      } else if ( T < 5) {
        T_an <- T 
      } else {
        T_an <- T
      }
      an <- anFormula.alt(phi,M,N,T_an)  #The an function according the the empirical regression
      print(anFormula(phi,M,N,T_an))
      # an <- 0.03
      print(N)
      print(T)
      print(mu)
      print(alpha)
      print(an)
      parlist = list(alpha = alpha, mubeta = mu, sigma=sigma, gam=NULL)
      result <- getEstimateDiffAn(Data,nrep,an,cl,M, parlist)


      result.l[r, count] <- result$nominal.size.l
      result.m[r, count] <- result$nominal.size.m
      result.h[r, count] <- result$nominal.size.h
      print(result$nominal.size.m)

      print(Sys.time() - t)
      }
    }
  }
}






# for(i in 1:count){
#   t.out <- Sys.time()
#   phi <- phi.data[[i]]$phi
#   Data = phi.data.pair$Data
#   an <- anFormula(phi,M,phi$N,phi$T)
#   result <- getEstimate(Data,nrep,an)
#   regression.data[i,] <- cbind(result$nominal.size,phi$N,phi$T)
#   print(cbind(result$nominal.size,phi$N,phi$T,phi$alpha,phi$mu,phi$sigma))
#   print(i)
#   print(Sys.time() - t.out)
#
# }
result.h <- result.h * 100
result.l <- result.l * 100
result.m <- result.m * 100

write.csv(result.h, file = "/home/haoyu/SizeTest/results/sizeTestM2SimH_1.csv")
write.csv(result.m, file = "/home/haoyu/SizeTest/results/sizeTestM2SimM_1.csv")
write.csv(result.l, file = "/home/haoyu/SizeTest/results/sizeTestM2SimL_1.csv")
