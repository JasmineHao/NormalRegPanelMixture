library(NormalRegPanelMixture)
library(foreach)

#Generate Data
M <- 2 #Number of Type
p <- 0 #Number of Z
q <- 0 #Number of X
k_max <- 5
nrep <- 2000
cl <- makeCluster(detectCores()-1)

set.seed(123456)
Nset <- c(200,400)
Tset <- c(3, 5, 8)

alphaset <- list(c(0.5,0.5),c(0.2,0.8))
muset <- list(c(-1,1),c(-0.5,0.5))
sigmaset <- list(c(0.8,1.2))

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


PerformEMtest <- function (data, an, m = M, z = NULL, parallel) {
  library(doParallel) # workers might need information
  library(NormalRegPanelMixture)# workers might need information
  
  out.h0 <- normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=m,vcov.method = "none")
  out.h1 <- normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an=an,parallel = FALSE)
  return(2 * max(out.h1$penloglik - out.h0$loglik))
}


getEstimate <- function(Data,nrep,an,cl,M, parlist){
  registerDoParallel(cl)
  results <- foreach (k = 1:nrep)%dopar% {
    library(NormalRegPanelMixture)
    data <- Data[,k]
    out.h0 <- NormalRegPanelMixture::normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none")
    out.h1 <- NormalRegPanelMixture::normalpanelmixPMLE.M1(y=data$Y,x=data$X,z = data$Z, parlist=out.h0$parlist,an=an)
    
    crit <- try(NormalRegPanelMixture::regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, parallel = FALSE, nrep=1000)$crit)
    if (class(crit) == "try-error"){
      crit <- NormalRegPanelMixture::regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, an = an, z = data$Z, parallel = FALSE)$crit
    }
    c(2 * max(out.h1$penloglik - out.h0$loglik), crit)
  }
  
  lr.estimate <- t(t(sapply(results, function(x) x[1])))
  lr.crit <- t(sapply(results, function(x) x[2:length(x)]))
  lr.size <- matrix(0.0,nr=nrep,ncol=1) #Nomimal size
  for ( k in 1:nrep){
    lr.size[k,] <- 1 * (lr.estimate[k,] > lr.crit[k,2])
  }
  
  return(list(crit = lr.crit,nominal.size = apply(lr.size,2,mean),
              lr.estimate = lr.estimate ))
}




#GeneratePhiDataPairs
phi.data <- list()
nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset)


NTset <- expand.grid(Nset,Tset)
Parset <- expand.grid(muset,alphaset)
nNT <- dim(NTset)[1]
nPar <- dim(Parset)[1]

result.PLR <- matrix(0,nr=(nNT),nc=nPar)
rownames(result.PLR) <- apply(NTset,1,paste,collapse = ",")
colnames(result.PLR) <- apply(Parset,1,paste,collapse = ",")


for (r in 1:nNT){
  N <-  NTset[r,1]
  T <-  NTset[r,2]
  for (count in 1:nPar){
    mu <- Parset[count,1][[1]]
    alpha <- Parset[count,2][[1]]
    sigma <- sigmaset[[1]]
    
    t <- Sys.time()
    phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = NULL, beta = NULL, N = N, T = T, M = M, p = p, q = q)
    
    phi.data.pair <- GenerateSample(phi,nrep)
    Data = phi.data.pair$Data
    phi = phi.data.pair$phi
    
    if (T > 5) {
      T_an <- T 
    } else if ( T < 5) {
      T_an <- T 
    } else {
      T_an <- T
    }
    an <- 10 * anFormula.alt(phi,M,N,T_an)  #The an function according the the empirical regression
    print(paste('N is' ,N))
    print(paste('T is' ,T))
    print(paste('mu is' , paste(mu)))
    print(paste('alpha is', paste(alpha)))
    print(paste('an is',anFormula(phi,M,N,T_an)))
    parlist = list(alpha = alpha, mubeta = mu, sigma=sigma, gam=NULL)
    result <- getEstimate(Data,nrep,an,cl,M, parlist)
    
    result.PLR[r, count] <- result$nominal.size
    print(result$nominal.size)
    
    print(Sys.time() - t)
    
  }
}




result.PLR <- result.PLR * 100

write.csv(result.PLR, file = paste("results/sizeTestM2SimPLR.csv", sep="_"))
