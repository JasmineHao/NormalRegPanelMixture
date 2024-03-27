library(NormalRegPanelMixture)
library(doParallel)
#Generate Data
M <- 2 #Number of Type
p <- 0 #Number of Z
q <- 1 #Number of X

nrep <- 500
cl <- makeCluster(12)

set.seed(123456)
Nset <- c(200,400)
Tset <- c(3, 5, 8)
Tset <- c(2)
alphaset <- list(c(0.5,0.5),c(0.2,0.8))
muset <- list(c(-1,1),c(-0.5,0.5))
sigmaset <- list(c(0.8,1.2))
betaset <- list(c(1,1))

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


getEstimateDiffAn <- function(Data,nrep,an,cl,M){
  registerDoParallel(cl)
  results <- foreach (k = 1:nrep)%dopar% {
    library(NormalRegPanelMixture)
    data <- Data[,k]
    out.h0 <- NormalRegPanelMixture::regpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none")
    out.h1.l <- NormalRegPanelMixture::regpanelmixMaxPhi(y=data$Y,x=data$X, parlist=out.h0$parlist,an=(0.1 * an),parallel = FALSE)
    out.h1.m <- NormalRegPanelMixture::regpanelmixMaxPhi(y=data$Y,x=data$X, parlist=out.h0$parlist,an=(an),parallel = FALSE)
    out.h1.h <- NormalRegPanelMixture::regpanelmixMaxPhi(y=data$Y,x=data$X, parlist=out.h0$parlist,an=(10 * an) ,parallel = FALSE)
    
    crit <- try(NormalRegPanelMixture::regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, parallel = FALSE, nrep=1000)$crit)
    if (class(crit) == "try-error"){
      crit <- NormalRegPanelMixture::regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, parallel = FALSE)$crit
    }
    c(2 * max(out.h1.l$penloglik - out.h0$loglik), 2 * max(out.h1.m$penloglik - out.h0$loglik), 2 * max(out.h1.h$penloglik - out.h0$loglik), crit)
    
  }
  
  lr.estimate.l <- t(t(sapply(results, function(x) x[1])))
  lr.estimate.m <- t(t(sapply(results, function(x) x[2])))
  lr.estimate.h <- t(t(sapply(results, function(x) x[3])))
  lr.crit <- t(sapply(results, function(x) x[4:length(x)]))
  
  lr.size.l <- matrix(0.0,nr=nrep,ncol=1) #Nomimal size
  lr.size.m <- matrix(0.0,nr=nrep,ncol=1) #Nomimal size
  lr.size.h <- matrix(0.0,nr=nrep,ncol=1) #Nomimal size
  
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
nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset) * length(betaset)

NTset <- expand.grid(Nset,Tset)
Parset <- expand.grid(muset,alphaset,betaset,sigmaset)
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
  print(paste(r,"/",nNT))
  count <- 0
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
    an <- anFormula(phi,M,N,T,q=1) #The an function according the the empirical regression
    print(an)
    result <- getEstimateDiffAn(Data,nrep,an,cl,M)
    result.l[r, count] <- result$nominal.size.l
    result.m[r, count] <- result$nominal.size.m
    result.h[r, count] <- result$nominal.size.h
    print(result$nominal.size.m)
    print(Sys.time() - t)
    
  }
  
}
result.m <- result.m * 100
result.l <- result.l * 100
result.h <- result.h * 100


write.csv(rbind(result.h,result.m,result.l), file="results/sizeTestM2Regressorm_HML_T2.csv")



