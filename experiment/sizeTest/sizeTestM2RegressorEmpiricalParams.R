library(NormalRegPanelMixture)
library(doParallel)
#Generate Data
M <- 2 #Number of Type
p <- 0 #Number of Z
q <- 1 #Number of X

nrep <- 500
cl <- makeCluster(16)

set.seed(123456)
Nset <- c(200,500)
Tset <- c(2,5)

# alpha, mu, beta, sigma
# Japan chemical
c(c(0.36242787290857,0.63757212709143),c(-1.12814915135624,0.637461534227848),c(-0.163875786207446,-0.129411535097118),c(0.728355020828508,0.400456330577787) )

# Japan electronics
c(c(0.358672018033109,0.641327981966891),c(-1.06723152982067,0.59185971988529),c(0.068764706487433,0.0132637820665338),c(0.940074055069804,0.290094393239493))
# Chile food
c(c(0.183677353641936,0.816322646358064),c(-1.1182498011832,0.263661752032259),c(0.0434714956383482,0.204016320924021),c(1.3545562508164,0.660769644804908))
# Chile fabricated
c(c(0.292476796491764,0.707523203508236),c(-0.919255975118669,0.379600878757378), c(0.144001683035264, 0.0705297739448472),c(1.06910901565551,0.658356667265215))

Parset = list(list(c(0.36242787290857,0.63757212709143),c(-1.12814915135624,0.637461534227848),c(-0.163875786207446,-0.129411535097118),c(0.728355020828508,0.400456330577787) ), list(c(0.358672018033109,0.641327981966891),c(-1.06723152982067,0.59185971988529),c(0.068764706487433,0.0132637820665338),c(0.940074055069804,0.290094393239493)), list(c(0.183677353641936,0.816322646358064),c(-1.1182498011832,0.263661752032259),c(0.0434714956383482,0.204016320924021),c(1.3545562508164,0.660769644804908)), list(c(0.292476796491764,0.707523203508236),c(-0.919255975118669,0.379600878757378), c(0.144001683035264, 0.0705297739448472),c(1.06910901565551,0.658356667265215)))

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
nset <- length(Nset) * length(Tset) * length(Parset)
NTset <- expand.grid(Nset,Tset)
# Parset <- expand.grid(muset,alphaset,betaset,sigmaset)
nNT <- dim(NTset)[1]
nPar <-  length(Parset)



result.l <- matrix(0,nr=(nNT),nc=nPar)
rownames(result.l) <- apply(NTset,1,paste,collapse = ",")
# colnames(result.l) <- apply(Parset,1,paste,collapse = ",")

result.m <- matrix(0,nr=(nNT),nc=nPar)
rownames(result.m) <- apply(NTset,1,paste,collapse = ",")
# colnames(result.m) <- apply(Parset,1,paste,collapse = ",")

result.h <- matrix(0,nr=(nNT),nc=nPar)
rownames(result.h) <- apply(NTset,1,paste,collapse = ",")
# colnames(result.h) <- apply(Parset,1,paste,collapse = ",")


for (r in 1:nNT){
  N <-  NTset[r,1]
  T <-  NTset[r,2]
  print(paste(r,"/",nNT))
  
  for (count in 1:nPar){
      t <- Sys.time()
      alpha = Parset[[count]][[1]]
      mu = Parset[[count]][[2]]
      beta = Parset[[count]][[3]]
      sigma = Parset[[count]][[4]]
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

write.csv(result.m,file="results/sizeTest/sizeTestM2RegressorEmpiricalParams.csv")

