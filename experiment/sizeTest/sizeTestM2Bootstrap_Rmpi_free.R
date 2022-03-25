# Generate table 3 and 4
library(NormalRegPanelMixture)
library(doParallel)

#Generate Data
M <- 2 #Number of Type
p <- 0 #Number of Z
q <- 0 #Number of X

set.seed(123456)
Nset <- c(100,500)
Tset <- c(2,5,10)
alphaset <- list(c(0.5,0.5),c(0.2,0.8))
muset <- list(c(-1,1),c(-0.5,0.5))
sigma <- c(0.8,1.2)


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
  # if (q != 0){
  #   parlist <- list('alpha' = alpha, 
  #                   'mubeta' = t(cbind(mu,beta)), 
  #                   'sigma' = sigma, 'gam' = gamma)
  # }else{
  #   parlist <- list('alpha' = alpha, 'mubeta' = mu, 'sigma' = sigma, 'gam' = gamma)
  # }
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
  
  crit <- regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, nbtsp = 999 ,parallel = TRUE,cl=cl)$crit
  return(crit)
}


getEstimateDiffAn <- function(Data,nrep,an,cl){
  lr.crit.l <- matrix(0.0,nr=nrep,ncol=3)
  lr.crit.m <- matrix(0.0,nr=nrep,ncol=3)
  lr.crit.h <- matrix(0.0,nr=nrep,ncol=3)
  lr.estimate <- matrix(0.0,nr=nrep,ncol=1)
  lr.size.l <- matrix(0.0,nr=nrep,ncol=1) #Nomimal size
  lr.size.m <- matrix(0.0,nr=nrep,ncol=1) #Nomimal size
  lr.size.h <- matrix(0.0,nr=nrep,ncol=1) #Nomimal size

  for (k in 1:nrep){

    data <- Data[,k]
    out.h0 <- normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none")
    out.h1 <- normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an= an ,update.alpha = 1,parallel = TRUE,cl=cl)
    # out.h1.h <- normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an= 10 *　an ,update.alpha = 1,parallel = TRUE,cl=cl)
    lr.estimate[k,] <- 2 * max(out.h1$penloglik - out.h0$loglik)
    # crit <- regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, nbtsp = 199 ,parallel = TRUE,cl=cl)$crit
    # lr.crit[k,] <- crit
  }
  crit.h <- regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, nbtsp = 199, an = 10 * an ,parallel = TRUE,cl=cl)$crit
  crit.m <- regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, nbtsp = 199, an = an ,parallel = TRUE,cl=cl)$crit
  crit.l <- regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, nbtsp = 199, an = 0.1 * an ,parallel = TRUE,cl=cl)$crit
  for ( k in 1:nrep){
    lr.crit.h[k,] <- crit.h
    lr.crit.m[k,] <- crit.m
    lr.crit.l[k,] <- crit.l
  }
  
  for ( k in 1:nrep){
    lr.size.l[k,] <- 1 * (lr.estimate[k,] > lr.crit.l[k,2])
    lr.size.m[k,] <- 1 * (lr.estimate[k,] > lr.crit.m[k,2])
    lr.size.h[k,] <- 1 * (lr.estimate[k,] > lr.crit.h[k,2])
  }

  return(list(nominal.size.l = apply(lr.size.l,2,mean), nominal.size.m = apply(lr.size.m,2,mean),  nominal.size.h = apply(lr.size.h,2,mean) ))
}



#GeneratePhiDataPairs
count <- 0
nrep <- 1000
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
      cl <- makeCluster(7)
      
      t <- Sys.time()
      phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
                 beta = beta, N = N, T = T, M = M, p = p, q = q)
      
      phi.data.pair <- GenerateSample(phi,nrep)
      count <- count + 1
      Data = phi.data.pair$Data
      phi = phi.data.pair$phi
      
      an <- anFormula(phi,M,N,T)  #The an function according the the empirical regression
      result <- getEstimateDiffAn(Data,nrep,an,cl)
      
      
      result.l[r, count] <- result$nominal.size.l
      result.h[r, count] <- result$nominal.size.h
      result.m[r, count] <- result$nominal.size.m
      
        
      print(Sys.time() - t)
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


write.csv(result.h,file="results/sizeTest/sizeTestM2BootH.csv")
write.csv(result.m,file="results/sizeTest/sizeTestM2BootM.csv")
write.csv(result.l,file="results/sizeTest/sizeTestM2BootL.csv")
# registerDoParallel(detectCores())
# t <- Sys.time()
# foreach(i=1:2, .combine = rbind)%dopar%{
#   library(NormalRegPanelMixture)
#   phi <- phi.data[[i]]$phi
#   Data = phi.data.pair$Data
#   an <- anFormula(phi,M,phi$N,phi$T)
#   result <- getEstimate(Data,nrep,an)
#   cbind(result$nominal.size,phi$N,phi$T,phi$alpha,phi$mu,phi$sigma)
# }
# print(Sys.time() - t)
