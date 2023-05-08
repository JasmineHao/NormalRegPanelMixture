library(NormalRegPanelMixture)
library(foreach)

set.seed(123456)

#Generate Data
M <- 3 #Number of Type
p <- 0 #Number of Z
q <- 0 #Number of X


Nset <- c(100,500)
Tset <- c(2,5)
alphaset <- list(c(1/3,1/3,1/3),c(0.25,0.5,0.25))
muset <- list(c(-1,0,1),c(-1.5,0,1.5),c(-1,0,2),c(-0.5,0,1.5))
sigmaset <- list(c(1,1,1),c(0.6,1.2,0.6),c(0.6,0.6,1.2))
#GeneratePhiDataPairs
count <- 0
nrep <- 2000
nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset) * length(sigmaset)

cl <- makeCluster(detectCores()-1)


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


getEstimate <- function(Data,nrep,cl,M, parlist){
  registerDoParallel(cl)
  results <- foreach (k = 1:nrep)%dopar% {
    library(NormalRegPanelMixture)
    data <- Data[,k]
    T <- dim(data$Y)[1]
    N <- dim(data$Y)[2]
    out.h0 <- NormalRegPanelMixture::normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none")
    an <- NormalRegPanelMixture::anFormula(out.h0$parlist,M,N,T) 
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



phi.data <- list()
power.data <- matrix(0,nr=(nset),nc=6)

for (mu in muset){
  for (sigma in sigmaset){
    for (alpha in alphaset){
      for (N in Nset){
        for (T in Tset){
          count <- count + 1
          
          t <- Sys.time()
          phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
                     beta = beta, N = N, T = T, M = M, p = p, q = q)
          
          phi.data.pair <- GenerateSample(phi,nrep)
          Data = phi.data.pair$Data
          phi = phi.data.pair$phi
          
          
          print(paste('N is' ,N))
          print(paste('T is' ,T))
          print(paste('mu is' , paste(mu)))
          print(paste('alpha is', paste(alpha)))
          
          parlist = list(alpha = alpha, mubeta = mu, sigma = sigma, gam = NULL)
          result <- getEstimate(Data,nrep,cl,M-1, parlist)
          
          
          power.data[count, 1:3] <- 
            cbind(100*(result$nominal.size),phi$T,phi$N)
          
          power.data[count,4] <- paste(mu,collapse=",")
          power.data[count,5] <- paste(alpha,collapse=",")
          power.data[count,6] <- paste(sigma,collapse=",")
          print(paste("Power is",100*(result$nominal.size)))
          print(Sys.time() - t)
        }
      }
    }
  }
}


write.csv(power.data,file=paste("results/powerTest/powerTestM2SimPLR.csv", sep="_"))

