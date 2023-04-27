library(NormalRegPanelMixture)
library(foreach)
set.seed(123456)
options(warn = -1)
# nrep <- 500
# cl <- makeCluster(64)
nrep <- 500
cl <- makeCluster(15)
M <- 2 #Number of Type
p <- 0 #Number of Z
q <- 1 #Number of X

Nset <- c(100,500)
Tset <- c(2,5)
alphaset <- list(c(0.2,0.8))
muset <- list(c(-0.1,0.1),c(-0.5,0.5))
sigmaset <- list(c(0.3,0.1),c(0.5,0.5))
betaset <- list(c(1,1),c(-1,1))
anset <- c(1e-3,1e-2,1e-1)


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

GetMisclTerm <- function(phi) {
  m <- phi$M
  
  if (m == 2)
  {
    omega.12  <- omega.12(phi)
    return (log(omega.12 /(1-omega.12)))
  }
  
  if (m == 3) # ln(omega_12 omega_23 / (0.5-omega_12)(0.5-omega_23))
  {
    omega.123 <- omega.123(phi)
    omega.12 <- omega.123[1]
    omega.23 <- omega.123[2]
    return (log(omega.12 * omega.23 / ((1-omega.12)*(1-omega.23))))
  }
  omega.1234 <- omega.1234(phi)
  omega.12 <- omega.1234[1]
  omega.23 <- omega.1234[2]
  omega.34 <- omega.1234[3]
  # (m == 4) # ln(omega_12 omega_23 omega_34 / (0.5-omega_12)(0.5-omega_23)(0.5-omega_34))
  return (log(omega.12 * omega.23 * omega.34 /
                ((1-omega.12)*(1-omega.23)*(1-omega.34))))
  
}


getEstimate <- function(Data,phi,nrep,an,m,parlist,cl){
  registerDoParallel(cl)
  # results <- foreach (k = 1:nrep)%dopar% {
  #   library(NormalRegPanelMixture)
  #   data <- Data[,k]
  #   out.h0 <- NormalRegPanelMixture::regpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=m,vcov.method = "none")
  #   out.h1 <- NormalRegPanelMixture::regpanelmixMaxPhi(y=data$Y,x=data$X, parlist=out.h0$parlist,an=an,parallel = FALSE)
  #   
  #   crit <- try(NormalRegPanelMixture::regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, parallel = FALSE, nrep=1000)$crit)
  #   if (class(crit) == "try-error"){
  #     crit <- NormalRegPanelMixture::regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, an = an, z = data$Z, parallel = FALSE)$crit
  #   }
  #   c(2 * max(out.h1$penloglik - out.h0$loglik),crit)
  # }
  # lr.estimate <- t(t(sapply(results, function(x) x[1])))
  # 
  # lr.crit <- t(sapply(results, function(x) x[2:length(x)]))
  # return(list(est = lr.estimate , crit = lr.crit,nominal.size = apply(lr.size,2,mean)))
  lr.estimate <- rep(0,nrep)
  lr.crit <- matrix(0,nrow=nrep,ncol=3)
  for (k in 1:nrep){
    data <- Data[,k]
    # print(k)
    out.h0 <- NormalRegPanelMixture::regpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=m,vcov.method = "none")
    out.h1 <- NormalRegPanelMixture::regpanelmixMaxPhi(y=data$Y,x=data$X, parlist=out.h0$parlist,an=an,parallel = TRUE, cl=cl)
    crit <- try(NormalRegPanelMixture::regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, parallel = TRUE, nrep=1000, cl=cl)$crit)
    if (class(crit) == "try-error"){
        crit <- NormalRegPanelMixture::regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, an = an, z = data$Z, parallel = TRUE, cl=cl)$crit
    }
    lr.estimate[k] <- 2 * max(out.h1$penloglik - out.h0$loglik)
    lr.crit[k,] <- crit
  }
  lr.size <- 1 * (lr.estimate > lr.crit[,2])
  return(list(est = lr.estimate , crit = lr.crit,nominal.size = mean(lr.size)))
}

#GeneratePhiDataPairs


phi.data <- list()
nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset) * length(sigmaset) * length(betaset)
regression.data <- matrix(0,nr=(nset*length(anset)),nc=5)

count <- 0
for (N in Nset){
  for (T in Tset){
    for (mu in muset){
      for (alpha in alphaset){
        for (sigma in sigmaset){
          for (beta in betaset){
            for (an in anset){
              t <- Sys.time()
              phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
                         beta = beta, N = N, T = T, M = M, p = p, q = q)
              
              count <- count + 1
              phi.data.pair <- GenerateSample(phi,nrep)
              Data = phi.data.pair$Data
              phi = phi.data.pair$phi
              
              # phi.data[[count]] <- phi.data.pair
              parlist = list(alpha = alpha, mubeta = mu, sigma=sigma, beta=beta, gam=NULL)
              result <- getEstimate(Data,phi,nrep,an,m=M,parlist=parlist,cl=cl)
              
              omega <- GetMisclTerm(phi)
              print(result$nominal.size)
              regression.data[count, ] <-
                cbind(result$nominal.size,an, phi$N,phi$T,omega)
              print(Sys.time() - t)
              print(paste(count,"/",(nset*length(anset))) )
            }
          }
        }
      }
    }
  }
}


colnames(regression.data) <- c("nom.size", "an" ,"N","T", "omega")

write.csv(regression.data,file="penaltyTestM2Regressor.csv",row.names=FALSE)
# write.csv(regression.data,file="/home/haoyu/results/penaltyTestM2Regressor.csv",row.names=FALSE)
