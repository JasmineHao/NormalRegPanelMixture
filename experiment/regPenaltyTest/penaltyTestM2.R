library(NormalRegPanelMixture)
library(foreach)

nrep <- 500
cl <- makeCluster(64)
set.seed(123456)
#Generate Data
Nset <- c(100,500)
Tset <- c(2,5,10)

alphaset <- list(c(0.5,0.5),c(0.2,0.8))
muset <- list(c(-1,1),c(-0.5,0.5), c(-0.8,0.8))
sigmaset <- list(c(1, 1), c(1.5, 0.75), c(0.8,1.2))
anset <- c(0.01, 0.05,0.1,0.2, 0.3, 0.4)


#The parameters that are fixed
M <- 2 #Number of Type
p <- 0 #Number of Z
q <- 0 #Number of X
gamma <- matrix(0)
beta <- matrix(0)

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
  results <- foreach (k = 1:nrep)%dopar% {
    library(NormalRegPanelMixture)
    data <- Data[,k]
    out.h0 <- NormalRegPanelMixture::normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=m,vcov.method = "none")
    out.h1 <- NormalRegPanelMixture::normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an=an,parallel = FALSE)
    crit <- try(NormalRegPanelMixture::regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, parallel = FALSE, nrep=1000)$crit)
    if (class(crit) == "try-error"){
      crit <- NormalRegPanelMixture::regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, an = an, z = data$Z, parallel = FALSE)$crit
    }
    c(2 * max(out.h1$penloglik - out.h0$loglik),crit)
  }
  lr.estimate <- t(t(sapply(results, function(x) x[1])))
  lr.crit <- t(sapply(results, function(x) x[2:length(x)]))
  lr.size <- 1 * (lr.estimate > lr.crit[,2])
  return(list(est = lr.estimate , crit = lr.crit,nominal.size = apply(lr.size,2,mean)))
}

#GeneratePhiDataPairs
count <- 0

phi.data <- list()
nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset) * length(sigmaset)
regression.data <- matrix(0,nr=(nset*length(anset)),nc=5)


# ====== BEGIN EXPERIMENT ======
## 1. Initialization
# Case when m = 2
for (N in Nset){
  for (T in Tset){
    for (mu in muset){
      for (alpha in alphaset){
        for (sigma in sigmaset){
          for (an in anset){
            t <- Sys.time()
            phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
                       beta = beta, N = N, T = T, M = M, p = p, q = q)

            count <- count + 1
            phi.data.pair <- GenerateSample(phi,nrep)
            Data = phi.data.pair$Data
            phi = phi.data.pair$phi

            # phi.data[[count]] <- phi.data.pair
            parlist = list(alpha = alpha, mubeta = mu, sigma=sigma, gam=NULL)
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


colnames(regression.data) <- c("nom.size", "an" ,"N","T", "omega")
# Begin estimation

stargazer(regression.data)
#The colnames are nom.size, an, T, N
write.csv(regression.data,file="/home/haoyu/results/penaltyTestM2.csv",row.names=FALSE)


