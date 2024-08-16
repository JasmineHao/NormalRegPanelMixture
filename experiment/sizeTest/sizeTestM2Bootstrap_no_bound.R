# Generate table 3 and 4
library(NormalRegPanelMixture)
library(foreach)

M <- 2 #Number of Type
p <- 0 #Number of Z
q <- 0 #Number of X
nrep <- 100

cl <- makeCluster(12)


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



getEstimateDiffAn <- function(Data,nrep,an,cl,M, parlist){
  registerDoParallel(cl)
  
  # Remove this: 
  #data <- Data[,1]
  #lr.crit.m <- NormalRegPanelMixture::regpanelmixCritBoot(y=data$Y, x=data$X, parlist=parlist, an=(an), z = data$Z, parallel = TRUE, nbtsp = 499, ninits = 10)$crit
  #lr.crit.l <- NormalRegPanelMixture::regpanelmixCritBoot(y=data$Y, x=data$X, parlist=parlist, an=(0.01 * an), z = data$Z, parallel = TRUE, nbtsp = 499, ninits = 10)$crit
  #lr.crit.h <- NormalRegPanelMixture::regpanelmixCritBoot(y=data$Y, x=data$X, parlist=parlist, an=(100 * an), z = data$Z, parallel = TRUE, nbtsp = 499, ninits = 10)$crit
  
  
  results <- foreach (k = 1:nrep)%dopar% {
    library(NormalRegPanelMixture)
    data <- Data[,k]
    out.h0 <- NormalRegPanelMixture::normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none")
    out.h1.m <- NormalRegPanelMixture::normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an=(an),parallel = FALSE,eps=0)
    out.h1.l <- NormalRegPanelMixture::normalpanelmixMaxPhi(y = data$Y, parlist = out.h0$parlist, an = ( 0.01 * an), parallel = FALSE,eps=0)
    out.h1.h <- NormalRegPanelMixture::normalpanelmixMaxPhi(y = data$Y, parlist = out.h0$parlist, an = (100 * an), parallel = FALSE,eps=0)
    
    crit.m <- NormalRegPanelMixture::regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, an=(an), z = data$Z, parallel = FALSE, nbtsp = 199, ninits = 10)$crit
    crit.l <- NormalRegPanelMixture::regpanelmixCritBoot(y = data$Y, x = data$X, parlist = out.h0$parlist, an = (0.01 * an), z = data$Z, parallel = FALSE, nbtsp = 199, ninits = 10)$crit
    crit.h <- NormalRegPanelMixture::regpanelmixCritBoot(y = data$Y, x = data$X, parlist = out.h0$parlist, an = (100 * an), z = data$Z, parallel = FALSE, nbtsp = 199, ninits = 10)$crit
    
    c(2 * max(out.h1.l$penloglik - out.h0$loglik), 2 * max(out.h1.m$penloglik - out.h0$loglik), 2 * max(out.h1.h$penloglik - out.h0$loglik) , crit.l, crit.m, crit.h)
    
  }
  
  
  lr.estimate.l <- t(t(sapply(results, function(x) x[1])))
  lr.estimate.m <- t(t(sapply(results, function(x) x[2])))
  lr.estimate.h <- t(t(sapply(results, function(x) x[3])))
  
  lr.crit.l <- t(sapply(results, function(x) x[4:6]))
  lr.crit.m <- t(sapply(results, function(x) x[7:9]))
  lr.crit.h <- t(sapply(results, function(x) x[10:12]))
  
  
  lr.size.l <- matrix(0.0,nr=nrep,ncol=1) #Nomimal size
  lr.size.m <- matrix(0.0,nr=nrep,ncol=1) #Nomimal size
  lr.size.h <- matrix(0.0,nr=nrep,ncol=1) #Nomimal size
  
  for ( k in 1:nrep){
    lr.size.l[k,] <- 1 * (lr.estimate.l[k,] > lr.crit.l[k,2])
    lr.size.m[k,] <- 1 * (lr.estimate.m[k,] > lr.crit.m[k,2])
    lr.size.h[k,] <- 1 * (lr.estimate.h[k,] > lr.crit.h[k,2])
  }
  
  return(list(nominal.size.l = apply(lr.size.l,2,mean), nominal.size.m = apply(lr.size.m,2,mean),  nominal.size.h = apply(lr.size.h,2,mean) ))
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


print("Simulation Started!!")

for (r in 1:nNT){
  N <-  NTset[r,1]
  T <-  NTset[r,2]
  
  
  for (count in 1:nPar){
    mu <- Parset[count,1][[1]]
    alpha <- Parset[count,2][[1]]
    sigma <- sigmaset[[1]]
    t <- Sys.time()
    phi = list(alpha = alpha, mu = mu, sigma = sigma, gamma = NULL, beta = NULL, N = N, T = T, M = M, p = p, q = q)
    
    phi.data.pair <- GenerateSample(phi,nrep)
    
    Data = phi.data.pair$Data
    phi = phi.data.pair$phi
    
    an <- anFormula.alt(phi,M,N,T)  #The an function according the the empirical regression
    print(N)
    print(T)
    print(mu)
    print(alpha)
    print(an)
    parlist = list(alpha = alpha, mubeta = mu, sigma=sigma, gam=NULL)
    result <- getEstimateDiffAn(Data,nrep,an,cl,M, parlist)
    print(r)
    print(count)
    print(result$nominal.size.m)
    
    result.l[r, count] <- result$nominal.size.l
    result.m[r, count] <- result$nominal.size.m
    result.h[r, count] <- result$nominal.size.h
    
    print(result.m[r, count])
    print(Sys.time() - t)
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

write.csv(rbind(result.h,result.m,result.l), file="results/sizeTestM2Boot_HML_T2_no_bound.csv")


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
