library(NormalRegPanelMixture)
library(foreach)

#Generate Data
M <- 2 #Number of Type
p <- 0 #Number of Z
q <- 1 #Number of X
p.0 <- 0
q.0 <- 0 
nrep <- 500
cl <- makeCluster(10)

set.seed(123456)
Nset <- c(200,400)
Tset <- c(3,5)

alphaset <- list(c(0.5,0.5),c(0.2,0.8))

muset <- list(c(-1,1),c(0.5,0.5))
mu0set <- list(c(-1,1))
betaset <- list(c(0.5,0.5))
sigmaset <- list(c(0.8,1.2))
sigma0set <- list(c(0.8,1.2))

anFormula.alt <- function(parlist, m, n, t){
  omega <- omega.12(parlist)
  omega <- pmin(pmax(omega, 1e-16), 0.5-1e-16)  # an becomes NaN if omega[j]=0 or 1
  omega.term <- log(omega /(1-omega))
  b <-   c(-0.8112790,  -0.2882271,   4.6374028,  -0.1012959,  -0.1973225)
  x <- (  b[1] + b[2]/t + b[3]/n + b[5] * omega.term ) / b[4]   # maxa=1
  an <- 1 / (1 + exp(x))
  an
}


# generateDataAR1(alpha=alpha,mu=mu,sigma=sigma,gamma=gamma,beta=beta,mu0=mu0,sigma0=sigma0,gamma0=gamma0,beta0=beta0,N = n, T = t,M=m,p=p,q=q,p.0=p.0,q.0=q.0,x=x,z=z, x0 = x0, z0=z0))
# tmp <- lapply(seq_len
GenerateSample <- function(phi,nrep){
  p = phi$p
  q = phi$q
  p.0 = phi$p.0
  q.0 = phi$q.0
  
  N = phi$N
  T = phi$T
  M = phi$M
  alpha = phi$alpha
  sigma = phi$sigma
  sigma0 = phi$sigma0
  mu = phi$mu
  mu0 = phi$mu0
  gamma = phi$gamma
  beta = phi$beta
  
  Data <- replicate(nrep,generateDataAR1(alpha,mu,sigma,gamma,beta,mu0,sigma0,gamma0,beta0,N , T ,M,p,q,p.0,q.0))
  return(list(phi=phi,Data=Data))
}



getEstimateDiffAn <- function(Data,nrep,an,an_0,cl,M, parlist){
  registerDoParallel(cl)
  results <- foreach (k = 1:nrep)%dopar% {
    library(NormalRegPanelMixture)
    data <- Data[,k]
    data.0 <- list(Y = data$Y0, X = data$X0, Z = data$Z0)
    out.h0 <- NormalRegPanelMixture::regpanelmixPMLE(y=data$Y,x=data$X, 
                                                     z = data$Z,m=M,vcov.method = "none", data.0 = data.0)
    out.h1.l <- NormalRegPanelMixture::regpanelmixMaxPhi(y=data$Y,x=data$X,
                                                         parlist=out.h0$parlist,an=(0.01 * an), 
                                                         an_0 = (0.01* an_0),parallel = FALSE, data.0 = data.0)
    out.h1.m <- NormalRegPanelMixture::regpanelmixMaxPhi(y=data$Y,x=data$X,
                                                         parlist=out.h0$parlist,an=(an), 
                                                         an_0 = (an_0),parallel = FALSE, data.0 = data.0)
    out.h1.h <- NormalRegPanelMixture::regpanelmixMaxPhi(y=data$Y,x=data$X,
                                                         parlist=out.h0$parlist,
                                                         an=(100 * an), an_0 = (100 * an_0),
                                                         parallel = FALSE, data.0 = data.0)
    
    crit.l <- NormalRegPanelMixture::regpanelmixCritBootAR1(y = data$Y, x = data$X,
                                                            parlist = out.h0$parlist, 
                                                            z = data$Z, parallel = FALSE, 
                                                            data.0 = data.0, an=(0.01 * an), 
                                                            an_0 = (0.01* an_0), nbtsp = 199)$crit
    crit.m <- NormalRegPanelMixture::regpanelmixCritBootAR1(y = data$Y, x = data$X, 
                                                            parlist = out.h0$parlist, 
                                                            z = data$Z, parallel = FALSE, 
                                                            data.0 = data.0, an=(an), 
                                                            an_0 = ( an_0), nbtsp = 199)$crit
    crit.h <- NormalRegPanelMixture::regpanelmixCritBootAR1(y = data$Y, x = data$X, 
                                                            parlist = out.h0$parlist, 
                                                            z = data$Z, parallel = FALSE, 
                                                            data.0 = data.0, an=(100 * an), 
                                                            an_0 = (100 * an_0), nbtsp = 199)$crit
    
    
    c(2 * max(out.h1.l$penloglik - out.h0$loglik), 2 * max(out.h1.m$penloglik - out.h0$loglik), 2 * max(out.h1.h$penloglik - out.h0$loglik), crit.l, crit.m, crit.h)
    
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
  
  print(lr.crit.l)
  
  for ( k in 1:nrep){
    lr.size.l[k,] <- 1 * (lr.estimate.l[k,] > lr.crit.l[k,2])
    lr.size.m[k,] <- 1 * (lr.estimate.m[k,] > lr.crit.m[k,2])
    lr.size.h[k,] <- 1 * (lr.estimate.h[k,] > lr.crit.h[k,2])
  }
  
  return(list(nominal.size.l = apply(lr.size.l,2,mean),
              nominal.size.m = apply(lr.size.m,2,mean), 
              nominal.size.h = apply(lr.size.h,2,mean) ))
}




#GeneratePhiDataPairs
count <- 0

phi.data <- list()
nset <- length(Nset) * length(Tset) * length(muset) * 
  length(alphaset) * length(mu0set) * length(betaset)


NTset <- expand.grid(Nset,Tset)
Parset <- expand.grid(muset,alphaset, betaset, mu0set)
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
  
  
  for (count in 1:nPar){
    mu <- Parset[count,1][[1]]
    alpha <- Parset[count,2][[1]]
    beta <- Parset[count,3][[1]]
    mu0 <- Parset[count,4][[1]]
    
    sigma <- sigmaset[[1]]
    sigma0 <- sigma0set[[1]]
    
    t <- Sys.time()
    phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = NULL, beta = beta, 
               mu0 = mu0, sigma0 = sigma0 , N = N, T = T, M = M, p = p, q = q, p.0 = p.0, q.0 = q.0)
    phi.0 <- list(alpha = alpha, gamma = NULL, beta = NULL, 
                 mu = mu0, sigma = sigma0)
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
    an <- anFormula.alt(phi,M,N,T_an)  #The an function according the the empirical regression
    an_0 <- anFormula.t0(phi.0, M, N, q = q.0)
    
    # an <- 0.03
    print(N)
    print(T)
    print(mu)
    print(alpha)
    print(anFormula(phi,M,N,T_an))
    print(an)
    parlist = list(alpha = alpha, mubeta = mu, sigma=sigma, gam=NULL)
    result <- getEstimateDiffAn(Data,nrep,an, an_0, cl,M, parlist)
    
    
    result.l[r, count] <- result$nominal.size.l
    result.m[r, count] <- result$nominal.size.m
    result.h[r, count] <- result$nominal.size.h
    print(result$nominal.size.m)
    
    print(Sys.time() - t)
    
  }
}




result.h <- result.h * 100
result.l <- result.l * 100
result.m <- result.m * 100

write.csv(rbind(result.h,result.m,result.l), file="results/sizeTestM2SimH_AR1_HML.csv")
