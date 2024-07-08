library(NormalRegPanelMixture)
library(foreach)

#Generate Data
M <- 2 #Number of Type
p <- 0 #Number of Z
q <- 3 #Number of X
p.0 <- round(( p - 1 ) / 2)
q.0 <- round(( q - 1 ) / 2) 
nrep <- 100
cl <- makeCluster(16)

set.seed(123456)
Nset <- c(200)
Tset <- c(3)


alphaset <- list(c(0.5,0.5))

rho <- c(0.5,0.5)
beta.r1 <- c(-1,1)
beta.r2 <- c(1,1)

mu.r <- c(-1,1)

muset <- list(mu.r * (1-rho))
mu0set <- list(mu.r / sqrt(1- rho**2))


betaset <- list( t(matrix(rbind(rho, beta.r1, -rho * beta.r1),nrow = q)), t(matrix(rbind(rho, beta.r2, -rho * beta.r2),nrow = q)) )
beta0set <- list( (beta.r1 / sqrt(1 - rho**2)),   (beta.r2 / sqrt(1 - rho**2)))


sigma <- c(0.75,0.75)
sigma0 <- sqrt(sigma**2 / (1 - rho**2))


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
  beta0 = phi$beta0
  
  Data <- replicate(nrep,generateDataAR1(
    alpha,mu,sigma,gamma,beta,mu0,sigma0,gamma0,
    beta0, N, T ,M,p,q,p.0,q.0))
  return(list(phi=phi,Data=Data))
}



getEstimateDiffAn <- function(Data,nrep,an,an_0,cl,M){
  registerDoParallel(cl)
  results <- foreach (k = 1:nrep)%dopar% {
    library(NormalRegPanelMixture)
    data <- Data[,k]
    data.0 <- list(Y = data$Y0, X = data$X0, Z = data$Z0)
    out.h0 <- NormalRegPanelMixture::regpanelmixPMLE(
      y=data$Y,x=data$X, z = data$Z,m=M-1,vcov.method = "none", data.0 = data.0)
    
    out.h1.m <- NormalRegPanelMixture::regpanelmixMaxPhi(y=data$Y,x=data$X,
                                                         parlist=out.h0$parlist,an=(an), 
                                                         an_0 = (an_0),parallel = FALSE, data.0 = data.0)
    
    
    
    
    crit.m <- NormalRegPanelMixture::regpanelmixCritBootAR1(y = data$Y, x = data$X, 
                                                            parlist = out.h0$parlist,  z = data$Z, parallel = FALSE,  data.0 = data.0, an=(an),  ninits = 1, an_0 = ( an_0), nbtsp = 199)$crit
    
    c( 2 * max(out.h1.m$penloglik - out.h0$loglik), crit.m)
    
  }
  
  
  
  lr.estimate.m <- t(t(sapply(results, function(x) x[1])))
  
  lr.crit.m <- t(sapply(results, function(x) x[2:4]))
  
  lr.size.m <- matrix(0.0,nr=nrep,ncol=1) #Nomimal size
  
  
  
  for ( k in 1:nrep){
    lr.size.m[k,] <- 1 * (lr.estimate.m[k,] > lr.crit.m[k,2])
  }
  
  return(list(nominal.size.m = apply(lr.size.m,2,mean)
  ))
}




#GeneratePhiDataPairs
count <- 1
r <- 1
phi.data <- list()
nset <- length(Nset) * length(Tset) * length(muset) * 
  length(alphaset) * length(mu0set) * length(betaset)


NTset <- expand.grid(Nset,Tset)
Parset <- expand.grid(muset,alphaset,betaset)
Par0set <- expand.grid(mu0set,alphaset,beta0set)
nNT <- dim(NTset)[1]
nPar <- dim(Parset)[1]

result.m <- matrix(0,nr=(nNT),nc=nPar)
rownames(result.m) <- apply(NTset,1,paste,collapse = ",")
colnames(result.m) <- apply(Parset,1,paste,collapse = ",")

for (r in 1:nNT){
  N <-  NTset[r,1]
  T <-  NTset[r,2]
  
  for (count in 1:nPar){
    mu <- Parset[count,1][[1]]
    alpha <- Parset[count,2][[1]]
    beta <- Parset[count,3][[1]]
    
    mu0  <- Par0set[count,1][[1]]
    beta0 <- Par0set[count,3][[1]]
    
    t <- Sys.time()
    phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = NULL, beta = beta, beta0 = beta0, 
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
    an <- anFormula(phi,M,N,T_an)  #The an function according the the empirical regression
    an_0 <- anFormula.t0(phi.0, M, N, q = q.0)
    
    # an <- 0.03
    print(N)
    print(T)
    print(mu)
    print(alpha)
    result <- getEstimateDiffAn(Data,nrep,an, an_0, cl,M)
    
    
    result.m[r, count] <- result$nominal.size.m
    
    print(result$nominal.size.m)
    
    print(Sys.time() - t)
    
  }
}


result.m <- result.m * 100

write.csv(result.m, file="results/powerTestM2SimRegressor_AR1_M.csv")
