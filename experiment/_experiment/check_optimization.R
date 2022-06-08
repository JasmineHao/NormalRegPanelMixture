library(nloptr)
library(NormalRegPanelMixture)
library(doParallel)

#Generate Data
M <- 1 #Number of Type
p <- 0 #Number of Z
q <- 0 #Number of X

nrep <- 1

cl <- makeCluster(2)


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

count <- 0
phi.data <- list()
nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset) * length(betaset)

NTset <- expand.grid(Nset,Tset)
Parset <- expand.grid(muset,alphaset,betaset,sigmaset)
nNT <- dim(NTset)[1]
nPar <- dim(Parset)[1]



N <- 200
T <- 5
gamma <- NULL
mu <- c(0)
sigma <- c(1)
beta <- c(1)
q <- 1
p <- 0
t <- Sys.time()
phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma, beta = beta, N = N, T = T, M = M, p = p, q = q, X=NULL)
          
phi.data.pair <- GenerateSample(phi,nrep)
Data = phi.data.pair$Data
phi = phi.data.pair$phi
an <- anFormula(phi,M,N,T,q=1) #The an function according the the empirical regression
print(an)

for(k in 1:nrep){
  data <- Data[,k]
  out.h0 <- NormalRegPanelMixture::regpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none")
  out.h1 <- NormalRegPanelMixture::regpanelmixMaxPhi(y=data$Y,x=data$X, parlist=out.h0$parlist,an=(an),update.alpha = 1,parallel = FALSE)
}

M_LN_SQRT_2PI = 0.918938533204672741780329736406 #followed the definition of Rmath.R

loglikelihood <- function(y,x=NULL,parlist=NULL){
  nt <- prod(dim(y))
  mubeta <- parlist$mubeta
  sigma  <- parlist$sigma
  alpha  <- parlist$alpha
  alp_sig <- alpha/sigma
  r <- as.vector(y) - x %*% mubeta[2:nrow(mubeta),] - mubeta[1,]
  r <- r / sigma
  r   <- 0.5 * r * r
  
  f <- exp(-r) / sqrt(2*pi) / sigma
  
  return(sum(log(f)))
}





obj <- function(theta){
  sum ( 0.5*(as.vector(y) - x %*% theta[2] - theta[1] )^2/theta[3] + 0.5* log(theta[3]) )
  
}


x0 <- out.h0$coefficients[2:4]
slsqp(x0,obj)

# compare the M_0 = 1 direct optimization
res <- optim(theta <- out.h0$coefficients[2:4], obj, hessian=TRUE)


# ----------------------------------------------------------------------------
q <- 0
p <- 0
t <- Sys.time()
phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma, N = N, T = T, M = M, p = p, q = q, X=NULL)

phi.data.pair <- GenerateSample(phi,nrep)
Data = phi.data.pair$Data
phi = phi.data.pair$phi
an <- anFormula(phi,M,N,T,q=1) #The an function according the the empirical regression
print(an)

for(k in 1:nrep){
  data <- Data[,k]
  out.h0 <- NormalRegPanelMixture::normalpanelmixPMLE(y=data$Y, z = data$Z,m=M,vcov.method = "none")
  out.h1 <- NormalRegPanelMixture::normalpanelmixMaxPhi(y=data$Y, parlist=out.h0$parlist,an=(an),update.alpha = 1,parallel = FALSE)
}

# ----------------------------------------------------------------------------
# 仿写c++ loglikelihood
obj_0 <- function(theta){
  ll <- -nt * M_LN_SQRT_2PI
  eps <- (as.vector(y) - theta[1])/ theta[2]
  r_t <- 0.5 * eps^2
  r_t <- matrix(r_t,nrow=T,ncol=N)
  r   <- colSums(r_t) + T * log(theta[2])
  l_j <- exp(-r)
  ll + sum(log(l_j))
}
  
# ----------------------------------------------------------------------------
obj <- function(theta){
  sum ( 0.5*(as.vector(y) - theta[1] )^2/theta[2] + 0.5* log(theta[2]) )
}

x <- data$X
y <- data$Y 
x0 <- out.h0$coefficients[2:3]
optim(theta <- out.h0$coefficients[2:3], obj, hessian=TRUE)

# ----------------------------------------------------------------------------
# 检查cpp file写的对不对
# ----------------------------------------------------------------------------
x <- data$X
y <- data$Y 
z = NULL; ninits = 25; epsilon = 1e-08; maxit = 2000;
epsilon.short = 1e-02; maxit.short = 500; binit = NULL; in.coefficient = NULL
t <- dim(y)[1] # Number of year
n <- dim(y)[2] # Number of firms
nt <- n * t
p <- 0
y <- as.vector(y)
gam <- NULL
ninits.short <- ninits * 10 * (1 + p) * m

ninits.short <- 1
tmp <- normalpanelmixPMLEinit(y = y, z = z, ninits = ninits.short, m = m)

sd0 <- sd(y) * sqrt((n - 1) / n)
an <- 1 / n # penalty term for variance
sigma0 <- rep(sd0, m)
mu0 <- double(m + 1) # dummy
ztilde <- matrix(0) # dummy
b0 <- as.matrix(rbind(tmp$alpha, tmp$mu, tmp$sigma, tmp$gam))
out.short <- cppnormalpanelmixPMLE(
  b0, y, ztilde, mu0, sigma0, m, p, t, an, maxit.short,
  ninits.short, epsilon.short
)

