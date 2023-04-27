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

# ----------------------------------------------------------------------------
# **********************************************
# case: q = 1, m = 2
# 
# **********************************************
# ----------------------------------------------------------------------------
M <- 2
N <- 200
T <- 5
gamma <- NULL
mu <- c(-1,1)
sigma <- c(0.5, 0.5)
beta <- c(0,0)
alpha <- c(0.5,0.5)
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
  out.h1 <- NormalRegPanelMixture::regpanelmixMaxPhi(y=data$Y,x=data$X, parlist=out.h0$parlist,an=(an),parallel = FALSE)
}

M_LN_SQRT_2PI = 0.918938533204672741780329736406 #followed the definition of Rmath.R


x <- data$X
y <- data$Y 
theta <- out.h0$coefficients
nt <- N * T
# compare the M_0 = 1 direct optimization
yrep <- cbind(y,y)
x1 <- cbind(1,x)
obj_0 <- function(theta){
  a  <- c(theta[1], 1- theta[1])
  b  <- matrix(theta[3:6],ncol=2,nrow=2)
  c  <- theta[7:8]
  ll <- -nt * M_LN_SQRT_2PI
  
  eps <- t(t(yrep - x1 %*% b)/ c)
  r_t <- 0.5 * eps^2
  r <- matrix(0,nrow=N,ncol=M)
  for (nn in 1:N){
    r[nn,] = colSums(r_t[((nn-1) * T + 1):(nn*T),])
  }
  r   <- t(t(r) + T * log(c))
  
  minr <- apply(r,1,min)
  l_j <- t( a * t(exp( minr-r)) )
  l_j_sum <- apply(l_j,1,sum)
  ll <- ll + sum( log(l_j_sum) - minr)
  -ll
}


density_check <- function(theta){
  a  <- c(theta[1], 1- theta[1])
  b  <- matrix(theta[3:6],ncol=2,nrow=2)
  c  <- theta[7:8]
  
  eps <- t(t(yrep - x1 %*% b)/ c)
  
  
  eps1 <- apply(dnorm(matrix(eps[,1],nrow=T)),2,prod)
  eps2 <- apply(dnorm(matrix(eps[,2],nrow=T)),2,prod)
  
  ll <- sum(log(a[1] * eps1 + a[2] * eps2))
}

tmp <- cppRegPanelmixPMLE(b0, y, x, ztilde, mu0, sigma0, M, p, t, an, maxit.short,
                   ninits.short, epsilon.short)

print(out.h0$loglik)

obj_0(theta)

optim_2 <- slsqp(theta,obj_0)
print(optim_2$par)
print(optim_2$value)
# print(optim(theta <- out.h0$coefficients, obj_0, hessian=TRUE)) cannot be solved

print(out.h0$coefficients)


# ----------------------------------------------------------------------------
# 检查cpp file写的对不对
# ----------------------------------------------------------------------------
x <- data$X
y <- data$Y 
m <- M
z = NULL; ninits = 25; epsilon = 1e-08; maxit = 2000;
epsilon.short = 1e-02; maxit.short = 500; binit = NULL; in.coefficient = NULL
t <- dim(y)[1] # Number of year
n <- dim(y)[2] # Number of firms
nt <- n * t
p <- 0
y <- as.vector(y)
gam <- NULL
x1  <- cbind(1, x)
q1   <- ncol(x1)

# ninits.short <- ninits * 10 * (1 + p) * m

ninits.short <- 1
tmp <- regpanelmixPMLEinit(y = y, x = x, z = z, ninits = ninits.short, m = M)

sd0 <- sd(y) * sqrt((n - 1) / n)
an <- 1 / n # penalty term for variance
sigma0 <- rep(sd0, M)
mu0 <- double(M + 1) # dummy
ztilde <- matrix(0) # dummy
b0 <- as.matrix(rbind(tmp$alpha, tmp$mu, tmp$sigma, tmp$gam))
out.short <- cppRegPanelmixPMLE(b0, y, x, ztilde, mu0, sigma0, M, p, t, an, maxit.short,
                                ninits.short, epsilon.short)

components <- order(out.short$penloglikset, decreasing = TRUE)[1:ninits]
b1 <- b0[ ,components] # b0 has been updated
if ((! is.null(in.coefficient)) & (dim(b1)[1]==length(in.coefficient)) ){
  b1[,components] <- in.coefficient
}

out <- cppRegPanelmixPMLE(b1, y, x, ztilde, mu0, sigma0, m, p, t, an, maxit, ninits, epsilon)
index     <- which.max(out$penloglikset)
alpha <- b1[1:m,index] # b0 has been updated

mubeta <- matrix(b1[(1+m):((q1+1)*m),index],nrow=q1,ncol=m)
# mubeta <- matrix(b1[3:10,index],nrow=q1,ncol=m)
sigma <- b1[(1+(q1+1)*m):((q1+2)*m),index]
if (!is.null(z)) {
  gam     <- b1[((q1+2)*m+1):((q1+2)*m+p),index]
}
penloglik <- out$penloglikset[index]
loglik    <- out$loglikset[index]

# ----------------------------------------------------------------------------
# **********************************************
# case: q = 1, m = 1
# 不知道为什么m = 1 的情况cpp算起来不太对
# 
# **********************************************
# ----------------------------------------------------------------------------

N <- 200
T <- 5
M <- 1
gamma <- NULL
alpha <- c(1)
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
  out.h1 <- NormalRegPanelMixture::regpanelmixMaxPhi(y=data$Y,x=data$X, parlist=out.h0$parlist,an=(an),parallel = FALSE)
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

x <- data$X
y <- data$Y 
theta <- out.h0$coefficients[2:4]
# compare the M_0 = 1 direct optimization
obj_0 <- function(theta){
  ll <- -nt * M_LN_SQRT_2PI
  eps <- (as.vector(y) - x %*% theta[2] - theta[1])/ theta[3]
  r_t <- 0.5 * eps^2
  r_t <- matrix(r_t,nrow=T,ncol=N)
  r   <- colSums(r_t) + T * log(theta[3])
  l_j <- exp(-r)
  - (ll + sum(log(l_j)))
}

print(optim(theta <- out.h0$coefficients[2:4], obj_0, hessian=TRUE))
print(out.h0$coefficients[2:4])
print(out.h0$loglik)

# ----------------------------------------------------------------------------
# 检查cpp file写的对不对
# ----------------------------------------------------------------------------
x <- data$X
y <- data$Y 
m <- M
z = NULL; ninits = 25; epsilon = 1e-08; maxit = 2000;
epsilon.short = 1e-02; maxit.short = 500; binit = NULL; in.coefficient = NULL
t <- dim(y)[1] # Number of year
n <- dim(y)[2] # Number of firms
nt <- n * t
p <- 0
y <- as.vector(y)
gam <- NULL
x1  <- cbind(1, x)
q1   <- ncol(x1)

ninits.short <- ninits * 10 * (1 + p) * m

# ninits.short <- 1
tmp <- regpanelmixPMLEinit(y = y, x = x, z = z, ninits = ninits.short, m = m)

sd0 <- sd(y) * sqrt((n - 1) / n)
an <- 1 / n # penalty term for variance
sigma0 <- rep(sd0, m)
mu0 <- double(m + 1) # dummy
ztilde <- matrix(0) # dummy
b0 <- as.matrix(rbind(tmp$alpha, tmp$mu, tmp$sigma, tmp$gam))
out.short <- cppRegPanelmixPMLE(b0, y, x, ztilde, mu0, sigma0, m, p, t, an, maxit.short,
                                ninits.short, epsilon.short)

components <- order(out.short$penloglikset, decreasing = TRUE)[1:ninits]
b1 <- b0[ ,components] # b0 has been updated
if ((! is.null(in.coefficient)) & (dim(b1)[1]==length(in.coefficient)) ){
  b1[,components] <- in.coefficient
}

out <- cppRegPanelmixPMLE(b1, y, x, ztilde, mu0, sigma0, m, p, t, an, maxit, ninits, epsilon)
index     <- which.max(out$penloglikset)
alpha <- b1[1:m,index] # b0 has been updated

mubeta <- matrix(b1[(1+m):((q1+1)*m),index],nrow=q1,ncol=m)
# mubeta <- matrix(b1[3:10,index],nrow=q1,ncol=m)
sigma <- b1[(1+(q1+1)*m):((q1+2)*m),index]
if (!is.null(z)) {
  gam     <- b1[((q1+2)*m+1):((q1+2)*m+p),index]
}
penloglik <- out$penloglikset[index]
loglik    <- out$loglikset[index]


# ----------------------------------------------------------------------------
# case: q = 0, m = 1
# ----------------------------------------------------------------------------
q <- 0
p <- 0
M <- 1
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
  out.h1 <- NormalRegPanelMixture::normalpanelmixMaxPhi(y=data$Y, parlist=out.h0$parlist,an=(an),parallel = FALSE)
}

# ----------------------------------------------------------------------------
# 仿写c++ loglikelihood
x <- data$X
y <- data$Y 
theta <- out.h0$coefficients[2:3]
obj_0 <- function(theta){
  ll <- -nt * M_LN_SQRT_2PI
  eps <- (as.vector(y) - theta[1])/ theta[2]
  r_t <- 0.5 * eps^2
  r_t <- matrix(r_t,nrow=T,ncol=N)
  r   <- colSums(r_t) + T * log(theta[2])
  l_j <- exp(-r)
  - (ll + sum(log(l_j)))
}
print(optim(theta <- out.h0$coefficients[2:3], obj_0, hessian=TRUE))
print(out.h0$coefficients[2:3])
print(out.h0$loglik)
# ----------------------------------------------------------------------------
# 检查cpp file写的对不对
# ----------------------------------------------------------------------------
x <- data$X
y <- data$Y 
m <- M
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

