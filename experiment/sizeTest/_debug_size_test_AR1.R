library(NormalRegPanelMixture)
library(foreach)

#Generate Data
M <- 2 #Number of Type
p <- 0 #Number of Z
q <- 3 #Number of X
p.0 <- round(( p - 1 ) / 2)
q.0 <- round(( q - 1 ) / 2) 
nrep <- 500
cl <- makeCluster(2)

set.seed(123456)
Nset <- c(200)
Tset <- c(3,5)

alphaset <- list(c(0.5,0.5),c(0.2,0.8))

rho <- c(0.5,0.5)
beta.r <- c(1,-1)
mu.r <- c(-1,1)

muset <- list(mu.r * (1-rho))
mu0set <- list(mu.r / sqrt(1- rho**2))


betaset <- list( t(matrix(rbind(rho, beta.r, -rho * beta.r),nrow = q)) )
beta0set <- list(beta.r / sqrt(1 - rho**2))


sigma <- c(0.8,1.2)
sigma0 <- sqrt(sigma**2 / (1 - rho**2))

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

N <- 200
T <- 3

mu <- muset[[1]]
alpha <- alphaset[[2]]
beta <- betaset[[1]]
mu0  <- mu0set[[1]]
beta0 <- beta0set[[1]]

phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = NULL, beta = beta, beta0 = beta0, 
           mu0 = mu0, sigma0 = sigma0 , N = N, T = T, M = M, p = p, q = q, p.0 = p.0, q.0 = q.0)
phi.0 <- list(alpha = alpha, gamma = NULL, beta = NULL, 
              mu = mu0, sigma = sigma0)
phi.data.pair <- GenerateSample(phi,nrep)
Data = phi.data.pair$Data
phi = phi.data.pair$phi

an <- anFormula.alt(phi,M,N,T)  #The an function according the the empirical regression
an_0 <- anFormula.t0(phi.0, M, N, q = q.0)

k <- 1
data <- Data[,k]
data.0 <- list(Y = data$Y0, X = data$X0, Z = data$Z0)

# 1. Debug for PMLE
# --------------------------------
# y = data$Y
# x = data$X
# m = M
# z = data$Z
# ninits = 10
# epsilon = 1e-08
# maxit = 2000
# epsilon.short = 1e-02
# maxit.short = 500
# binit = NULL
# in.coefficient=NULL
# 
# t  <- nrow(y)
# n  <- ncol(y)
# nt <- n*t
# 
# y   <- as.vector(y)
# x   <- as.matrix(x)   # n by (q1-1) matrix
# 
# if (nrow(x) != nt) { stop("y and x must have the same number of rows.") }
# x1  <- cbind(1, x)
# q1   <- ncol(x1)
# 
# p       <- 0
# gam   <- NULL
# ninits.short <- ninits*(q1+p)*m
# 
# y0 <- data.0$Y
# x0 <- data.0$X
# z0 <- data.0$Z
# z.init <- cbind(y0,x0,z0) # use to determine the mixture probability
# z.init <- scale(z.init)
# q1.0 <- ncol(z.init) + 1 # the dimension of gamma_0
# 
# npar    <- m-1 + (q1+1)*m + p + q1.0 * (m-1)  # number of parameters
# 
# xz      <- cbind(x, z)
# ls.out  <- lsfit(xz, y)
# sd0     <- sqrt(mean(ls.out$residuals^2))
# tmp <- regpanelmixPMLEinit(y = y, x = x, z = z, ninits = ninits.short,
#                            m = m, model.ar1=TRUE, z.init =z.init)
# 
# 
# h       <- 0    # setting h=0 gives PMLE
# tau     <- 0.5  # setting tau=0.5 gives PMLE
# k <- 0 # setting k=0 gives PMLE
# 
# an <- 1 / n # penalty term for variance
# # an_0 <- 1 / n # penalty term for variance
# an_0 <- 0.3 # use the default from KS 15
# 
# sigma.0 <- rep(sd0, m)
# mu.0 <- double(m + 1) # dummy
# 
# # # Begin of AR1 case
# z.init <- as.matrix(z.init)
# # short EM
# b0 <- rbind(tmp$alpha, tmp$mubeta, tmp$sigma, tmp$gam, tmp$gamma0)
# if (!is.null(binit)) {
#   b0[, 1] <- binit
# }
# #
# if (is.null(z)){
# ztilde <- matrix(0) # dummy
# } else {
#   ztilde <- z
# }
# out.short <- cppRegPanelmixPMLEAR1(
#   b0, y, x, ztilde, cbind(1,z.init), mu.0, sigma.0, m, p, t, an, maxit.short,
#   ninits.short, epsilon.short
# )
# #
# # #
# # # long EM
# components <- order(out.short$penloglikset, decreasing = TRUE)[1:ninits]
# b1 <- b0[, components] # b0 has been updated
# if ((!is.null(in.coefficient)) & (dim(b1)[1] == length(in.coefficient))) {
#   b1[, components] <- in.coefficient
# }
# out <- cppRegPanelmixPMLEAR1(
#   b1, y, x, ztilde, cbind(1,z.init), mu.0, sigma.0, m, p, t, an, maxit,
#   ninits, epsilon
# )

# 
# 
t = Sys.time() 
out.h0 <- regpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none", data.0 = data.0)

# 2. Debug for MaxPhi
# --------------------------------
# #   
# parlist = out.h0$parlist
# z = NULL
# tauset = c(0.1, 0.3, 0.5)
# ninits = 10
# epsilon.short = 1e-02
# epsilon = 1e-08
# maxit.short = 500
# maxit = 2000
# verb = FALSE
# parallel = FALSE
# cl = NULL
# k_max = 3
# 
# 
# y <- data$Y
# warn <- options(warn = -1) # Turn off warnings
# t <- nrow(y)
# n <- ncol(y)
# q1 <- ncol(x) + 1
# m <- length(parlist$alpha)
# 
# if (!is.null(z)) {
#   z <- as.matrix(z)
#   p <- ncol(z)
# } else {
#   p <- 0
# }
# 
# # determine if the model is AR1
# if (is.null(data.0)){
#   model.ar1 <- FALSE
#   q1.0 <- 0
# 
# } else {
#   model.ar1 <- TRUE
#   y0 <- data.0$Y
#   x0 <- data.0$X
#   z0 <- data.0$Z
#   z.init <- cbind(y0,x0,z0) # use to determine the mixture probability
#   q1.0 <- ncol(z.init) + 1 # the dimension of gamma_0
#   if (!is.null(x0)) {
#     q1.0 <- ncol(x0) + 1
#   } else {
#     q1.0 <- 1
#   }
#   if (!is.null(z0)) {
#     p.0 <- ncol(z)
#   } else {
#     p.0 <- 0
#   }
# }
# 
# coefficient.all <- matrix(0, nrow = m * length(tauset), ncol = ((q1 + q1.0 + 3) * (m + 1) + p + p.0))
# 
# 
# mubeta0 <- parlist$mubeta
# mu0h <- c(-1e+10,mubeta0[1,],1e+10)        # m+2 by 1
# 
# sigma0 <- parlist$sigma
# sigma0h <- c(sigma0[1:h],sigma0[h:m])  # m+1 by 1
# gam0 <- parlist$gam
# 
# if (model.ar1){
#   # for AR1 model
#   mubeta0.0 <- parlist$mubeta0
#   if (q1.0 > 1){
#     mubeta0.0h <- c(-1e+10,mubeta0.0[1,],1e+10)        # m+2 by 1
#   }
#   else{
#     mubeta0.0h <- c(-1e+10,mubeta0.0,1e+10)        # m+2 by 1
#   }
#   sigma0.0 <- parlist$sigma0
#   sigma0.0h <- c(sigma0.0[1:h],sigma0.0[h:m])  # m+1 by 1
# }
# 
# h <- 1
# m1 <- m + 1
# tmp <- regpanelmixPhiInit(y = y, x = x, z = z, parlist = out.h0$parlist, h = h, tau =tau, ninits = ninits.short, model.ar1=TRUE, z.init =z.init)
# 
# 
# z.init <- as.matrix(z.init)
# # short EM
# b0 <- rbind(tmp$alpha, tmp$mubeta, tmp$sigma, tmp$gam, tmp$gamma0)
# 
# 
# out.short <- cppRegPanelmixPMLEAR1(b0, y, x, ztilde, cbind(1,z.init), mu0h, sigma0h, m1, p, t, an, maxit, ninits, epsilon, tau, h, k)
# 
# 
# out.short <- cppRegPanelmixPMLEAR1(b0, y, x, ztilde, cbind(1,z.init), mu0h, sigma0h, m1, p, t, an, maxit.short, ninits.short, epsilon.short, tau, h, k)
# # 

# Finalize if the procedure work
out.h1.h <- regpanelmixMaxPhi(y=data$Y,x=data$X,
                  parlist=out.h0$parlist,
                  an=(an), an_0 = (an_0),
                  parallel = FALSE, data.0 = data.0)
  
print(Sys.time() - t)

# t <- Sys.time()
# regpanelmixCritBootAR1(y = data$Y, x = data$X,
#                        parlist = out.h0$parlist,
#                        z = data$Z, parallel = FALSE,
#                        data.0 = data.0, an=( an),
#                        an_0 = ( an_0), nbtsp = 5)$crit
# print(Sys.time() - t)
# 

# 2. Debug for regpanelmixCritBootAR1
# --------------------------------



#   Logistic regression
# --------------------------------

# Example usage:
# Generate some example data
# set.seed(123)
# X <- matrix(rnorm(100 * 2), ncol = 2)
# y <- ifelse(X[, 1] + X[, 2] + rnorm(100) > 0, 1, 0)
# 
# # Add column of ones for the intercept term
# X <- cbind(1, X)
# 
# # Fit logistic regression
# alpha <- 0.01
# iterations <- 1000
# theta <- logistic_regression(X, y, alpha, iterations)
# 
# # Print the coefficients
# print(theta)

