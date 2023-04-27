library(stargazer)
library(ggplot2)
library(reshape)
# library(normalregMix)
library(foreign)
library(NormalRegPanelMixture)
options('nloptr.show.inequality.warning'=FALSE)
options(warn = -1)

cl <- makeCluster(6)
# df <- readRDS("/home/haoyu/NormalRegPanelMixture/data/ChileanClean.rds")

df <- readRDS("data/ChileanClean.rds")

each.code <- 311

ind.each <- subset(df,ciiu_3d==each.code)
ind.name <- ind.each$ciiu3d_descr[1]
ind.each$lny <- log(ind.each$GO)
ind.each$lnm <- log(ind.each$WI)
ind.each$lnl <- log(ind.each$L)
ind.each$lnk <- log(ind.each$K)

m.share <- cast(ind.each,id ~ year,value="si")#Collapse the dataframe into panel form , year against firm id
row.names(m.share) <- m.share$id
m.share <- m.share[,!(colnames(m.share)=="id")]
m.share <- m.share[complete.cases(m.share),]
T.cap <- dim(m.share)[2]

estimate.df <- matrix(0,nr=5,nc=5)
AIC.df <- matrix(0,nr=5,nc=5)
crit.df <- matrix(0,nr=5,nc=5)
result.df <- matrix(0,nr=5,nc=5)

T <- 2

t.start <- T.cap-T+1
t.seq <- seq(from=t.start,to=t.start+T-1)
m.share.t <- m.share[,t.seq]
data <- list(Y = t(m.share.t[complete.cases(m.share.t),]), X = NULL,  Z = NULL)
N <- dim(data$Y)[2]

h1.coefficient = NULL
M <- 2
{
# Estimate the null model
out.h0 <- normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none",in.coefficient=h1.coefficient)
an <- anFormula(out.h0$parlist,M,N,T)
print("-----------------------------------------")
print(paste("T=",T,"M = ",M,"an=",an))
if (is.na(an)){
  an <- 1.0
}
# # Estimate the alternative model
out.h1 <- normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an=an)
h1.parlist = out.h1$parlist

lr.estimate <- 2 * max(out.h1$penloglik - out.h0$loglik)

# Simulate the asymptotic distribution
lr.crit <- try(regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, cl=cl,parallel = TRUE)$crit)
if (class(lr.crit) == "try-error"){
  lr.crit <- regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, cl=cl,parallel = TRUE)$crit
}

# Store the estimation results
estimate.df[T,M] <- paste('$',round(lr.estimate,2),'^{',paste(rep('*',sum(lr.estimate > lr.crit)),  collapse = ""),'}','$', sep = "")
AIC.df[T,M] <- out.h0$aic
crit.df[T,M] <- paste(round(lr.crit,2),collapse = ",")
# If fail to reject the test, break the loop
if (sum(lr.estimate > lr.crit) < 1) break

print(lr.estimate)
print(lr.crit)
print(out.h0$aic)
print(out.h0$loglik)
}



#-------------------------------------------------

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
  b  <- matrix(theta[3:4],ncol=1,nrow=2)
  c  <- theta[5:6]
  ll <- -nt * M_LN_SQRT_2PI
  
  x1 <- matrix(rep(1,length(y)))
  x1b <- cbind(x1*b[1],x1*b[2])
  
  eps <- t(t(yrep - (x1b) )/ c)
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
  b  <- matrix(theta[3:4],ncol=1,nrow=2)
  c  <- theta[5:6]
  
  eps <- t(t(yrep - (x1b) )/ c)
  
  
  eps1 <- apply(dnorm(matrix(eps[,1],nrow=T)),2,prod)
  eps2 <- apply(dnorm(matrix(eps[,2],nrow=T)),2,prod)
  
  ll <- sum(log(a[1] * eps1 / (c[1] **T)   + a[2] * eps2 / (c[2]**T)))
  ll
  }
obj_0(theta)
density_check(theta)

# ----------------------------------------------------------------------------
# 检查cpp file写的对不对: q = 0, M = 2
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
tmp <- normalpanelmixPMLEinit(y = y, z = z, ninits = ninits.short, m = M)

sd0 <- sd(y) * sqrt((n - 1) / n)
an <- 1 / n # penalty term for variance
sigma0 <- rep(sd0, M)
mu0 <- double(M + 1) # dummy
ztilde <- matrix(0) # dummy
b0 <- as.matrix(rbind(tmp$alpha, tmp$mu, tmp$sigma, tmp$gam))
out.short <- cppnormalpanelmixPMLE(b0, y, ztilde, mu0, sigma0, M, p, t, an, maxit.short,
                                ninits.short, epsilon.short)

components <- order(out.short$penloglikset, decreasing = TRUE)[1:ninits]
b1 <- b0[ ,components] # b0 has been updated
if ((! is.null(in.coefficient)) & (dim(b1)[1]==length(in.coefficient)) ){
  b1[,components] <- in.coefficient
}

out <- cppnormalpanelmixPMLE(b1, y, ztilde, mu0, sigma0, m, p, t, an, maxit, ninits, epsilon)
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


