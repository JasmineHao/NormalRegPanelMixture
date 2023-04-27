library(stargazer)
library(ggplot2)
library(reshape)
# library(normalregMix)
library(foreign)
library(NormalRegPanelMixture)
options('nloptr.show.inequality.warning'=FALSE)
options(warn = -1)

# df <- readRDS("/home/haoyu/NormalRegPanelMixture/data/ChileanClean.rds")

df <- readRDS("data/ChileanClean.rds")
cl <- makeCluster(8)
nrep <- 500
M <- 2 #Number of Type
p <- 0 #Number of Z
q <- 1 #Number of X

ind.code <- c(311,381,321,322,331,356,342,382,352,369,324)
ind.code <- c(311,381,321)
ind.names <- c()
for (each.code in ind.code){
  ind.each <- subset(df,ciiu_3d==each.code)
  ind.name <- ind.each$ciiu3d_descr[1]
  ind.names <- append(ind.names,ind.name)
}
ind.count <- length(ind.code)

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
  Data <- replicate(nrep,generateData(alpha,mu,sigma,gamma,beta,N,T,M,p,q,X=X))
  return(list(phi=phi,Data=Data))
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
    print(k)
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

# estimate.LR.df <- matrix(0,nr=(2*length(ind.code)),nc=5)
# rownames(estimate.LR.df) <- as.vector(rbind(ind.names,ind.names))
# colnames(estimate.LR.df) <- c("T=1","T=2","T=3","T=4","T=5")

estimate.LR.df <- matrix(0,nr=length(ind.code),nc=5)
rownames(estimate.LR.df) <- ind.names
colnames(estimate.LR.df) <- c("T=1","T=2","T=3","T=4","T=5")

count = 0
for (each.code in ind.code){
  t <- Sys.time()
  ind.each <- subset(df,ciiu_3d==each.code)
  ind.name <- ind.each$ciiu3d_descr[1]
  ind.each$lny <- log(ind.each$GO)
  ind.each$lnm <- log(ind.each$WI)
  ind.each$lnl <- log(ind.each$L)
  ind.each$lnk <- log(ind.each$K)
  year.list <- sort(unique(ind.each$year))
  T.cap <- max(year.list)
  count <- count+ 1
  for (T in 2:2){
    t.start <- T.cap-T+1
    #Reshape the data so that I can apply the test
    ind.each.t <- ind.each[ind.each$year>=t.start,]
    ind.each.t <- ind.each.t[complete.cases(ind.each.t),]
    ind.each.y <- cast(ind.each.t[,c("id","year","si")],id ~ year,value="si")
    id.list    <- ind.each.y[complete.cases(ind.each.y),"id"]
    #Remove the incomplete data, need balanced panel
    ind.each.t <- ind.each.t[ind.each.t$id %in% id.list,]
    ind.each.t <- ind.each.t[order(ind.each.t$id,ind.each.t$year),]
    #Reshape the Y 
    ind.each.y <- cast(ind.each.t[,c("id","year","si")],id ~ year,value="si")
    ind.each.y <- ind.each.y[,colnames(ind.each.y)!="id"]
    N <- dim(ind.each.y)[1]
    
    data <- list(Y = t(ind.each.y), X = matrix(ind.each.t$lnm),  Z = NULL)
    
    out.h0 <- regpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none")
    an <- 0.0001 * anFormula(out.h0$parlist,M,N,T,q=1)
    print("-----------------------------------------")
    print(paste("T=",T,"M = ",M,"an=",an))
    
    phi = list(alpha = out.h0$parlist$alpha,mu = out.h0$parlist$mubeta[1,],sigma = out.h0$parlist$sigma, gamma = out.h0$parlist$gam,  beta = out.h0$parlist$mubeta[2:(q+1),], N = N, T = T, M = M, p = p, q = q, X=NULL)
    GetMisclTerm(phi)
    phi.data.pair <- GenerateSample(phi,nrep)
    result <- getEstimate(phi.data.pair$Data,phi,nrep,an,M,parlist,cl)
    nominal.size <- result$nominal.size
    print(Sys.time() - t)
    
    
    estimate.LR.df[count,T] <- nominal.size
    
  }
}
write.csv(estimate.LR.df,file="/home/haoyu/results/Chile/ChileSizeTest.csv")    