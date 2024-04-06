library(stargazer)
library(ggplot2)
library(reshape)
# library(normalregMix)
library(foreign)
library(NormalRegPanelMixture)
options('nloptr.show.inequality.warning'=FALSE)
options(warn = -1)

set.seed(123)

#df <- readRDS("/home/haoyu/NormalRegPanelMixture/data/ChileanClean.rds")

df <- readRDS("data/ChileanClean.rds")
df <- df[df$WL > 0, ]
df <- df[df$K > 0, ]
df['si_sl'] <- df$WI /(df$WI + df$WL)
df['si_sl'] <- log(df['si_sl'])

cl <- makeCluster(16)

ind.code <- c(311,381,321,322,331,356,342,382,352,369,324)
ind.code <- c(311,381,321)

ind.names <- c()
for (each.code in ind.code){
  ind.each <- subset(df,ciiu_3d==each.code)
  ind.name <- ind.each$ciiu3d_descr[1]
  ind.names <- append(ind.names,ind.name)
}
ind.count <- length(ind.code)



estimate.LR.df.2 <- matrix(0,nr=length(ind.code),nc=15)
rownames(estimate.LR.df.2) <- ind.names
colnames(estimate.LR.df.2) <- c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
estimate.LR.df.3 <- matrix(0,nr=length(ind.code),nc=15)
rownames(estimate.LR.df.3) <- ind.names
colnames(estimate.LR.df.3) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
estimate.LR.df.4 <- matrix(0,nr=length(ind.code),nc=15)
rownames(estimate.LR.df.4) <- ind.names
colnames(estimate.LR.df.4) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
estimate.LR.df.5 <- matrix(0,nr=length(ind.code),nc=15)
rownames(estimate.LR.df.5) <- ind.names
colnames(estimate.LR.df.5) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")

AIC.df.2 <- matrix(0,nr=length(ind.code),nc=15)
rownames(AIC.df.2) <- ind.names
colnames(AIC.df.2) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
AIC.df.3 <- matrix(0,nr=length(ind.code),nc=15)
rownames(AIC.df.3) <- ind.names
colnames(AIC.df.3) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
AIC.df.4 <- matrix(0,nr=length(ind.code),nc=15)
rownames(AIC.df.4) <- ind.names
colnames(AIC.df.4) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
AIC.df.5 <- matrix(0,nr=length(ind.code),nc=15)
rownames(AIC.df.5) <- ind.names
colnames(AIC.df.5) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")


BIC.df.2 <- matrix(0,nr=length(ind.code),nc=15)
rownames(BIC.df.2) <- ind.names
colnames(BIC.df.2) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
BIC.df.3 <- matrix(0,nr=length(ind.code),nc=15)
rownames(BIC.df.3) <- ind.names
colnames(BIC.df.3) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
BIC.df.4 <- matrix(0,nr=length(ind.code),nc=15)
rownames(BIC.df.4) <- ind.names
colnames(BIC.df.4) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
BIC.df.5 <- matrix(0,nr=length(ind.code),nc=15)
rownames(BIC.df.5) <- ind.names
colnames(BIC.df.5) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")



crit.LR.df.2 <- matrix(0,nr=length(ind.code),nc=15)
rownames(crit.LR.df.2) <- ind.names
colnames(crit.LR.df.2) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
crit.LR.df.3 <- matrix(0,nr=length(ind.code),nc=15)
rownames(crit.LR.df.3) <- ind.names
colnames(crit.LR.df.3) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
crit.LR.df.4 <- matrix(0,nr=length(ind.code),nc=15)
rownames(crit.LR.df.4) <- ind.names
colnames(crit.LR.df.4) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
crit.LR.df.5 <- matrix(0,nr=length(ind.code),nc=15)
rownames(crit.LR.df.5) <- ind.names
colnames(crit.LR.df.5) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")


######################################################
#For panel data
######################################################

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
  
  coef.df <- matrix(0,nr=5,nc=15)
  estimate.df <- matrix(0,nr=5,nc=15)
  AIC.df <- matrix(0,nr=5,nc=15)
  BIC.df <- matrix(0,nr=5,nc=15)
  crit.df <- matrix(0,nr=5,nc=15)
  result.df <- matrix(0,nr=5,nc=15)
  ######################################################
  #For panel data
  ######################################################
  
  
  for (T in 2:5){
    t.start <- T.cap-T+1
    t.seq <- seq(from=t.start,to=t.start+T-1)
    
    ind.each.t <- ind.each[ind.each$year >= t.start,]
    ind.each.y <- cast(ind.each.t[,c("id","year","si_sl")],id ~ year,value="si_sl")
    id.list    <- ind.each.y[complete.cases(ind.each.y),"id"]
    #Remove the incomplete data, need balanced panel
    ind.each.t <- ind.each.t[ind.each.t$id %in% id.list,]
    ind.each.t <- ind.each.t[order(ind.each.t$id,ind.each.t$year),]
    #Reshape the Y 
    ind.each.y <- cast(ind.each.t[,c("id","year","si_sl")],id ~ year,value="si_sl")
    ind.each.y <- ind.each.y[,colnames(ind.each.y)!="id"]
    
    ind.each.y <- (ind.each.y - mean(ind.each.t$si_sl))/(sd(ind.each.t$si_sl))
    ind.each.xk <- (ind.each.t$lnk - mean(ind.each.t$lnk))/(sd(ind.each.t$lnk))
    ind.each.xl <- (ind.each.t$lnl - mean(ind.each.t$lnl))/(sd(ind.each.t$lnl))
    
    data <- list(Y = t(ind.each.y), X = data.frame(col1=ind.each.xk,col2=ind.each.xl),  Z = NULL)
    
    
    N <- dim(data$Y)[2]
    
    h1.coefficient = NULL
    estimate.crit <- 1
    for (M in 1:10){
      # Estimate the null model
      out.h0 <- regpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none",in.coefficient=h1.coefficient)
      an <- anFormula(out.h0$parlist,M,N,T,q=2)
      print("-----------------------------------------")
      print(paste("T=",T,"M = ",M,"an=",an))
      if (is.na(an)){
        an <- 1.0
      }
      # Estimate the alternative model
      out.h1 <- regpanelmixMaxPhi(y=data$Y,x=data$X, z = data$Z,parlist=out.h0$parlist,an=an)
      h1.parlist = out.h1$parlist
      
      lr.estimate <- 2 * max(out.h1$penloglik - out.h0$loglik)
      
      # Simulate the asymptotic distribution
      if (estimate.crit == 1){
        lr.crit <- try(regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, cl=cl,parallel = TRUE)$crit)
        if (class(lr.crit) == "try-error"){
          #lr.crit <- c(0,0,0)
          lr.crit <- regpanelmixCritBoot(y=data$Y, x=data$X,  parlist=out.h0$parlist, z = data$Z, cl=cl,parallel = TRUE)$crit
        }
      }
      # Store the estimation results
      coef.df[T,M] <- paste(paste(names(out.h0$coefficients), collapse = ","), paste(out.h0$coefficients, collapse = ","), collapse = ",")
      estimate.df[T,M] <- paste('$',round(lr.estimate,2),'^{',paste(rep('*',sum(lr.estimate > lr.crit)),  collapse = ""),'}','$', sep = "")
      AIC.df[T,M] <- round(out.h0$aic,2)
      BIC.df[T,M] <- round(out.h0$bic,2)
      
      
      crit.df[T,M] <- paste(round(lr.crit,2),collapse = ",")
      # If fail to reject the test, break the loop
      print(lr.estimate)
      print(lr.crit)
      
      if (sum(lr.estimate > lr.crit) < 3){
        estimate.crit <- 0
      }
    }
  }
  ###################################################################
  #     Output
  ###################################################################
  count = count + 1
  print("*************************************")
  print(paste("Finished", ind.name))
  print( Sys.time() - t)
  print("*************************************")
  
  estimate.LR.df.2[count,] <- estimate.df[2,]
  estimate.LR.df.3[count,] <- estimate.df[3,]
  estimate.LR.df.4[count,] <- estimate.df[4,]
  estimate.LR.df.5[count,] <- estimate.df[5,]
  AIC.df.2[count,] <- AIC.df[2,]
  AIC.df.3[count,] <- AIC.df[3,]
  AIC.df.4[count,] <- AIC.df[4,]
  AIC.df.5[count,] <- AIC.df[5,]
  
  BIC.df.2[count,] <- BIC.df[2,]
  BIC.df.3[count,] <- BIC.df[3,]
  BIC.df.4[count,] <- BIC.df[4,]
  BIC.df.5[count,] <- BIC.df[5,]
  crit.LR.df.2[count,] <- crit.df[2,]
  crit.LR.df.3[count,] <- crit.df[3,]
  crit.LR.df.4[count,] <- crit.df[4,]
  crit.LR.df.5[count,] <- crit.df[5,]
  
  colnames(estimate.df) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
  rownames(estimate.df) <- c("T=1","T=2","T=3","T=4","T=5")
  
  colnames(crit.df) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
  rownames(crit.df) <- c("T=1","T=2","T=3","T=4","T=5")
  
}

# write.csv(cbind(estimate.LR.df.3,AIC.df.3),file="/home/haoyu/results/Japan/resultLR3_regressor.csv")
#     write.csv(cbind(estimate.LR.df.4,AIC.df.4),file="/home/haoyu/results/Japan/resultLR4_regressor.csv")
# write.csv(cbind(estimate.LR.df.5,AIC.df.5),file="/home/haoyu/results/Japan/resultLR5_regressor.csv")

count <- length(ind.names)
df.2 <- data.frame(matrix('-',nrow=3*length(ind.names),ncol=15))
df.2[ 3* 1:count -2,] <- estimate.LR.df.2
df.2[ 3* 1:count -1,] <- AIC.df.2
df.2[ 3* 1:count ,] <- BIC.df.2
rownames(df.2)[ 3* 1:count -2] <- rownames(estimate.LR.df.2)
colnames(df.2) <- colnames(estimate.LR.df.2)


df.3 <- data.frame(matrix('-',nrow=3*length(ind.names),ncol=15))
df.3[ 3* 1:count -2,] <- estimate.LR.df.3
df.3[ 3* 1:count -1,] <- AIC.df.3
df.3[ 3* 1:count ,] <- BIC.df.3
rownames(df.3)[ 3* 1:count -2] <- rownames(estimate.LR.df.3)
colnames(df.3) <- colnames(estimate.LR.df.3)

df.4 <- data.frame(matrix('-',nrow=3*length(ind.names),ncol=15))
df.4[ 3* 1:count -2,] <- estimate.LR.df.4
df.4[ 3* 1:count -1,] <- AIC.df.4
df.4[ 3* 1:count ,] <- BIC.df.4
rownames(df.4)[ 3* 1:count -2] <- rownames(estimate.LR.df.4)
colnames(df.4) <- colnames(estimate.LR.df.4)


df.5 <- data.frame(matrix('-',nrow=3*length(ind.names),ncol=15))
df.5[ 3* 1:count -2,] <- estimate.LR.df.5
df.5[ 3* 1:count -1,] <- AIC.df.5
df.5[ 3* 1:count ,] <- BIC.df.5
rownames(df.5)[ 3* 1:count -2] <- rownames(estimate.LR.df.5)
colnames(df.5) <- colnames(estimate.LR.df.5)


write.csv(df.2,file="results/Chile/resultLR2_beta_ratio_normed_kl.csv")
write.csv(df.3,file="results/Chile/resultLR3_beta_ratio_normed_kl.csv")
write.csv(df.4,file="results/Chile/resultLR4_beta_ratio_normed_kl.csv")
write.csv(df.5,file="results/Chile/resultLR5_beta_ratio_normed_kl.csv")

