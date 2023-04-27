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

ind.code <- c(311,381,321,322,331,356,342,382,352,369,324)
ind.code <- c(311,381,321)
ind.names <- c()
for (each.code in ind.code){
  ind.each <- subset(df,ciiu_3d==each.code)
  ind.name <- ind.each$ciiu3d_descr[1]
  ind.names <- append(ind.names,ind.name)
}
ind.count <- length(ind.code)


estimate.LR.df.2 <- matrix(0,nr=length(ind.code),nc=5)
rownames(estimate.LR.df.2) <- ind.names
colnames(estimate.LR.df.2) <- c("M=1","M=2","M=3","M=4","M=5")
estimate.LR.df.3 <- matrix(0,nr=length(ind.code),nc=5)
rownames(estimate.LR.df.3) <- ind.names
colnames(estimate.LR.df.3) <- c("M=1","M=2","M=3","M=4","M=5")
estimate.LR.df.4 <- matrix(0,nr=length(ind.code),nc=5)
rownames(estimate.LR.df.4) <- ind.names
colnames(estimate.LR.df.4) <- c("M=1","M=2","M=3","M=4","M=5")
estimate.LR.df.5 <- matrix(0,nr=length(ind.code),nc=5)
rownames(estimate.LR.df.5) <- ind.names
colnames(estimate.LR.df.5) <- c("M=1","M=2","M=3","M=4","M=5")

AIC.df.2 <- matrix(0,nr=length(ind.code),nc=5)
rownames(AIC.df.2) <- ind.names
colnames(AIC.df.2) <- c("M=1","M=2","M=3","M=4","M=5")
AIC.df.3 <- matrix(0,nr=length(ind.code),nc=5)
rownames(AIC.df.3) <- ind.names
colnames(AIC.df.3) <- c("M=1","M=2","M=3","M=4","M=5")
AIC.df.4 <- matrix(0,nr=length(ind.code),nc=5)
rownames(AIC.df.4) <- ind.names
colnames(AIC.df.4) <- c("M=1","M=2","M=3","M=4","M=5")
AIC.df.5 <- matrix(0,nr=length(ind.code),nc=5)
rownames(AIC.df.5) <- ind.names
colnames(AIC.df.5) <- c("M=1","M=2","M=3","M=4","M=5")

crit.LR.df.2 <- matrix(0,nr=length(ind.code),nc=5)
rownames(crit.LR.df.2) <- ind.names
colnames(crit.LR.df.2) <- c("M=1","M=2","M=3","M=4","M=5")
crit.LR.df.3 <- matrix(0,nr=length(ind.code),nc=5)
rownames(crit.LR.df.3) <- ind.names
colnames(crit.LR.df.3) <- c("M=1","M=2","M=3","M=4","M=5")
crit.LR.df.4 <- matrix(0,nr=length(ind.code),nc=5)
rownames(crit.LR.df.4) <- ind.names
colnames(crit.LR.df.4) <- c("M=1","M=2","M=3","M=4","M=5")
crit.LR.df.5 <- matrix(0,nr=length(ind.code),nc=5)
rownames(crit.LR.df.5) <- ind.names
colnames(crit.LR.df.5) <- c("M=1","M=2","M=3","M=4","M=5")
# for (each.code in ind.code){
#   t <- Sys.time()
#   ind.each <- subset(df,ciiu_3d==each.code)
#   ind.name <- ind.each$ciiu3d_descr[1]
#   ind.each$lny <- log(ind.each$GO)
#   ind.each$lnm <- log(ind.each$WI)
#   ind.each$lnl <- log(ind.each$L)
#   ind.each$lnk <- log(ind.each$K)
#   ######################################################
#   # Describe the data
#   ######################################################

#   desc.each <- ind.each[ind.each$L != 0 ,c("si","lny","lnm","lnl","lnk")]
#   desc.each <- desc.each[complete.cases(desc.each),]



#   ######################################################
#   # Select the data out
#   ######################################################
#   m.share <- cast(ind.each,id ~ year,value="si")#Collapse the dataframe into panel form , year against firm id
#   row.names(m.share) <- m.share$id
#   m.share <- m.share[,!(colnames(m.share)=="id")]
#   T.cap <- dim(m.share)[2]

#   estimate.df <- matrix(0,nr=5,nc=5)
#   crit.df <- matrix(0,nr=5,nc=5)
#   result.df <- matrix(0,nr=5,nc=5)
#   # crit.df.boot <- matrix(0,nr=5,nc=5)

#   ######################################################
#   #For cross-sectional data
#   ######################################################
#   m.share.t <- m.share[,T.cap]
#   m.share.t <- m.share.t[complete.cases(m.share.t)]
#   N <- length(m.share.t)
#   for (M in 1:5){

#     out.h0 <- normalmixPMLE(y = m.share.t,m=M)
#     an <- normalregMix::anFormula(out.h0$parlist,M,N)
#     print(paste("T=",1,"M = ",M,"an=",an))
#     out.h1 <- normalmixMaxPhi(y=m.share.t,parlist = out.h0$parlist,an=an)
#     #Obtain critical value
#     lr.crit <- normalmixCrit(y=m.share.t,parlist = out.h0$parlist)$crit
#     crit.df[1,M] <- paste(round(lr.crit,2),collapse = ",")
#     # lr.crit.boot <- normalmixCritBoot(y=m.share.t,parlist = out.h0$parlist,parallel = FALSE,nbtsp = 200)$crit
#     # crif.df.boot[1,M] <- paste(round(lr.crit.boot,2),collapse = ",")
#     estimate.df[1,M] <- 2 * max(out.h1$penloglik - out.h0$loglik)
#     result.df[1,M] <- (2 * max(out.h1$penloglik - out.h0$loglik) > lr.crit[2])

#   }

count = 0

for (each.code in ind.code){
  t <- Sys.time()
  ind.each <- subset(df,ciiu_3d==each.code)
  ind.name <- ind.each$ciiu3d_descr[1]
  ind.each$lny <- log(ind.each$GO)
  ind.each$lnm <- log(ind.each$WI)
  ind.each$lnl <- log(ind.each$L)
  ind.each$lnk <- log(ind.each$K)
  ######################################################
  # Describe the data
  ######################################################

  desc.each <- ind.each[ind.each$L != 0 ,c("si","lny","lnm","lnl","lnk")]
  # desc.each <- desc.each[complete.cases(desc.each),]

  ######################################################
  # Select the data out
  ######################################################
  m.share <- cast(ind.each,id ~ year,value="si")#Collapse the dataframe into panel form , year against firm id
  row.names(m.share) <- m.share$id
  m.share <- m.share[,!(colnames(m.share)=="id")]
  m.share <- m.share[complete.cases(m.share),]
  T.cap <- dim(m.share)[2]

  estimate.df <- matrix(0,nr=5,nc=5)
  AIC.df <- matrix(0,nr=5,nc=5)
  crit.df <- matrix(0,nr=5,nc=5)
  result.df <- matrix(0,nr=5,nc=5)
  ######################################################
  #For panel data
  ######################################################


  for (T in 2:5){
    t.start <- T.cap-T+1
    t.seq <- seq(from=t.start,to=t.start+T-1)
    m.share.t <- m.share[,t.seq]
    data <- list(Y = t(m.share.t[complete.cases(m.share.t),]), X = NULL,  Z = NULL)
    N <- dim(data$Y)[2]

    h1.coefficient = NULL
    for (M in 1:5){
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
      AIC.df[T,M] <- paste(round(out.h0$aic,2),round(out.h0$bic,2))
      crit.df[T,M] <- paste(round(lr.crit,2),collapse = ",")
      # If fail to reject the test, break the loop
      if (sum(lr.estimate > lr.crit) < 1) break

      print(lr.estimate)
      print(lr.crit)
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
  
  crit.LR.df.2[count,] <- crit.df[2,]
  crit.LR.df.3[count,] <- crit.df[3,]
  crit.LR.df.4[count,] <- crit.df[4,]
  crit.LR.df.5[count,] <- crit.df[5,]

  colnames(estimate.df) <- c("M=1","M=2","M=3","M=4","M=5")
  rownames(estimate.df) <- c("T=1","T=2","T=3","T=4","T=5")

  colnames(crit.df) <- c("M=1","M=2","M=3","M=4","M=5")
  rownames(crit.df) <- c("T=1","T=2","T=3","T=4","T=5")

}



# colnames(crit.df.boot) <- c("M=1","M=2","M=3","M=4","M=5")
# rownames(crit.df.boot) <- c("T=1","T=2","T=3","T=4","T=5")


df.2 <- data.frame(matrix('-',nrow=2*length(ind.names),ncol=5))
df.2[ 2* 1:count -1,] <- estimate.LR.df.2
df.2[ 2* 1:count,] <- AIC.df.2
rownames(df.2)[ 2* 1:count -1] <- rownames(estimate.LR.df.2)
colnames(df.2) <- colnames(estimate.LR.df.2)

df.3 <- data.frame(matrix('-',nrow=2*length(ind.names),ncol=5))
df.3[ 2* 1:count -1,] <- estimate.LR.df.3
df.3[ 2* 1:count,] <- AIC.df.3
rownames(df.3)[ 2* 1:count -1] <- rownames(estimate.LR.df.3)
colnames(df.3) <- colnames(estimate.LR.df.3)

df.4 <- data.frame(matrix('-',nrow=2*length(ind.names),ncol=5))
df.4[ 2* 1:count -1,] <- estimate.LR.df.4
df.4[ 2* 1:count,] <- AIC.df.4
rownames(df.4)[ 2* 1:count -1] <- rownames(estimate.LR.df.4)
colnames(df.4) <- colnames(estimate.LR.df.4)

df.5 <- data.frame(matrix('-',nrow=2*length(ind.names),ncol=5))
df.5[ 2* 1:count -1,] <- estimate.LR.df.5
df.5[ 2* 1:count,] <- AIC.df.5
rownames(df.5)[ 2* 1:count -1] <- rownames(estimate.LR.df.5)
colnames(df.5) <- colnames(estimate.LR.df.5)

write.csv(df.2,file="results/Chile/resultLR2.csv")
write.csv(df.3,file="results/Chile/resultLR3.csv")
write.csv(df.4,file="results/Chile/resultLR4.csv")
write.csv(df.5,file="results/Chile/resultLR5.csv")

# write.csv(df.2,file="/home/haoyu/results/Chile/resultLR2.csv")
# write.csv(df.3,file="/home/haoyu/results/Chile/resultLR3.csv")
# write.csv(df.4,file="/home/haoyu/results/Chile/resultLR4.csv")
# write.csv(df.5,file="/home/haoyu/results/Chile/resultLR5.csv")

library(xtable)
# stargazer(crit.df,title=paste("estimate",ind.name))
sink("/home/haoyu/results/Chile/result_text.txt")
xtable(estimate.LR.df.2)
xtable(crit.LR.df.2)
xtable(estimate.LR.df.3)
xtable(crit.LR.df.3)
xtable(estimate.LR.df.4)
xtable(crit.LR.df.4)
xtable(estimate.LR.df.5)
xtable(crit.LR.df.5)
sink()

sink("/home/haoyu/results/Chile/result.txt")
sink("results/Chile/result.txt")
stargazer(estimate.LR.df.2)
stargazer(AIC.df.2)
stargazer(crit.LR.df.2)

stargazer(estimate.LR.df.3)
stargazer(AIC.df.3)
stargazer(crit.LR.df.3)

stargazer(estimate.LR.df.4)
stargazer(AIC.df.4)
stargazer(crit.LR.df.4)

stargazer(estimate.LR.df.5)
stargazer(AIC.df.5)
stargazer(crit.LR.df.5)
sink()

