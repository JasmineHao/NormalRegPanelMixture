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
cl <- makeCluster(12)

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


count = 0
for (each.code in ind.code){
  t <- Sys.time()
  ind.each <- subset(df,ciiu_3d==each.code)
  ind.name <- ind.each$ciiu3d_descr[1]
  ind.each$y <- log(ind.each$GO)
  ind.each$lnm <- log(ind.each$WI)
  ind.each$lnl <- log(ind.each$L)
  ind.each$lnk <- log(ind.each$K)
  ######################################################
  # Describe the data
  ######################################################

  desc.each <- ind.each[ind.each$L != 0 ,c("si","y","lnm","lnl","lnk")]
  # desc.each <- desc.each[complete.cases(desc.each),]
  year.list <- sort(unique(ind.each$year))
  T.cap <- max(year.list)
  
  coef.df <- matrix(0,nr=5,nc=5)
  estimate.df <- matrix(0,nr=5,nc=5)
  AIC.df <- matrix(0,nr=5,nc=5)
  crit.df <- matrix(0,nr=5,nc=5)
  result.df <- matrix(0,nr=5,nc=5)
  ######################################################
  #For panel data
  ######################################################


  for (T in 3:5){
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
    
    h1.coefficient = NULL
    for (M in 1:5){
      # Estimate the null model
      data <- list(Y = t(ind.each.y), X = matrix(ind.each.t$lnm),  Z = NULL)
      
      out.h0 <- regpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none",in.coefficient=h1.coefficient)
      an <- anFormula(out.h0$parlist,M,N,T,q=1)
      print("-----------------------------------------")
      print(paste("T=",T,"M = ",M,"an=",an))
      if (is.na(an)){
        an <- 1.0
      }
      # Estimate the alternative model
      out.h1 <- regpanelmixMaxPhi(y=data$Y,x=data$X,parlist=out.h0$parlist,an=an)
      h1.parlist = out.h1$parlist

      lr.estimate <- 2 * max(out.h1$penloglik - out.h0$loglik)

      # Simulate the asymptotic distribution
      lr.crit <- regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, cl=cl,parallel = TRUE)$crit
      
      # Store the estimation results
      coef.df[T,M] <- paste(paste(names(out.h0$coefficients), collapse = ","), paste(out.h0$coefficients, collapse = ","))
      
      estimate.df[T,M] <- paste('$',round(lr.estimate,2),'^{',paste(rep('*',sum(lr.estimate > lr.crit)),  collapse = ""),'}','$', sep = "")
      AIC.df[T,M] <- out.h0$aic
      crit.df[T,M] <- paste(round(lr.crit,2),collapse = ",")
      # If fail to reject the test, break the loop
      print(lr.estimate)
      print(lr.crit)
      if (sum(lr.estimate > lr.crit) < 1) break


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

count <- length(ind.names)
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

# Combine the data frames
combined_df <- rbind(
  cbind(df.3, original_df = "df.3")
)

write.csv(combined_df,file="results/Empirical/Chile_combined_regressorBoot.csv")


# sink("/home/haoyu/results/Empirical/Chile_result_regressor.txt")
sink("results/Empirical/Chile_result_regressor_boot.txt")
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
