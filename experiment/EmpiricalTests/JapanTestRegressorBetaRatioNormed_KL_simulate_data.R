# This script replicates table 7, 10, 11 - 14 in the paper.
library(stargazer)
library(reshape)
library(NormalRegPanelMixture)
# library(Hmisc)
set.seed(123)
nrep <- 10
# library(normalregMix)

options('nloptr.show.inequality.warning'=FALSE)
cl <- makeCluster(16)
##################################################
#1. food; 2. textile; 3. wood; 4. paper; 5. chemical; 6. petro;
#7.plastic; 8. ceramics; 9. steel; 10. othermetal;
#11. metal product; 12. machine; 13. electronics;
#14. transportation
#equipment; 15. precision instrument; 16. other;
##################################################

ind_list <- c("food","textile", "wood","paper", "chemical",
              "petro","plastic","ceramics","steel","othermetal",
              "metal product","machine","electronics",
              "transportation equipment","precision instrument",
              "other")


# df <- readRDS("/home/haoyu/NormalRegPanelMixture/data/JapanClean.rds")
df <- read.csv(file="data/data_production_function_missing2zero_ver3.csv")
df[df==0] <- NA
#Function
#source("C:/Users/Jasmine/Dropbox/GNR/R/productionEstimation.R")
df$t <- df$year
df <- df[order(df$id,df$t),]
df['si_sl'] <- exp(df$lnmY_it) /(exp(df$lnmY_it) + exp(df$lnlY_it))
df <- df[df$si_sl > 0, ]
df <- df[df$si_sl < 1, ]
df['si_sl'] = log(df['si_sl'])



ind.code <- c(5,13,12,14,1)
ind.code <- c(5,13,12)
ind.names <- c()
for (each.code in ind.code){
  ind.name <- ind_list[each.code]
  ind.names <- append(ind.names,ind.name)
}

ind.count <- length(ind.code)

######################################################
#Define data generating function for empirical parameters
######################################################

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


######################################################
#For panel data
######################################################

estimate.LR.df.H <- matrix(0,nr=length(ind.code)*5,nc=15)
colnames(estimate.LR.df.H) <- c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10", "M=11","M=12","M=13","M=14","M=15")
estimate.LR.df.M <- matrix(0,nr=length(ind.code)*5,nc=15)
colnames(estimate.LR.df.M) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10", "M=11","M=12","M=13","M=14","M=15")
estimate.LR.df.L <- matrix(0,nr=length(ind.code)*5,nc=15)
colnames(estimate.LR.df.L) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10", "M=11","M=12","M=13","M=14","M=15")

AIC.df.M <- matrix(0,nr=length(ind.code)*5,nc=15)
colnames(AIC.df.M) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10", "M=11","M=12","M=13","M=14","M=15")


BIC.df.M <- matrix(0,nr=length(ind.code)*5,nc=15)
colnames(BIC.df.M) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10", "M=11","M=12","M=13","M=14","M=15")



count = 0
for (each.code in ind.code){
  t <- Sys.time()
  ind.each <- subset(df,industry_2==each.code)
  ind.each <- ind.each[,c("id","year","si_sl","k_it","l_it")]
  ind.each <- ind.each[complete.cases(ind.each),]
  ind.each['lnk'] <- ind.each['k_it']
  ind.each['lnl'] <- ind.each['l_it']
  each.name <- ind_list[each.code]
  
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
  
  T <- 5
  # for (T in 2:5){
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
  
  data_0 <- list(Y = t(ind.each.y), X = data.frame(col1=ind.each.xk,col2=ind.each.xl),  Z = NULL)
  
  
  N <- dim(data$Y)[2]
  
  h1.coefficient = NULL
  M_0 <- 5
  out.h0_0 <- regpanelmixPMLE(y=data_0$Y,x=data_0$X, z = data_0$Z,m=M_0,vcov.method = "none",in.coefficient=h1.coefficient)
  phi <- out.h0_0$parlist
  phi$T <- dim(data_0$Y)[1]
  phi$N <- dim(data_0$Y)[2]
  phi$M <- 5
  phi$p <- 0
  phi$q <- 2
  phi$gamma <- NULL
  phi$mu <- phi$mubeta[1,]
  phi$beta <- t(phi$mubeta[2:3,])
  phi$X <- data_0$X
  data.phi.pair <- GenerateSample(phi,3)
  Data <- data.phi.pair$Data
    
  for (i in 1:3){
    data <- Data[,i]
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
      out.h1.M <- regpanelmixMaxPhi(y=data$Y,x=data$X, z = data$Z,parlist=out.h0$parlist,an=an)
      lr.estimate.M <- 2 * max(out.h1.M$penloglik - out.h0$loglik)
      
      out.h1.H <- regpanelmixMaxPhi(y=data$Y,x=data$X, z = data$Z,parlist=out.h0$parlist,an=100*an)
      lr.estimate.H <- 2 * max(out.h1.H$penloglik - out.h0$loglik)
      
      out.h1.L <- regpanelmixMaxPhi(y=data$Y,x=data$X, z = data$Z,parlist=out.h0$parlist,an=0.01*an)
      lr.estimate.L <- 2 * max(out.h1.L$penloglik - out.h0$loglik)
      
      # Simulate the asymptotic distribution
      if (estimate.crit == 1){
        lr.crit <- try(regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, cl=cl,parallel = TRUE)$crit)
        if (class(lr.crit) == "try-error"){
          #lr.crit <- c(0,0,0)
          lr.crit <- regpanelmixCritBoot(y=data$Y, x=data$X,  parlist=out.h0$parlist, z = data$Z, cl=cl,parallel = TRUE)$crit
        }
      }
      # Store the estimation results
      estimate.LR.df.M[(i+count*3),M] <- paste('$',round(lr.estimate.M,2),'^{',paste(rep('*',sum(lr.estimate.M > lr.crit)),  collapse = ""),'}','$', sep = "")
      estimate.LR.df.H[(i+count*3),M] <- paste('$',round(lr.estimate.H,2),'^{',paste(rep('*',sum(lr.estimate.H > lr.crit)),  collapse = ""),'}','$', sep = "")
      estimate.LR.df.L[(i+count*3),M] <- paste('$',round(lr.estimate.L,2),'^{',paste(rep('*',sum(lr.estimate.L > lr.crit)),  collapse = ""),'}','$', sep = "")
      
      AIC.df.M[(i+count*3) ,M] <- round(out.h0$aic,2)
      BIC.df.M[(i+count*3),M] <- round(out.h0$bic,2)
      
      
      print(lr.estimate.M)
      print(lr.crit)
      
      if (sum(lr.estimate.L > lr.crit) < 3){
         estimate.crit <- 0
      }
    }

  }
  count <- count + 1
  print("*************************************")
  print(paste("Finished", each.name))
  print( Sys.time() - t)
  print("*************************************")
}


###################################################################
#     Output
###################################################################

# write.csv(cbind(estimate.LR.df.3,AIC.df.3),file="/home/haoyu/results/Japan/resultLR3_regressor.csv")
#     write.csv(cbind(estimate.LR.df.4,AIC.df.4),file="/home/haoyu/results/Japan/resultLR4_regressor.csv")
# write.csv(cbind(estimate.LR.df.5,AIC.df.5),file="/home/haoyu/results/Japan/resultLR5_regressor.csv")

write.csv(estimate.LR.df.H,file="results/Japan/BetaRatioNormed_KL_simulate_data/estimate.LR.df.H.csv")
write.csv(estimate.LR.df.M,file="results/Japan/BetaRatioNormed_KL_simulate_data/estimate.LR.df.M.csv")
write.csv(estimate.LR.df.L,file="results/Japan/BetaRatioNormed_KL_simulate_data/estimate.LR.df.L.csv")
write.csv(AIC.df.M,file="results/Japan/BetaRatioNormed_KL_simulate_data/AIC.df.M.csv")
write.csv(BIC.df.M,file="results/Japan/BetaRatioNormed_KL_simulate_data/BIC.df.M.csv")

