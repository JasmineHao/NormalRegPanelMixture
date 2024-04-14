# This script replicates table 7, 10, 11 - 14 in the paper. 
library(stargazer)
library(ggplot2)
library(reshape)
library(NormalRegPanelMixture)
# library(normalregMix)

options('nloptr.show.inequality.warning'=FALSE)

##################################################
#1. food; 2. textile; 3. wood; 4. paper; 5. chemical; 6. petro; 
#7.plastic; 8. ceramics; 9. steel; 10. othermetal;
#11. metal product; 12. machine; 13. electronics; 
#14. transportation
#equipment; 15. precision instrument; 16. other;
##################################################

ind_list <- c("Food","Textile", "Wood","Paper","Petro" , "Chemical","Plastic","Ceramics","Steel","Othermetal",
              "Metal product","Machine","Electronics",
              "Transportation equipment","Precision instrument",
              "Other")


df <- read.csv(file="data/data_production_function_missing2zero.csv")
df[df==0] <- NA
#Function
#source("C:/Users/Jasmine/Dropbox/GNR/R/productionEstimation.R")
df$t <- df$year
df <- df[order(df$id,df$t),]
ind12 <- subset(df,industry_2==12)
ind12 <- ind12[order(ind12$id,ind12$t),]


#Descriptive data for ind12
stargazer(ind12,type="text")

ind.code <- c(5,13,12)


estimate.LR.df.1 <- matrix(0,nr=length(ind.code),nc=10)
rownames(estimate.LR.df.1) <- ind_list[ind.code]
colnames(estimate.LR.df.1) <- c("M=1","M=2","M=3","M=4","M=5","M=6","M=7","M=8","M=9","M=10")
estimate.LR.df.2 <- matrix(0,nr=length(ind.code),nc=10)
rownames(estimate.LR.df.2) <- ind_list[ind.code]
colnames(estimate.LR.df.2) <- c("M=1","M=2","M=3","M=4","M=5","M=6","M=7","M=8","M=9","M=10")
estimate.LR.df.3 <- matrix(0,nr=length(ind.code),nc=10)
rownames(estimate.LR.df.3) <- ind_list[ind.code]
colnames(estimate.LR.df.3) <- c("M=1","M=2","M=3","M=4","M=5","M=6","M=7","M=8","M=9","M=10")
estimate.LR.df.4 <- matrix(0,nr=length(ind.code),nc=10)
rownames(estimate.LR.df.4) <- ind_list[ind.code]
colnames(estimate.LR.df.4) <- c("M=1","M=2","M=3","M=4","M=5","M=6","M=7","M=8","M=9","M=10")
estimate.LR.df.5 <- matrix(0,nr=length(ind.code),nc=10)
rownames(estimate.LR.df.5) <- ind_list[ind.code]
colnames(estimate.LR.df.5) <- c("M=1","M=2","M=3","M=4","M=5","M=6","M=7","M=8","M=9","M=10")
count = 1

######################################################
#For panel data
######################################################

for (each.code in ind.code){
  t <- Sys.time()
  ind.each <- subset(df,industry_2==each.code)
  ind.name <- ind_list[each.code]
  
  m.share <- cast(ind.each,id ~ year,value="lnmY_it")
  
  row.names(m.share) <- m.share$id 
  m.share <- m.share[,!(colnames(m.share)=="id")] 
  T.cap <- dim(m.share)[2]
  
  estimate.df <- matrix(0,nr=5,nc=10)
  crit.df <- matrix(0,nr=5,nc=10)
  crit.df.boot <- matrix(0,nr=5,nc=10)
  for (T in 3:3){
    t.start <- T.cap-T+1
    t.seq <- seq(from=t.start,to=t.start+T-1)
    m.share.t <- m.share[,t.seq]
    data <- list(Y = t(m.share.t[complete.cases(m.share.t),]), X = NULL,  Z = NULL)
    N <- dim(data$Y)[2]
    
    for (M in 1:10){
      out.h0 <- normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none")
      # phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
      #            beta = beta, N = N, T = T, M = M, p = p, q = q)
      an <-  anFormula(out.h0$parlist,M,N,T) 
      out.h1 <- normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an=an)
      
      
      lr.estimate <- 2 * max(out.h1$penloglik - out.h0$loglik)
      
      
      
      lr.crit <- try(regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z,cl=NULL , parallel = TRUE)$crit)
      if (class(lr.crit) == "try-error"){
        lr.crit <- c(0,0,0) 
        
      }
      estimate.df[T, M] <- paste("$", round(lr.estimate, 2), "^{", paste(rep("*", sum(lr.estimate > lr.crit)), collapse = ""), "}", "$", sep = "")
      crit.df[T,M] <- paste(round(lr.crit,2),collapse = ",")
      
      # lr.crit <- regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z,cl=NULL , parallel = TRUE,nrep=1000)$crit
      #lr.crit.boot <- regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, nbtsp = 199 ,parallel = FALSE)$crit
      #crit.df.boot[T,M] <- paste(round(lr.crit.boot ,2),collapse = ",")
      
      print(lr.estimate)
      print(lr.crit)
    }
  }
  colnames(estimate.df) <- c("M=1","M=2","M=3","M=4","M=5","M=6","M=7","M=8","M=9","M=10")
  rownames(estimate.df) <- c("T=1","T=2","T=3","T=4","T=5")
  
  colnames(crit.df) <- c("M=1","M=2","M=3","M=4","M=5","M=6","M=7","M=8","M=9","M=10")
  rownames(crit.df) <- c("T=1","T=2","T=3","T=4","T=5")
  
  # colnames(crit.df.boot) <- c("M=1","M=2","M=3","M=4","M=5")
  # rownames(crit.df.boot) <- c("T=1","T=2","T=3","T=4","T=5")
  
  
  estimate.LR.df.1[count,] <- estimate.df[1,]
  estimate.LR.df.2[count,] <- estimate.df[2,]
  estimate.LR.df.3[count,] <- estimate.df[3,]
  estimate.LR.df.4[count,] <- estimate.df[4,]
  estimate.LR.df.5[count,] <- estimate.df[5,]
  count = count + 1
  
  print("*************************************")
  print(paste("Finished", ind.name))
  print( Sys.time() - t)
  print("*************************************")
  
  sink(paste("results/Empirical/Japan_CritBoot",ind.name,".txt"))
  
  stargazer(ind.each,type="latex",title=paste("Descriptive data for ",ind.name, " industry in Japan"))
  print(paste("Estimate LR for ",ind.name))
  print(estimate.df)
  stargazer(crit.df,title=paste("Critical Values Asymptotics",each.code))
  # regpanelmixMEMtest(y = data$Y,x=NULL,t=5,m=2,crit.method="none")
  #stargazer(crit.df.boot,title=paste("Critical Values Bootstrapped",each.code))
  sink()
}






rownames(estimate.LR.df.2) <- ind_list[ind.code]
rownames(estimate.LR.df.3) <- ind_list[ind.code]
rownames(estimate.LR.df.4) <- ind_list[ind.code]
rownames(estimate.LR.df.5) <- ind_list[ind.code]


sink("results/Empirical/Japan_resultBoot.txt")
stargazer(estimate.LR.df.1)
stargazer(estimate.LR.df.2)
stargazer(estimate.LR.df.3)
stargazer(estimate.LR.df.4)
stargazer(estimate.LR.df.5)
sink()
