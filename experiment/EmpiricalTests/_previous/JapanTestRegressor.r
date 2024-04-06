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


estimate.LR.df.1 <- matrix(0,nr=length(ind.code),nc=5)
rownames(estimate.LR.df.1) <- ind_list[ind.code]
colnames(estimate.LR.df.1) <- c("M=1","M=2","M=3","M=4","M=5")
estimate.LR.df.2 <- matrix(0,nr=length(ind.code),nc=5)
rownames(estimate.LR.df.2) <- ind_list[ind.code]
colnames(estimate.LR.df.2) <- c("M=1","M=2","M=3","M=4","M=5")
estimate.LR.df.3 <- matrix(0,nr=length(ind.code),nc=5)
rownames(estimate.LR.df.3) <- ind_list[ind.code]
colnames(estimate.LR.df.3) <- c("M=1","M=2","M=3","M=4","M=5")
estimate.LR.df.4 <- matrix(0,nr=length(ind.code),nc=5)
rownames(estimate.LR.df.4) <- ind_list[ind.code]
colnames(estimate.LR.df.4) <- c("M=1","M=2","M=3","M=4","M=5")
estimate.LR.df.5 <- matrix(0,nr=length(ind.code),nc=5)
rownames(estimate.LR.df.5) <- ind_list[ind.code]
colnames(estimate.LR.df.5) <- c("M=1","M=2","M=3","M=4","M=5")
count = 1

######################################################
#For panel data
######################################################

for (each.code in ind.code){
  t <- Sys.time()

  ind.each <- subset(df,industry_2==each.code)
  ind.each <- ind.each[complete.cases(ind.each),]
  ind.each$lnm <- log(ind.each$m_it)
  each.name <- ind_list[each.code]
  
  year.list <- sort(unique(ind.each$year))
  ######################################################
  # Select the data out
  ######################################################
  T.cap <- max(year.list) 
  
  m.share <- cast(ind.each,id ~ year,value="lnmY_it")
  
  row.names(m.share) <- m.share$id 
  m.share <- m.share[,!(colnames(m.share)=="id")] 
  
  desc.each <- matrix(0, nr = 5, nc = 5)
  estimate.df <- matrix(0,nr=5,nc=5)
  crit.df <- matrix(0,nr=5,nc=5)
  crit.df.boot <- matrix(0,nr=5,nc=5)
  for (T in 3:3){
    t.start <- T.cap-T+1
    t.seq <- seq(from=t.start,to=t.start+T-1)
    
    ind.each.t <- ind.each[ind.each$year >= t.start,]
    ind.each.t <- ind.each.t[complete.cases(ind.each.t),]
    ind.each.y <- cast(ind.each.t[,c("id","year","lnmY_it")],id ~ year,value="lnmY_it")
    id.list    <- ind.each.y[complete.cases(ind.each.y),"id"]
    #Remove the incomplete data, need balanced panel
    ind.each.t <- ind.each.t[ind.each.t$id %in% id.list,]
    ind.each.t <- ind.each.t[order(ind.each.t$id,ind.each.t$year),]
    #Reshape the Y 
    ind.each.y <- cast(ind.each.t[,c("id","year","lnmY_it")],id ~ year,value="lnmY_it")
    ind.each.y <- ind.each.y[,colnames(ind.each.y)!="id"]
    
    
    data <- list(Y = t(ind.each.y), X = matrix(ind.each.t$lnm),  Z = NULL)
    N <- dim(data$Y)[2]
    
    for (M in 1:5){
      out.h0 <- regpanelmixPMLE(y=data$Y,x=data$X, z = NULL,m=M,vcov.method = "none")
      # phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
      #            beta = beta, N = N, T = T, M = M, p = p, q = q)
      an <-  anFormula(out.h0$parlist,M,N,T) 
      out.h1 <- regpanelmixMaxPhi(y=data$Y,x=data$X, z = NULL,parlist=out.h0$parlist,an=an)      
      
      lr.estimate <- 2 * max(out.h1$penloglik - out.h0$loglik)
      
      
      
      lr.crit <- try(regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist ,parallel = TRUE)$crit)
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
   
     desc.each[T, ] <- c(each.name, each.code, dim(m.share)[1], round(mean(as.matrix(ind.each.y)), 2), round(sd(as.matrix(ind.each.y)), 2))
  }
  colnames(estimate.df) <- c("M=1","M=2","M=3","M=4","M=5")
  rownames(estimate.df) <- c("T=1","T=2","T=3","T=4","T=5")
  
  colnames(crit.df) <- c("M=1","M=2","M=3","M=4","M=5")
  rownames(crit.df) <- c("T=1","T=2","T=3","T=4","T=5")
  
  # colnames(crit.df.boot) <- c("M=1","M=2","M=3","M=4","M=5")
  # rownames(crit.df.boot) <- c("T=1","T=2","T=3","T=4","T=5")
  
  
  estimate.LR.df.1[count,] <- estimate.df[1,]
  estimate.LR.df.2[count,] <- estimate.df[2,]
  estimate.LR.df.3[count,] <- estimate.df[3,]
  estimate.LR.df.4[count,] <- estimate.df[4,]
  estimate.LR.df.5[count,] <- estimate.df[5,]
  count = count + 1
  colnames(desc.each) <- c("Industry Name", "code", "n", "mean", "sd")
  rownames(desc.each) <- c("T=1", "T=2", "T=3", "T=4", "T=5")
  print("*************************************")
  print(paste("Finished", each.name))
  print( Sys.time() - t)
  print("*************************************")
  
  sink(paste("results/Empirical/Japan_regressorCrit",each.name,".txt"))
  stargazer(desc.each, type = "text", title = paste("Descriptives for ", each.name, each.code))
  stargazer(ind.each,type="latex",title=paste("Descriptive data for ",each.name, " industry in Japan"))
  print(paste("Estimate LR for ",each.name))
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


sink("results/Empirical/Japan_result_regressor.txt")
stargazer(estimate.LR.df.1)
stargazer(estimate.LR.df.2)
stargazer(estimate.LR.df.3)
stargazer(estimate.LR.df.4)
stargazer(estimate.LR.df.5)
sink()
