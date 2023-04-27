library(stargazer)
library(ggplot2)
library(reshape)
# library(normalregMix)
library(foreign)
library(NormalRegPanelMixture)
options('nloptr.show.inequality.warning'=FALSE)
options(warn = -1)



df <- readRDS("data/ChileanClean.rds")

# ind.code <- ind.code[1:3] #For test purpose, only use the first three industries
# ind.code <- c(352,342,369,381,321,313,341,322,390,311,351,324,356,312)

ind.code <- c(311,381,321,322,331,356,342,382,352,369,324,332,384,312,313,341,351,383)
ind.names <- c("Wood products, except furniture","Machinery, except electrical","Manufacture of furniture and fixtures, except primarily of metal","Transport equipment","Other chemicals","Printing and publishing","Other non-metallic mineral products","Fabricated metal products","Textiles","Beverages","Paper and products","Wearing apparel, except footwear","Other manufactured products","Food products","Industrial chemicals","Footwear, except rubber or plastic","Plastic products","Animal feeds, etc")

ind.count <- length(ind.code)
cl <- makeCluster(detectCores())


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

# for (each.code in ind.code){
#   t <- Sys.time()
#   ind.each <- subset(df,ciiu_3d==each.code)
#   ind.name <- ind.each$ciiu3d_descr[1]
#   ind.each$lny <- log(ind.each$GO)
#   ind.each$lnm <- log(ind.each$WI)
#   ind.each$lnl <- log(ind.each$L)
#   ind.each$lnk <- log(ind.each$K)
#   ######################################################
#   # Refine the data
#   ######################################################
#   ind.each <- ind.each[ind.each$L != 0 ,c("id","year","si","lny","lnm","lnl","lnk")]
#   year.list <- sort(unique(ind.each$year))
#   ######################################################
#   # Select the data out
#   ######################################################
#   T.cap <- max(year.list)
  
#   estimate.m.df <- matrix(0,nr=5,nc=5)
#   crit.m.df <- matrix(0,nr=5,nc=5)
#   result.m.df <- matrix(0,nr=5,nc=5)
  
#   estimate.l.df <- matrix(0,nr=5,nc=5)
#   crit.l.df <- matrix(0,nr=5,nc=5)
#   result.l.df <- matrix(0,nr=5,nc=5)
  
#   estimate.k.df <- matrix(0,nr=5,nc=5)
#   crit.k.df <- matrix(0,nr=5,nc=5)
#   result.k.df <- matrix(0,nr=5,nc=5)
  
#   estimate.all.df <- matrix(0,nr=5,nc=5)
#   crit.all.df <- matrix(0,nr=5,nc=5)
#   result.all.df <- matrix(0,nr=5,nc=5)
#   # crit.df.boot <- matrix(0,nr=5,nc=5)
  
  # ######################################################
  # #For cross-sectional data
  # ######################################################
  # ind.each.t <- ind.each[ind.each$year==T.cap,]
  # ind.each.t <- ind.each.t[complete.cases(ind.each.t),]
  # N <- length(ind.each.t)
  # for (M in 1:5){
  #   ############################################
  #   #Estimate with regressor lnM
  #   ############################################
  #   out.h0 <- normalmixPMLE(y = ind.each.t$si, x = ind.each.t$lnm ,m=M)
  #   an <- normalregMix::anFormula(out.h0$parlist,M,N)
  #   print(paste("T=",1,"M = ",M,"an=",an))
  #   out.h1 <- regmixMaxPhi(y=ind.each.t$si,x = as.matrix(ind.each.t$lnm) ,parlist = out.h0$parlist,an=an)
  #     #Obtain critical value
  #   lr.crit <- regmixCrit(y=ind.each.t$si,x = as.matrix(ind.each.t$lnm),parlist = out.h0$parlist)$crit
  #   crit.m.df[1,M] <- paste(round(lr.crit,2),collapse = ",")
  #   # lr.crit.boot <- normalmixCritBoot(y=m.share.t,parlist = out.h0$parlist,parallel = FALSE,nbtsp = 200)$crit
  #   # crif.df.boot[1,M] <- paste(round(lr.crit.boot,2),collapse = ",")
  #   estimate.m.df[1,M] <- 2 * max(out.h1$penloglik - out.h0$loglik)
  #   result.m.df[1,M] <- (2 * max(out.h1$penloglik - out.h0$loglik) > lr.crit[2])
    
  #   ############################################
  #   #Estimate with regressor lnK
  #   ############################################
  #   out.h0 <- normalmixPMLE(y = ind.each.t$si, x = ind.each.t$lnk ,m=M)
  #   an <- normalregMix::anFormula(out.h0$parlist,M,N)
  #   print(paste("T=",1,"M = ",M,"an=",an))
  #   out.h1 <- regmixMaxPhi(y=ind.each.t$si,x = as.matrix( ind.each.t$lnk) ,parlist = out.h0$parlist,an=an)
  #   #Obtain critical value
  #   lr.crit <- regmixCrit(y=ind.each.t$si,x = as.matrix( ind.each.t$lnk),parlist = out.h0$parlist)$crit
  #   crit.k.df[1,M] <- paste(round(lr.crit,2),collapse = ",")
  #   # lr.crit.boot <- normalmixCritBoot(y=m.share.t,parlist = out.h0$parlist,parallel = FALSE,nbtsp = 200)$crit
  #   # crif.df.boot[1,M] <- paste(round(lr.crit.boot,2),collapse = ",")
  #   estimate.k.df[1,M] <- 2 * max(out.h1$penloglik - out.h0$loglik)
  #   result.k.df[1,M] <- (2 * max(out.h1$penloglik - out.h0$loglik) > lr.crit[2])
  #   ############################################
  #   #Estimate with regressor lnL
  #   ############################################
  #   out.h0 <- normalmixPMLE(y = ind.each.t$si, x = ind.each.t$lnl ,m=M)
  #   an <- normalregMix::anFormula(out.h0$parlist,M,N)
  #   print(paste("T=",1,"M = ",M,"an=",an))
  #   out.h1 <- regmixMaxPhi(y=ind.each.t$si,x = as.matrix( ind.each.t$lnl) ,parlist = out.h0$parlist,an=an)
  #   #Obtain critical value
  #   lr.crit <- regmixCrit(y=ind.each.t$si,x = as.matrix( ind.each.t$lnl),parlist = out.h0$parlist)$crit
  #   crit.l.df[1,M] <- paste(round(lr.crit,2),collapse = ",")
  #   # lr.crit.boot <- normalmixCritBoot(y=m.share.t,parlist = out.h0$parlist,parallel = FALSE,nbtsp = 200)$crit
  #   # crif.df.boot[1,M] <- paste(round(lr.crit.boot,2),collapse = ",")
  #   estimate.l.df[1,M] <- 2 * max(out.h1$penloglik - out.h0$loglik)
  #   result.l.df[1,M] <- (2 * max(out.h1$penloglik - out.h0$loglik) > lr.crit[2])
    
  # }
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
  # Refine the data
  ######################################################
  ind.each <- ind.each[ind.each$L != 0 ,c("id","year","si","lny","lnm","lnl","lnk")]
  year.list <- sort(unique(ind.each$year))
  ######################################################
  # Select the data out
  ######################################################
  T.cap <- max(year.list)
  
  estimate.m.df <- matrix(0,nr=5,nc=5)
  crit.m.df <- matrix(0,nr=5,nc=5)
  result.m.df <- matrix(0,nr=5,nc=5)
  
  estimate.l.df <- matrix(0,nr=5,nc=5)
  crit.l.df <- matrix(0,nr=5,nc=5)
  result.l.df <- matrix(0,nr=5,nc=5)
  
  estimate.k.df <- matrix(0,nr=5,nc=5)
  crit.k.df <- matrix(0,nr=5,nc=5)
  result.k.df <- matrix(0,nr=5,nc=5)
  
  desc.each <-matrix(0,nr=5,nc=5)
  
  estimate.all.df <- matrix(0,nr=5,nc=5)
  crit.all.df <- matrix(0,nr=5,nc=5)
  result.all.df <- matrix(0,nr=5,nc=5)
  
  
  # crit.df.boot <- matrix(0,nr=5,nc=5)  
  ######################################################
  #For panel data
  ######################################################
  
  
  for (T in 2:5){
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
    
    for (M in 1:5){
      ############################################
      #Estimate with regressor lnM
      ############################################
      data <- list(Y = t(ind.each.y), X = matrix(ind.each.t$lnm),  Z = NULL)
      out.h0 <- regpanelmixPMLE(y=data$Y,x=data$X, z = NULL,m=M,vcov.method = "none")
      # phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
      #            beta = beta, N = N, T = T, M = M, p = p, q = q)
      an <- anFormula(out.h0$parlist,M,N,T,q=1) 
      print(paste("T=",T,"M = ",M,"an=",an))
      
      if (is.na(an)){ 
        an <- 1.0
      }
      out.h1 <- regpanelmixMaxPhi(y=data$Y,x=data$X, z = NULL,parlist=out.h0$parlist,an=an)
      lr.estimate <- 2 * max(out.h1$penloglik - out.h0$loglik)
      
      lr.crit <- try(regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist ,parallel = TRUE)$crit)
      if (class(lr.crit) == "try-error"){
        lr.crit <- c(0,0,0) 
      } 
      estimate.m.df[T, M] <- paste("$", round(lr.estimate, 2), "^{", paste(rep("*", sum(lr.estimate > lr.crit)), collapse = ""), "}", "$", sep = "")
      
      crit.m.df[T,M] <- paste(round(lr.crit,2),collapse = ",")
      # ############################################
      # #Estimate with regressor lnK
      # ############################################
      # data <- list(Y = t(ind.each.y), X = matrix(ind.each.t$lnk),  Z = NULL)
      # out.h0 <- regpanelmixPMLE(y=data$Y,x=data$X, z = NULL,m=M,vcov.method = "none")
      # # phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
      # #            beta = beta, N = N, T = T, M = M, p = p, q = q)
      # an <- anFormula(out.h0$parlist,M,N,T,q=1) 
      # print(paste("T=",T,"M = ",M,"an=",an))
      # if (is.na(an)){ 
      #   an <- 1.0
      # }
      # out.h1 <- regpanelmixMaxPhi(y=data$Y,x=data$X, z = NULL,parlist=out.h0$parlist,an=an)
      # lr.estimate <- 2 * max(out.h1$penloglik - out.h0$loglik)
      # estimate.k.df[T,M] <- lr.estimate
      # lr.crit <- try(regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist ,parallel = TRUE)$crit)
      # if (class(lr.crit) == "try-error"){
      #   lr.crit <- c(0,0,0) 
      # } 
      # crit.k.df[T,M] <- paste(round(lr.crit,2),collapse = ",")
      # ############################################
      # #Estimate with regressor lnL
      # ############################################
      # data <- list(Y = t(ind.each.y), X = matrix(ind.each.t$lnl),  Z = NULL)
      # out.h0 <- regpanelmixPMLE(y=data$Y,x=data$X, z = NULL,m=M,vcov.method = "none")
      # # phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
      # #            beta = beta, N = N, T = T, M = M, p = p, q = q)
      # an <- anFormula(out.h0$parlist,M,N,T,q=1) 
      # print(paste("T=",T,"M = ",M,"an=",an))
      # if (is.na(an)){ 
      #   an <- 1.0
      # }
      # out.h1 <- regpanelmixMaxPhi(y=data$Y,x=data$X, z = NULL,parlist=out.h0$parlist,an=an)
      # lr.estimate <- 2 * max(out.h1$penloglik - out.h0$loglik)
      # estimate.l.df[T,M] <- lr.estimate
      # lr.crit <- try(regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist ,parallel = TRUE)$crit)
      # if (class(lr.crit) == "try-error"){
      #   lr.crit <- c(0,0,0) 
      # } 
      # crit.l.df[T,M] <- paste(round(lr.crit,2),collapse = ",")
    }
    desc.each[T,] <- c(ind.name,each.code,dim(ind.each.y)[1],round(mean(as.matrix(ind.each.y)),2),round(sd(as.matrix(ind.each.y)),2))
  }
  ###################################################################
  #     Output
  ###################################################################

  
  estimate.LR.df.2[count,] <- estimate.m.df[2,]
  estimate.LR.df.3[count,] <- estimate.m.df[3,]
  estimate.LR.df.4[count,] <- estimate.m.df[4,]
  estimate.LR.df.5[count,] <- estimate.m.df[5,]
  
  colnames(desc.each) <- c("Industry Name","code","n","mean","sd") 
  rownames(desc.each) <- c("T=1","T=2","T=3","T=4","T=5")
  
  colnames(estimate.l.df) <- c("M=1","M=2","M=3","M=4","M=5")
  rownames(estimate.l.df) <- c("T=1","T=2","T=3","T=4","T=5")
  colnames(estimate.m.df) <- c("M=1","M=2","M=3","M=4","M=5")
  rownames(estimate.m.df) <- c("T=1","T=2","T=3","T=4","T=5")
  colnames(estimate.k.df) <- c("M=1","M=2","M=3","M=4","M=5")
  rownames(estimate.k.df) <- c("T=1","T=2","T=3","T=4","T=5")
  
  colnames(crit.l.df) <- c("M=1","M=2","M=3","M=4","M=5")
  rownames(crit.l.df) <- c("T=1","T=2","T=3","T=4","T=5")
  
  colnames(crit.m.df) <- c("M=1","M=2","M=3","M=4","M=5")
  rownames(crit.m.df) <- c("T=1","T=2","T=3","T=4","T=5")
  
  colnames(crit.k.df) <- c("M=1","M=2","M=3","M=4","M=5")
  rownames(crit.k.df) <- c("T=1","T=2","T=3","T=4","T=5")
  
  count = count + 1
  print(Sys.time()-t)
  print(paste(count,"/",ind.count))
  
  sink(paste("results/Chile/regressorCrit",ind.name,".txt"))
  
  stargazer(desc.each,type="text",title=paste("Descriptives for ",ind.name,each.code))
  
  print(paste("Chilean Producer Data: Estimated (lnM) for ",ind.name))
  print(estimate.m.df)
  # stargazer(estimate.m.df,type='latex',title = paste("Chilean Producer Data: Estimated (lnM) for ",ind.name))
  stargazer(crit.m.df,type="latex",title=paste("Simulated crit for(lnM)",ind.name,each.code))

  # stargazer(estimate.l.df,type='text',title = paste("Chilean Producer Data: Estimated LR(lnL) for ",ind.name))
  # stargazer(crit.l.df,type="text",title=paste("Simulated crit for(lnL)",ind.name,each.code))
  # 
  # stargazer(estimate.k.df,type='text',title = paste("Chilean Producer Data: Estimated LR(lnK) for ",ind.name))
  # stargazer(crit.k.df,type="text",title=paste("Simulated crit for(lnK)",ind.name,each.code))
  # stargazer(result.df,title = ind.name)
  
  sink()
}



# colnames(crit.df.boot) <- c("M=1","M=2","M=3","M=4","M=5")
# rownames(crit.df.boot) <- c("T=1","T=2","T=3","T=4","T=5")


# stargazer(crit.df,title=paste("estimate",ind.name))
rownames(estimate.LR.df.2) <- ind.names
rownames(estimate.LR.df.3) <- ind.names
rownames(estimate.LR.df.4) <- ind.names
rownames(estimate.LR.df.5) <- ind.names

sink("results/Chile/resultRegressors.txt")
stargazer(estimate.LR.df.2)
stargazer(estimate.LR.df.3)
stargazer(estimate.LR.df.4)
stargazer(estimate.LR.df.5)
sink()

