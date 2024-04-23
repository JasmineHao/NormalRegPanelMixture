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
cl <- makeCluster(7)
count = 0

estimate.LR.df.2 <- matrix(0,nr=length(ind.code),nc=10)
rownames(estimate.LR.df.2) <- ind.names
colnames(estimate.LR.df.2) <- c("M=1","M=2","M=3","M=4","M=5","M=6","M=7","M=8","M=9","M=10")
estimate.LR.df.3 <- matrix(0,nr=length(ind.code),nc=10)
rownames(estimate.LR.df.3) <- ind.names
colnames(estimate.LR.df.3) <- c("M=1","M=2","M=3","M=4","M=5","M=6","M=7","M=8","M=9","M=10")
estimate.LR.df.4 <- matrix(0,nr=length(ind.code),nc=10)
rownames(estimate.LR.df.4) <- ind.names
colnames(estimate.LR.df.4) <- c("M=1","M=2","M=3","M=4","M=5","M=6","M=7","M=8","M=9","M=10")
estimate.LR.df.5 <- matrix(0,nr=length(ind.code),nc=10)
rownames(estimate.LR.df.5) <- ind.names
colnames(estimate.LR.df.5) <- c("M=1","M=2","M=3","M=4","M=5","M=6","M=7","M=8","M=9","M=10")
# for (each.code in ind.code){
#   t <- Sys.time()
#   ind.each <- subset(df,ciiu_3d==each.code)
#   ind.name <- ind.each$ciiu3d_descr[1]
#   ind.each$y <- log(ind.each$GO)
#   ind.each$lnm <- log(ind.each$WI)
#   ind.each$lnl <- log(ind.each$L)
#   ind.each$lnk <- log(ind.each$K)
#   ######################################################
#   # Describe the data
#   ######################################################

#   desc.each <- ind.each[ind.each$L != 0 ,c("si","y","lnm","lnl","lnk")]
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
  desc.each <- desc.each[complete.cases(desc.each),]
  
  ######################################################
  # Select the data out
  ######################################################
  m.share <- cast(ind.each,id ~ year,value="si")#Collapse the dataframe into panel form , year against firm id
  row.names(m.share) <- m.share$id 
  m.share <- m.share[,!(colnames(m.share)=="id")] 
  T.cap <- dim(m.share)[2]
  
  estimate.df <- matrix(0,nr=5,nc=10)
  crit.df <- matrix(0,nr=5,nc=10)
  result.df <- matrix(0,nr=5,nc=10)
  ######################################################
  #For panel data
  ######################################################
  
  
  for (T in 3:3){
    t.start <- T.cap-T+1
    t.seq <- seq(from=t.start,to=t.start+T-1)
    m.share.t <- m.share[,t.seq]
    data <- list(Y = t(m.share.t[complete.cases(m.share.t),]), X = NULL,  Z = NULL)
    N <- dim(data$Y)[2]
    
    for (M in 1:10){
      tryCatch({
      out.h0 <- normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none")
      # phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
      #            beta = beta, N = N, T = T, M = M, p = p, q = q)
      an <- anFormula(out.h0$parlist,M,N,T) 
      print(paste("T=",T,"M = ",M,"an=",an))
      
      if (is.na(an)){ 
        an <- 1.0
      }
      out.h1 <- normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an=an)
      lr.estimate <- 2 * max(out.h1$penloglik - out.h0$loglik)
      
      
      lr.crit <- try(regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, cl=cl,parallel = TRUE)$crit)
      
      if (class(lr.crit) == "try-error"){
        
        lr.crit <- c(0,0,0) 
        
      }
      estimate.df[T,M] <- paste('$',round(lr.estimate,2),'^{',paste(rep('*',sum(lr.estimate > lr.crit)),  collapse = ""),'}','$', sep = "")
      crit.df[T,M] <- paste(round(lr.crit,2),collapse = ",")
    }
      ,
      error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      #lr.crit.boot <- regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist,cl=cl, nbtsp = 199 ,parallel = TRUE)$crit
      # crit.df.boot[T,M] <- paste(round(lr.crit.boot,2),collapse = ",")
    }
  }
  ###################################################################
  #     Output
  ###################################################################
  count = count + 1
  print(Sys.time()-t)
  print(paste(count,"/",ind.count))
  
  estimate.LR.df.2[count,] <- estimate.df[2,]
  estimate.LR.df.3[count,] <- estimate.df[3,]
  estimate.LR.df.4[count,] <- estimate.df[4,]
  estimate.LR.df.5[count,] <- estimate.df[5,]
  
  
  colnames(estimate.df) <- c("M=1","M=2","M=3","M=4","M=5","M=6","M=7","M=8","M=9","M=10")
  rownames(estimate.df) <- c("T=1","T=2","T=3","T=4","T=5")
  
  colnames(crit.df) <- c("M=1","M=2","M=3","M=4","M=5","M=6","M=7","M=8","M=9","M=10")
  rownames(crit.df) <- c("T=1","T=2","T=3","T=4","T=5")
  
  sink(paste("results/Empirical/Chile_crit_1-10",ind.name,".txt"))
  
  stargazer(desc.each,type="text",title=paste("Descriptive data for Chilean Industry: ",ind.name))
  print(paste("Columbian Producer Data: Estimated LR for ",ind.name))
  stargazer(estimate.df)
  stargazer(crit.df,type="text",title=paste("Simulated crit for ",ind.name,each.code))
  # stargazer(crit.df,type="text",title=paste("Bootstrapped crit for ",ind.name,each.code))
  sink()
}


# stargazer(crit.df,title=paste("estimate",ind.name))
rownames(estimate.LR.df.2) <- ind.names
rownames(estimate.LR.df.3) <- ind.names
rownames(estimate.LR.df.4) <- ind.names
rownames(estimate.LR.df.5) <- ind.names

# colnames(crit.df.boot) <- c("M=1","M=2","M=3","M=4","M=5")
# rownames(crit.df.boot) <- c("T=1","T=2","T=3","T=4","T=5")


# stargazer(crit.df,title=paste("estimate",ind.name))

sink("results/Empirical/Chile_resultBoot_1-10.txt")
stargazer(estimate.LR.df.2)
stargazer(estimate.LR.df.3)
stargazer(estimate.LR.df.4)
stargazer(estimate.LR.df.5)
sink()



