library(stargazer)
library(ggplot2)
library(reshape)
# library(normalregMix)
library(foreign)
library(NormalRegPanelMixture)
library(dplyr)

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

estimate.LR.df.3 <- matrix(0,nr=length(ind.code),nc=10)
rownames(estimate.LR.df.3) <- ind.names
colnames(estimate.LR.df.3) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")

AIC.df.3 <- matrix(0,nr=length(ind.code),nc=10)
rownames(AIC.df.3) <- ind.names
colnames(AIC.df.3) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")

BIC.df.3 <- matrix(0,nr=length(ind.code),nc=10)
rownames(BIC.df.3) <- ind.names
colnames(BIC.df.3) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")

Nonpar.df.3 <- matrix(0, nr = length(ind.code), nc = 10)
rownames(Nonpar.df.3) <- ind.names
colnames(Nonpar.df.3) <- c("M=1", "M=2", "M=3", "M=4", "M=5", "M=6", "M=7", "M=8", "M=9", "M=10")

crit.LR.df.3 <- matrix(0,nr=length(ind.code),nc=10)
rownames(crit.LR.df.3) <- ind.names
colnames(crit.LR.df.3) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")

count = 0
for (each.code in ind.code){
  t <- Sys.time()
  ind.each <- df %>%
    filter(ciiu_3d == each.code) %>%
    mutate(
      y = log(GO),
      lnm = log(WI),
      lnl = log(L),
      lnk = log(K)
    ) %>%
    group_by(id) %>%
    arrange(id, year) %>%
    mutate(
      y_l1 = lag(si, n = 1, default = NA),
      lnk_l1 = lag(lnk, n = 1, default = NA),
      lnl_l1 = lag(lnl, n = 1, default = NA),
      lnm_l1 = lag(lnm, n = 1, default = NA)
    ) %>%
    ungroup()
  
  ind.name <- ind.each$ciiu3d_descr[1]
  
  ######################################################
  # Describe the data
  ######################################################
  
  desc.each <- ind.each[ind.each$L != 0 ,c("si","y","lnm","lnl","lnk", "y_l1","lnm_l1","lnl_l1","lnk_l1")]
  # desc.each <- desc.each[complete.cases(desc.each),]
  year.list <- sort(unique(ind.each$year))
  T.cap <- max(year.list)
  
  coef.df <- matrix(0,nr=5,nc=10)
  estimate.df <- matrix(0,nr=5,nc=10)
  AIC.df <- matrix(0,nr=5,nc=10)
  BIC.df <- matrix(0, nr = 5, nc = 10)
  Nonpar.df <- matrix(0, nr = 5, nc = 10)
  crit.df <- matrix(0, nr = 5, nc = 10)
  result.df <- matrix(0,nr=5,nc=10)
  ######################################################
  #For panel data
  ######################################################
  
  
  for (T in 3:3){
    t.start <- T.cap - T + 1
    # Reshape the data so that I can apply the test
    ind.each.t <- ind.each[ind.each$year >= t.start, ]
    ind.each.t <- ind.each.t[complete.cases(ind.each.t), ]
    ind.each.y <- cast(ind.each.t[, c("id", "year", "si")], id ~ year, value = "si")
    id.list <- ind.each.y[complete.cases(ind.each.y), "id"]
    # Remove the incomplete data, need balanced panel
    ind.each.t <- ind.each.t[ind.each.t$id %in% id.list, ]
    ind.each.t <- ind.each.t[order(ind.each.t$id, ind.each.t$year), ]
    
    #Order the range of IDs
    ind.each.t <- ind.each.t %>%
      arrange(id)
    
    # Reshape the Y
    ind.each.y <- cast(ind.each.t[, c("id", "year", "y")], id ~ year, value = "y")
    ind.each.y <- ind.each.y[, colnames(ind.each.y) != "id"]
    
    data <- list(
      Y = t(ind.each.y), 
      X = data.frame(ind.each.t[,c("y_l1")]), 
      Z = NULL)
    
    data.0 <- list(
      Y = ind.each.t[ind.each.t$year==t.start,]$y_l1,
      X = NULL, 
      Z = NULL) # for the initial condition
    
    N <- dim(ind.each.y)[1]
    
    h1.coefficient = NULL
    estimate.crit <- 1
    
    for (M in 1:5){
      # Estimate the null model
      out.h0 <- regpanelmixPMLE(y=data$Y, x=data$X, z = data$Z, m=M,vcov.method = "none", in.coefficient=h1.coefficient, data.0= data.0)
      
      an <- anFormula(out.h0$parlist,M,N,T,q=1)
      an_0 <- anFormula.t0(out.h0$parlist0,M,N,q=1)
      print("-----------------------------------------")
      print(paste("T =",T,"M = ",M,"an =",an, "an_0 =",an_0))
      if (is.na(an)){
        an <- 1.0
      }
      # Estimate the alternative model
      out.h1 <- regpanelmixMaxPhi(y=data$Y,x=data$X, z = data$Z,parlist=out.h0$parlist,an=an, an_0=an_0, data.0 = data.0)
      h1.parlist = out.h1$parlist
      
      lr.estimate <- 2 * max(out.h1$penloglik - out.h0$loglik)
      
      # Simulate the asymptotic distribution, for AR1 only use Bootstrapped
      if (estimate.crit){
        lr.crit <- regpanelmixCritBootAR1(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, cl=cl,parallel = TRUE, data.0=data.0, an = an, an_0 = an_0, ninits=1)$crit
        estimate.df[T,M] <- paste('$',round(lr.estimate,2),'^{',paste(rep('*',sum(lr.estimate > lr.crit)),  collapse = ""),'}','$', sep = "")
      } else{
        lr.crit <- c(0,0,0)
        estimate.df[T,M] <- paste('$',round(lr.estimate,2),'$', sep = "")
      }
      
      
      coef.df[T,M] <- paste(paste(names(out.h0$coefficients), collapse = ","), paste(out.h0$coefficients, collapse = ","), collapse = ",")
      
      AIC.df[T,M] <- round(out.h0$aic,2)
      BIC.df[T,M] <- round(out.h0$bic,2)
      
      
      crit.df[T,M] <- paste(round(lr.crit,2),collapse = ",")
      # If fail to reject the test, break the loop
      print(lr.estimate)
      print(lr.crit)
      
      if (sum(lr.estimate > lr.crit) < 2){
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
  estimate.LR.df.3[count,] <- estimate.df[3,]
  AIC.df.3[count,] <- AIC.df[3,]
  BIC.df.3[count, ] <- BIC.df[3, ]
  Nonpar.df.3[count, ] <- Nonpar.df[3, ]
  crit.LR.df.3[count,] <- crit.df[3,]
  colnames(estimate.df) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
  rownames(estimate.df) <- c("T=1","T=2","T=3","T=4","T=5")
  
  colnames(crit.df) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
  rownames(crit.df) <- c("T=1","T=2","T=3","T=4","T=5")
  
  #sink(paste("/home/haoyu/results/Empirical/Chile_crit",ind.name,"_regressor.txt"))
  sink(paste("results/Empirical/Chile_Crit_", ind.name, "_AR1.txt"))
  stargazer(as.data.frame(desc.each),type="text",summary=TRUE,title=paste("Descriptive data for Chilean Industry: ",ind.name))
  print(paste("Chilean Producer Data: Estimated LR for",ind.name))
  print(coef.df)
  print(estimate.df)
  stargazer(crit.df,type="text",title=paste("Simulated crit for ",ind.name,each.code))
  sink()
}



# colnames(crit.df.boot) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
# rownames(crit.df.boot) <- c("T=1","T=2","T=3","T=4","T=5")

count <- length(ind.names)
df.3 <- data.frame(matrix('-',nrow=3*length(ind.names),ncol=10))
df.3[ 3* 1:count -2,] <- estimate.LR.df.3
df.3[ 3* 1:count -1,] <- AIC.df.3
df.3[ 3* 1:count ,] <- BIC.df.3
rownames(df.3)[ 3* 1:count -2] <- rownames(estimate.LR.df.3)
colnames(df.3) <- colnames(estimate.LR.df.3)


# Combine the data frames
combined_df <- rbind(
  cbind(df.3, original_df = "df.3")
)


# Write the combined data frame to a single file
write.csv(combined_df, file = "results/Empirical/Chile_AR1.csv")


