# This script replicates table 7, 10, 11 - 14 in the paper.
library(stargazer)
library(ggplot2)
library(reshape)
library(NormalRegPanelMixture)
# library(Hmisc)
set.seed(123)
library(dplyr)

# library(normalregMix)

options('nloptr.show.inequality.warning'=FALSE)
cl <- makeCluster(12)
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
df <- readRDS("data/JapanClean.rds")
df[df==0] <- NA
#Function
#source("C:/Users/Jasmine/Dropbox/GNR/R/productionEstimation.R")
df$t <- df$year
df <- df[order(df$id,df$t),]



ind.code <- c(5,13,12,14,1)
ind.code <- c(5,13,12)
ind.names <- c()
for (each.code in ind.code){
  ind.name <- ind_list[each.code]
  ind.names <- append(ind.names,ind.name)
}

ind.count <- length(ind.code)


estimate.LR.df.2 <- matrix(0,nr=length(ind.code),nc=10)
rownames(estimate.LR.df.2) <- ind.names
colnames(estimate.LR.df.2) <- c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
estimate.LR.df.3 <- matrix(0,nr=length(ind.code),nc=10)
rownames(estimate.LR.df.3) <- ind.names
colnames(estimate.LR.df.3) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
estimate.LR.df.4 <- matrix(0,nr=length(ind.code),nc=10)
rownames(estimate.LR.df.4) <- ind.names
colnames(estimate.LR.df.4) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
estimate.LR.df.5 <- matrix(0,nr=length(ind.code),nc=10)
rownames(estimate.LR.df.5) <- ind.names
colnames(estimate.LR.df.5) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")

AIC.df.2 <- matrix(0,nr=length(ind.code),nc=10)
rownames(AIC.df.2) <- ind.names
colnames(AIC.df.2) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
AIC.df.3 <- matrix(0,nr=length(ind.code),nc=10)
rownames(AIC.df.3) <- ind.names
colnames(AIC.df.3) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
AIC.df.4 <- matrix(0,nr=length(ind.code),nc=10)
rownames(AIC.df.4) <- ind.names
colnames(AIC.df.4) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
AIC.df.5 <- matrix(0,nr=length(ind.code),nc=10)
rownames(AIC.df.5) <- ind.names
colnames(AIC.df.5) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")


BIC.df.2 <- matrix(0,nr=length(ind.code),nc=10)
rownames(BIC.df.2) <- ind.names
colnames(BIC.df.2) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
BIC.df.3 <- matrix(0,nr=length(ind.code),nc=10)
rownames(BIC.df.3) <- ind.names
colnames(BIC.df.3) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
BIC.df.4 <- matrix(0,nr=length(ind.code),nc=10)
rownames(BIC.df.4) <- ind.names
colnames(BIC.df.4) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
BIC.df.5 <- matrix(0,nr=length(ind.code),nc=10)
rownames(BIC.df.5) <- ind.names
colnames(BIC.df.5) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")



crit.LR.df.2 <- matrix(0,nr=length(ind.code),nc=10)
rownames(crit.LR.df.2) <- ind.names
colnames(crit.LR.df.2) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
crit.LR.df.3 <- matrix(0,nr=length(ind.code),nc=10)
rownames(crit.LR.df.3) <- ind.names
colnames(crit.LR.df.3) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
crit.LR.df.4 <- matrix(0,nr=length(ind.code),nc=10)
rownames(crit.LR.df.4) <- ind.names
colnames(crit.LR.df.4) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
crit.LR.df.5 <- matrix(0,nr=length(ind.code),nc=10)
rownames(crit.LR.df.5) <- ind.names
colnames(crit.LR.df.5) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")


######################################################
#For panel data
######################################################

count = 0

for (each.code in ind.code){
  t <- Sys.time()
  ind.each <- subset(df,industry_2==each.code)
  ind.each <- ind.each[,c("id","year","lnmY_it","k_it")]
  ind.each <- ind.each[complete.cases(ind.each),]
  ind.each['lnk'] <- ind.each['k_it']
  ind.name <- ind_list[each.code]
  
  ind.each <- ind.each %>%
    group_by(id) %>%
    arrange(id, year) %>%
    mutate(
      y_l1 = lag(lnmY_it, n = 1, default = NA),
      lnk_l1 = lag(lnk, n = 1, default = NA),
    ) %>%
    ungroup()
  
  
  year.list <- sort(unique(ind.each$year))
  T.cap <- max(year.list) 
  
  coef.df <- matrix(0,nr=5,nc=10)
  estimate.df <- matrix(0,nr=5,nc=10)
  AIC.df <- matrix(0,nr=5,nc=10)
  BIC.df <- matrix(0,nr=5,nc=10)
  crit.df <- matrix(0,nr=5,nc=10)
  result.df <- matrix(0,nr=5,nc=10)
  ######################################################
  #For panel data
  ######################################################
  
  
  for (T in 3:3){
    t.start <- T.cap-T+1
    t.seq <- seq(from=t.start,to=t.start+T-1)
    
    ind.each.t <- ind.each[ind.each$year >= t.start,]
    ind.each.y <- cast(ind.each.t[,c("id","year","lnmY_it")],id ~ year,value="lnmY_it")
    id.list    <- ind.each.y[complete.cases(ind.each.y),"id"]
    #Remove the incomplete data, need balanced panel
    ind.each.t <- ind.each.t[ind.each.t$id %in% id.list,]
    ind.each.t <- ind.each.t[order(ind.each.t$id,ind.each.t$year),]
    
    # normalize data
    ind.each.t$y <- (ind.each.t$lnmY_it - mean(ind.each.t$lnmY_it)) / (sd(ind.each.t$lnmY_it))
    ind.each.t$x <- (ind.each.t$lnk - mean(ind.each.t$lnk))/(sd(ind.each.t$lnk))
    
    ind.each.t$y1 <- (ind.each.t$y_l1 - mean(ind.each.t$y_l1)) / (sd(ind.each.t$y_l1))
    ind.each.t$x1 <- (ind.each.t$lnk_l1 - mean(ind.each.t$lnk_l1))/(sd(ind.each.t$lnk_l1))
    

    #Reshape the Y 
    ind.each.y <- cast(ind.each.t[, c("id", "year", "y")], id ~ year, value = "y")
    ind.each.y <- ind.each.y[, colnames(ind.each.y) != "id"]

    ind.each.x <- ind.each.t$x

    ind.each.y1 <- ind.each.t$y1

    ind.each.y10 <- cast(ind.each.t[, c("id", "year", "y1")], id ~ year, value = "y1")
    ind.each.y10 <- ind.each.y10[, colnames(ind.each.y10) != "id"][, 1]

    ind.each.x1 <- cast(ind.each.t[, c("id", "year", "x1")], id ~ year, value = "x1")
    ind.each.x1 <- ind.each.x1[, colnames(ind.each.x1) != "id"][, 1]

    data <- list(Y = t(ind.each.y), X = data.frame(col1 = ind.each.y1), Z = NULL)
    data.0 <- list(Y = ind.each.y10, X = NULL, Z = NULL) # for the initial
    
    N <- dim(data$Y)[2]
    
    h1.coefficient = NULL
    estimate.crit <- 1
    for (M in 1:10){
      # Estimate the null model
      out.h0 <- regpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none",in.coefficient=h1.coefficient, data.0 = data.0)
      an <- anFormula(out.h0$parlist,M,N,T,q=1)
      print("-----------------------------------------")
      print(paste("T=",T,"M = ",M,"an=",an))
      if (is.na(an)){
        an <- 1.0
      }
      # Estimate the alternative model
      out.h1 <- regpanelmixMaxPhi(y=data$Y,x=data$X, z = data$Z,parlist=out.h0$parlist,an=an, data.0 = data.0)
      h1.parlist = out.h1$parlist
      
      lr.estimate <- 2 * max(out.h1$penloglik - out.h0$loglik)
      
      # Simulate the asymptotic distribution, for AR1 only use Bootstrapped
      if (estimate.crit){
        lr.crit <- regpanelmixCritBootAR1(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, cl=cl,parallel = TRUE, data.0=data.0)$crit
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
  
  # sink(paste("/home/haoyu/results/Empirical/Japan_Crit",ind.name,"_regressor.txt"))
  sink(paste("results/Empirical/Japan_Crit_", ind.name, "_AR1.txt"))
  stargazer(ind.each,type="latex",title=paste("Descriptive data for ",ind.name, " industry in Japan"))
  print(paste("Estimate LR for ",ind.name))
  print(coef.df)
  print(estimate.df)
  stargazer(crit.df,title=paste("Critical Values Asymptotics",each.code))
  # regpanelmixMEMtest(y = data$Y,x=NULL,t=5,m=2,crit.method="none")
  #stargazer(crit.df.boot,title=paste("Critical Values Bootstrapped",each.code))
  sink()
}

# write.csv(cbind(estimate.LR.df.3,AIC.df.3),file="/home/haoyu/results/Empirical/Japan_resultLR3_regressor.csv")
#     write.csv(cbind(estimate.LR.df.4,AIC.df.4),file="/home/haoyu/results/Empirical/Japan_resultLR4_regressor.csv")
# write.csv(cbind(estimate.LR.df.5,AIC.df.5),file="/home/haoyu/results/Empirical/Japan_resultLR5_regressor.csv")

count <- length(ind.names)
df.2 <- data.frame(matrix('-',nrow=3*length(ind.names),ncol=10))
df.2[ 3* 1:count -2,] <- estimate.LR.df.2
df.2[ 3* 1:count -1,] <- AIC.df.2
df.2[ 3* 1:count ,] <- BIC.df.2
rownames(df.2)[ 3* 1:count -2] <- rownames(estimate.LR.df.2)
colnames(df.2) <- colnames(estimate.LR.df.2)


df.3 <- data.frame(matrix('-',nrow=3*length(ind.names),ncol=10))
df.3[ 3* 1:count -2,] <- estimate.LR.df.3
df.3[ 3* 1:count -1,] <- AIC.df.3
df.3[ 3* 1:count ,] <- BIC.df.3
rownames(df.3)[ 3* 1:count -2] <- rownames(estimate.LR.df.3)
colnames(df.3) <- colnames(estimate.LR.df.3)

df.4 <- data.frame(matrix('-',nrow=3*length(ind.names),ncol=10))
df.4[ 3* 1:count -2,] <- estimate.LR.df.4
df.4[ 3* 1:count -1,] <- AIC.df.4
df.4[ 3* 1:count ,] <- BIC.df.4
rownames(df.4)[ 3* 1:count -2] <- rownames(estimate.LR.df.4)
colnames(df.4) <- colnames(estimate.LR.df.4)


df.5 <- data.frame(matrix('-',nrow=3*length(ind.names),ncol=10))
df.5[ 3* 1:count -2,] <- estimate.LR.df.5
df.5[ 3* 1:count -1,] <- AIC.df.5
df.5[ 3* 1:count ,] <- BIC.df.5
rownames(df.5)[ 3* 1:count -2] <- rownames(estimate.LR.df.5)
colnames(df.5) <- colnames(estimate.LR.df.5)


# Combine the data frames
combined_df <- rbind(
  cbind(df.2, original_df = "df.2"),
  cbind(df.3, original_df = "df.3"),
  cbind(df.4, original_df = "df.4"),
  cbind(df.5, original_df = "df.5")
)

# Write the combined data frame to a single file
write.csv(combined_df, file = "results/Empirical/Japan_AR1.csv")


