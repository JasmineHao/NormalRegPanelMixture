library(stargazer)
library(reshape)
# library(normalregMix)
library(foreign)
library(NormalRegPanelMixture)
library(haven)
library(dplyr)
library(splines)

# library(Hmisc)
set.seed(123)

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



estimate.LR.df.3 <- matrix(0,nr=length(ind.code),nc=10)
rownames(estimate.LR.df.3) <- ind.names
colnames(estimate.LR.df.3) <- c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")

AIC.df.3 <- matrix(0,nr=length(ind.code),nc=10)
rownames(AIC.df.3) <- ind.names
colnames(AIC.df.3) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")

BIC.df.3 <- matrix(0,nr=length(ind.code),nc=10)
rownames(BIC.df.3) <- ind.names
colnames(BIC.df.3) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")

count <- 0

bs_degree <- 2

model <- "k_dummy"

# candidate models
# c("k_dummy")

for (each.code in ind.code){
  t <- Sys.time()
  print(each.code)
  ind.name <- ind_list[each.code]
  
  ind.each <- subset(df,industry_2==each.code)
  
  ind.each$y <- ind.each$lnmY_it
  ind.each$lnl <- ind.each$l_it
  ind.each$lnk <- ind.each$k_it
  # ind.each$export <- ind.each$export # dummy for export 
  # ind.each$import <- ind.each$import # dummy for import
  
  # List of column names to apply the bs transformation
  # columns_to_transform <- c("lnk", "lnk_l1", "lnl", "lnl_l1" , "y_l1") 
  # columns_to_transform <- c("lnk", "lnl")
  # 
  # # Loop through each column name
  # for (col_name in columns_to_transform) {
  #   # Compute quantiles excluding NA values
  #   cleaned_data <- na.omit(ind.each[[col_name]])
  #   quantiles <- c(quantile(cleaned_data, probs = 0.33), quantile(cleaned_data, probs = 0.67))
  #   
  #   # Apply the bs function and create new columns
  #   bs_columns <- bs(ind.each[[col_name]], degree = bs_degree, knots = quantiles)
  #   new_colnames <- paste0(col_name, "_bs_", seq_len(ncol(bs_columns)))
  #   colnames(bs_columns) <- new_colnames
  #   
  #   # Bind the new columns to the original dataframe
  #   ind.each <- cbind(ind.each, bs_columns)
  # }
  # 
  # 
  coef.df <- matrix(0,nr=5,nc=10)
  estimate.df <- matrix(0,nr=5,nc=10)
  AIC.df <- matrix(0,nr=5,nc=10)
  BIC.df <- matrix(0, nr = 5, nc = 10)
  Nonpar.df <- matrix(0, nr = 5, nc = 10)
  crit.df <- matrix(0, nr = 5, nc = 10)
  result.df <- matrix(0,nr=5,nc=10)
  
  T <- 3
  year.list <- sort(unique(ind.each$year))
  T.cap <- max(year.list)
  
  t.start <- T.cap-T+1
  #Reshape the data so that I can apply the test
  ind.each.t <- ind.each[ind.each$year>=t.start,]
  
  ind.each.y <- cast(ind.each.t[,c("id","year","lnmY_it")],id ~ year,value="lnmY_it")
  id.list    <- ind.each.y[complete.cases(ind.each.y),"id"]
  #Remove the incomplete data, need balanced panel
  ind.each.t <- ind.each.t[ind.each.t$id %in% id.list,]
  ind.each.t <- ind.each.t[order(ind.each.t$id,ind.each.t$year),]
  #Reshape the Y 
  ind.each.y <- cast(ind.each.t[,c("id","year","lnmY_it")],id ~ year,value="lnmY_it")
  ind.each.y <- ind.each.y[,colnames(ind.each.y)!="id"]
  ind.each.y <- (ind.each.y - mean(ind.each.t$y))/(sd(ind.each.t$y))
  
  ind.each.x <- (ind.each.t$lnk - mean(ind.each.t$lnk))/(sd(ind.each.t$lnk))
  
  #ind.each.ciiu_dummy <- model.matrix(~ ciiu - 1, data = ind.each.t)[,-1]
  #ind.each.export <- ind.each.t$export
  #ind.each.import <- ind.each.t$import
  ind.each.k_hl   <- ind.each.x > median(ind.each.x) # indicator for high capital v.s. low capital
  
  
  
  # data <- list(Y = t(ind.each.y), X = cbind(ind.each.x)  ,  Z = ind.each.ciiu_dummy)
  # Find all columns in ind.each that contain '_bs_'
  # bs_columns <- grep("_bs_", colnames(ind.each.t), value = TRUE)
  
  data <- list(Y = t(ind.each.y), X = cbind(ind.each.k_hl*1)  ,  Z = NULL)
  N <- dim(ind.each.y)[1]
  
  h1.coefficient = NULL
  
  estimate.crit <- 1
  for (M in 1:10){
    # Estimate the null model
    out.h0 <- regpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none",in.coefficient=h1.coefficient)
    an <- anFormula(out.h0$parlist,M,N,T,q=1)
    print("-----------------------------------------")
    print(paste("T=",T,"M = ",M,"an=",an))
    if (is.na(an)){
      an <- 1.0
    }
    # Estimate the alternative model
    out.h1 <- regpanelmixMaxPhi(y=data$Y,x=data$X, z = data$Z,parlist=out.h0$parlist,an=an,ninits = 10)
    h1.parlist = out.h1$parlist
    
    lr.estimate <- 2 * max(out.h1$penloglik - out.h0$loglik)
    
    # Simulate the asymptotic distribution
    if (estimate.crit == 1){
      lr.crit <- try(regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, cl=cl,parallel = TRUE)$crit)
      if (class(lr.crit) == "try-error"){
        lr.crit <- regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, cl=cl,parallel = TRUE)$crit
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
    
    if (sum(lr.estimate > lr.crit) < 2){
      estimate.crit <- 0
    }
  }
  count = count + 1
  print("*************************************")
  print(paste("Finished", ind.name))
  print( Sys.time() - t)
  print("*************************************")
  estimate.LR.df.3[count,] <- estimate.df[3,]
  AIC.df.3[count,] <- AIC.df[3,]
  BIC.df.3[count, ] <- BIC.df[3, ]
  
  sink(paste("results/Empirical/Additional/Japan_Crit_", ind.name,"_", model ,".txt",sep=""))
  print(paste("Japan Producer Data: Estimated LR for", ind.name))
  print(coef.df)
  print(estimate.df)
  stargazer(crit.df, type = "text", title = paste("Simulated crit for ", ind.name, each.code))
  sink()
  
}


df.3 <- data.frame(matrix('-',nrow=3*length(ind.names),ncol=10))
df.3[ 3* 1:count -2,] <- estimate.LR.df.3
df.3[ 3* 1:count -1,] <- AIC.df.3
df.3[ 3* 1:count ,] <- BIC.df.3
rownames(df.3)[ 3* 1:count -2] <- rownames(estimate.LR.df.3)
colnames(df.3) <- colnames(estimate.LR.df.3)



write.csv(df.3,file=paste("results/Empirical/Additional/Japan_regressior_", model,".csv",sep=""))
