# This script replicates table 7, 10, 11 - 14 in the paper.
library(stargazer)
library(ggplot2)
library(reshape)
library(NormalRegPanelMixture)
# library(Hmisc)
set.seed(123)
library(dplyr)
library(splines)
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

estimate.LR.df.3 <- matrix(0,nr=length(ind.code),nc=10)
rownames(estimate.LR.df.3) <- ind.names
colnames(estimate.LR.df.3) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")

AIC.df.3 <- matrix(0,nr=length(ind.code),nc=10)
rownames(AIC.df.3) <- ind.names
colnames(AIC.df.3) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")

BIC.df.3 <- matrix(0,nr=length(ind.code),nc=10)
rownames(BIC.df.3) <- ind.names
colnames(BIC.df.3) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")

crit.LR.df.3 <- matrix(0,nr=length(ind.code),nc=10)
rownames(crit.LR.df.3) <- ind.names
colnames(crit.LR.df.3) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")

######################################################
#For panel data
######################################################

count = 0

bs_degree <- 2

for (each.code in ind.code){
  t <- Sys.time()
  ind.each <- subset(df,industry_2==each.code)
  ind.each <- ind.each[,c("id","year","lnmY_it","k_it","l_it")]
  ind.each <- ind.each[complete.cases(ind.each),]
  ind.each['lnk'] <- (ind.each$k_it - mean(ind.each$k_it) )/(sd(ind.each$k_it))
  ind.each['lnl'] <- (ind.each$l_it - mean(ind.each$l_it) )/(sd(ind.each$l_it))
  ind.each['y'] <- (ind.each$lnmY_it - mean(ind.each$lnmY_it) )/(sd(ind.each$lnmY_it))
  
  ind.name <- ind_list[each.code]
  # first create splines
  
  ind.each <- ind.each %>%
    group_by(id) %>%
    arrange(id, year) %>%
    mutate(
      y_l1 = lag(y, n = 1, default = NA),
      lnk_l1 = lag(lnk, n = 1, default = NA),
      lnl_l1 = lag(lnl, n = 1, default = NA),
    ) %>%
    ungroup()
  
  
  # List of column names to apply the bs transformation
  # columns_to_transform <- c("lnk", "lnk_l1", "lnl", "lnl_l1" , "y_l1")
  
  columns_to_transform <- c("lnk", "lnl")
  
  # Loop through each column name
  for (col_name in columns_to_transform) {
    # Compute quantiles excluding NA values
    cleaned_data <- na.omit(ind.each[[col_name]])
    quantiles <- c(quantile(cleaned_data, probs = 0.33), quantile(cleaned_data, probs = 0.67))
    
    # Apply the bs function and create new columns
    bs_columns <- bs(ind.each[[col_name]], degree = bs_degree, knots = quantiles)
    new_colnames <- paste0(col_name, "_bs_", seq_len(ncol(bs_columns)))
    colnames(bs_columns) <- new_colnames
    
    # Bind the new columns to the original dataframe
    ind.each <- cbind(ind.each, bs_columns)
  }
  
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
    ind.each.y <- cast(ind.each.t[,c("id","year","y")],id ~ year,value="y")
    id.list    <- ind.each.y[complete.cases(ind.each.y),"id"]
    #Remove the incomplete data, need balanced panel
    ind.each.t <- ind.each.t[ind.each.t$id %in% id.list,]
    ind.each.t <- ind.each.t[order(ind.each.t$id,ind.each.t$year),]
    #Reshape the Y 
    ind.each.y <- cast(ind.each.t[,c("id","year","y")],id ~ year,value="y")
    ind.each.y <- ind.each.y[,colnames(ind.each.y)!="id"]
    
    # Find all columns in ind.each that contain '_bs_'
    bs_columns <- grep("_bs_", colnames(ind.each.t), value = TRUE)
    
    # Create a data frame with all the bs_columns
    X_bs <- ind.each.t[, bs_columns]
    
    # Create the list with Y, X (with bs columns), and Z
    data <- list(
      Y = t(ind.each.y),
      X = X_bs,
      Z = NULL
    )
    
    
    N <- dim(data$Y)[2]
    
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
      out.h1 <- regpanelmixMaxPhi(y=data$Y,x=data$X, z = data$Z,parlist=out.h0$parlist,an=an)
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
      
      AIC.df[T,M] <- round(out.h0$aic,2)
      BIC.df[T,M] <- round(out.h0$bic,2)
      crit.df[T,M] <- paste(round(lr.crit,2),collapse = ",")
      if (estimate.crit == 1){
        estimate.df[T,M] <- paste('$',round(lr.estimate,2),'^{',paste(rep('*',sum(lr.estimate > lr.crit)),  collapse = ""),'}','$', sep = "")
        
      }
      else{
        estimate.df[T,M] <- paste('$',round(lr.estimate,2),'$', sep = "")
        
      }
      # If fail to reject the test, break the loop
      print(lr.estimate)
      print(lr.crit)
      
      if (sum(lr.estimate > lr.crit) < 1){
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
  BIC.df.3[count,] <- BIC.df[3,]
  crit.LR.df.3[count,] <- crit.df[3,]
  colnames(estimate.df) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
  rownames(estimate.df) <- c("T=1","T=2","T=3","T=4","T=5")
  
  colnames(crit.df) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
  rownames(crit.df) <- c("T=1","T=2","T=3","T=4","T=5")
  
  # sink(paste("/home/haoyu/results/Empirical/Japan_Crit",ind.name,"_regressor.txt"))
  sink(paste("results/Empirical/Japan_Crit_", ind.name, "_KL_Spline"))
  
  stargazer(ind.each,type="latex",title=paste("Descriptive data for ",ind.name, " industry in Japan"))
  print(paste("Estimate LR for ",ind.name))
  print(coef.df)
  print(estimate.df)
  stargazer(crit.df,title=paste("Critical Values Asymptotics",each.code))
  
  sink()
}

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
write.csv(combined_df, file = "results/Empirical/Japan_KL_Spline.csv")


