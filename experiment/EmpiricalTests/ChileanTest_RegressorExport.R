library(stargazer)
library(ggplot2)
library(reshape)
# library(normalregMix)
library(foreign)
library(NormalRegPanelMixture)
library(haven)
library(dplyr)

cl <- makeCluster(12)

options('nloptr.show.inequality.warning'=FALSE)
options(warn = -1)

set.seed(123)

df <- readRDS("data/ChileanClean.rds")

# Merge data with export 
df.0 <- read_dta("data/DATA_KL_EXIM5_2023.dta")
df.0$id = df.0$padron

# Check if 'padron' column exists and is not empty
if (!'padron' %in% colnames(df.0)) {
  stop("Column 'padron' does not exist in df.0")
}
if (all(is.na(df.0$padron)) || nrow(df.0) == 0) {
  stop("Column 'padron' in df.0 is empty or all values are NA")
}


# Ensure 'id' and 'year' columns have the same data type in both DataFrames
df.0 <- df.0 %>%
  mutate(id = as.numeric(padron),
         year = as.numeric(paste0('19',as.character(year))),
         ciiu = as.character(ciiu))

df <- df %>%
  mutate(id = as.numeric(id),
         year = as.numeric(year))

# Select relevant columns from df.0
df.0 <- df.0 %>%
  select(id, year, ciiu)

# Perform the merge
df <- merge(df, df.0, by = c("id", "year"))


ind.code <- c(311,381,321)
ind.names <- c()
for (each.code in ind.code){
  ind.each <- subset(df,ciiu_3d==each.code)
  ind.name <- ind.each$ciiu3d_descr[1]
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

for (each.code in ind.code){
  t <- Sys.time()
  ind.each <- subset(df,ciiu_3d==each.code)
  ind.name <- ind.each$ciiu3d_descr[1]
  ind.each$y <- log(ind.each$GO)
  ind.each$lnm <- log(ind.each$WI)
  ind.each$lnl <- log(ind.each$L)
  ind.each$lnk <- log(ind.each$K)
  ind.each$export <- ind.each$export # dummy for export 
  
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
  ind.each.t <- ind.each.t[complete.cases(ind.each.t),]
  ind.each.y <- cast(ind.each.t[,c("id","year","si")],id ~ year,value="si")
  id.list    <- ind.each.y[complete.cases(ind.each.y),"id"]
  #Remove the incomplete data, need balanced panel
  ind.each.t <- ind.each.t[ind.each.t$id %in% id.list,]
  ind.each.t <- ind.each.t[order(ind.each.t$id,ind.each.t$year),]
  #Reshape the Y 
  ind.each.y <- cast(ind.each.t[,c("id","year","si")],id ~ year,value="si")
  ind.each.y <- ind.each.y[,colnames(ind.each.y)!="id"]
  
  ind.each.y <- (ind.each.y - mean(ind.each.t$si))/(sd(ind.each.t$si))
  ind.each.x <- (ind.each.t$lnk - mean(ind.each.t$lnk))/(sd(ind.each.t$lnk))
  
  ind.each.ciiu_dummy <- model.matrix(~ ciiu - 1, data = ind.each.t)[,-1]
  ind.each.export <- ind.each.t$export
  ind.each.import <- ind.each.t$import
  
  data <- list(Y = t(ind.each.y), X = cbind(ind.each.x)  ,  Z = ind.each.ciiu_dummy)
  
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
}


df.3 <- data.frame(matrix('-',nrow=3*length(ind.names),ncol=10))
df.3[ 3* 1:count -2,] <- estimate.LR.df.3
df.3[ 3* 1:count -1,] <- AIC.df.3
df.3[ 3* 1:count ,] <- BIC.df.3
rownames(df.3)[ 3* 1:count -2] <- rownames(estimate.LR.df.3)
colnames(df.3) <- colnames(estimate.LR.df.3)



write.csv(df.3,file="results/Empirical/Chile_regressior_lnk_ciiu.csv")











