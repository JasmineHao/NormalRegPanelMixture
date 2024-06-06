library(stargazer)
library(ggplot2)
library(reshape)
# library(normalregMix)
library(foreign)
library(NormalRegPanelMixture)
options('nloptr.show.inequality.warning'=FALSE)
options(warn = -1)

set.seed(123)

#df <- readRDS("/home/haoyu/NormalRegPanelMixture/data/ChileanClean.rds")

df <- readRDS("data/ChileanClean.rds")
cl <- makeCluster(8)

ind.code <- c(311,381,321,322,331,356,342,382,352,369,324)
ind.code <- c(311,381,321)

ind.names <- c()
for (each.code in ind.code){
  ind.each <- subset(df,ciiu_3d==each.code)
  ind.name <- ind.each$ciiu3d_descr[1]
  ind.names <- append(ind.names,ind.name)
}

ind.count <- length(ind.code)

estimate.LR.df <- matrix(0,nr=length(ind.code),nc=10)
rownames(estimate.LR.df) <- ind.names
colnames(estimate.LR.df) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")

AIC.df <- matrix(0,nr=length(ind.code),nc=10)
rownames(AIC.df) <- ind.names
colnames(AIC.df) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")

BIC.df <- matrix(0,nr=length(ind.code),nc=10)
rownames(BIC.df) <- ind.names
colnames(BIC.df) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")

Nonpar.df <- matrix(0, nr = length(ind.code), nc = 10)
rownames(Nonpar.df) <- ind.names
colnames(Nonpar.df) <- c("M=1", "M=2", "M=3", "M=4", "M=5", "M=6", "M=7", "M=8", "M=9", "M=10")

Nonpar.pval.df <- matrix(0, nr = length(ind.code), nc = 10)
rownames(Nonpar.pval.df) <- ind.names
colnames(Nonpar.pval.df) <- c("M=1", "M=2", "M=3", "M=4", "M=5", "M=6", "M=7", "M=8", "M=9", "M=10")

crit.LR.df <- matrix(0,nr=length(ind.code),nc=10)
rownames(crit.LR.df) <- ind.names
colnames(crit.LR.df) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")

count = 0
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
  # desc.each <- desc.each[complete.cases(desc.each),]
  year.list <- sort(unique(ind.each$year))
  T.cap <- max(year.list)
  
  coef.df <- matrix(0,nr=5,nc=10)
  estimate.df <- matrix(0,nr=5,nc=10)
  AIC.df.each <- matrix(0,nr=5,nc=10)
  BIC.df.each <- matrix(0, nr = 5, nc = 10)
  Nonpar.df.each <- matrix(0, nr = 5, nc = 10)
  Nonpar.pval.df.each <- matrix(0, nr = 5, nc = 10)
  crit.df <- matrix(0, nr = 5, nc = 10)
  result.df <- matrix(0,nr=5,nc=10)
  ######################################################
  #For panel data
  ######################################################
  
  
  T <- 4  
  
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
  
  data <- list(Y = t(ind.each.y), X = NULL,  Z = NULL)
  
  # Create a sequence from 1 to T
  # --------------------------------
  T.sequence <- 1:T
  # Partition the sequence into even and odd numbers
  T.even <- T.sequence[T.sequence %% 2 == 0]
  T.odd <- T.sequence[T.sequence %% 2 == 1]
  
  
  # --------------------------------
  
  N <- dim(ind.each.y)[1]
  
  h1.coefficient = NULL
  estimate.crit <- 1
  estimate.crit.nonpar <- 1
  for (M in 1:10){
    # Estimate the null model
    out.h0 <- normalpanelmixPMLE(y=data$Y, z = data$Z,m=M,vcov.method = "none",in.coefficient=h1.coefficient)
    an <- anFormula(out.h0$parlist,M,N,T,q=1)
    print("-----------------------------------------")
    print(paste("T=",T,"M = ",M,"an=",an))
    if (is.na(an)){
      an <- 1.0
    }
    # Estimate the alternative model
    out.h1 <- normalpanelmixMaxPhi(y=data$Y, z = data$Z,parlist=out.h0$parlist,an=an)
    h1.parlist = out.h1$parlist
    
    lr.estimate <- 2 * max(out.h1$penloglik - out.h0$loglik)
    
    # Simulate the asymptotic distribution
    if (estimate.crit == 1){
      lr.crit <- try(regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, cl=cl,parallel = TRUE)$crit)
      if (class(lr.crit) == "try-error"){
        print("BOOTSTRAP")
        lr.crit <- regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, cl=cl,parallel = TRUE)$crit
      }
    }
    # Store the estimation results
    coef.df[T,M] <- paste(paste(names(out.h0$coefficients), collapse = ","), paste(out.h0$coefficients, collapse = ","), collapse = ",")
    estimate.df[T,M] <- paste('$',round(lr.estimate,2),'^{',paste(rep('*',sum(lr.estimate > lr.crit)),  collapse = ""),'}','$', sep = "")
    AIC.df.each[T,M] <- round(out.h0$aic,2)
    BIC.df.each[T,M] <- round(out.h0$bic,2)
    
    # --------------------------------
    # Nonpar
    if (estimate.crit.nonpar){
      data_P_W <- calculate_W_P(data,T.even, T.odd, n.grid=(ceiling(sqrt(M)) + 1), BB=199, type="indicator")
      Nonpar.df.each[T,M] <- construct_stat_KP(data_P_W$P_c,data_P_W$W_c, M, N)
      Nonpar.pval.df.each[T,M] <- return_p_val(data_P_W, M, N)
      if (Nonpar.pval.df.each[T,M] > 0.05){
        estimate.crit.nonpar <- 0
      }
    }
    
    # --------------------------------
    
    crit.df[T,M] <- paste(round(lr.crit,2),collapse = ",")
    # If fail to reject the test, break the loop
    print(lr.estimate)
    print(lr.crit)
    
    if (sum(lr.estimate > lr.crit) < 2){
      estimate.crit <- Inf
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
  estimate.LR.df[count,] <- estimate.df[4,]
  AIC.df[count,] <- AIC.df.each[4,]
  BIC.df[count, ] <- BIC.df.each[4, ]
  Nonpar.df[count, ] <- Nonpar.df.each[4, ]
  Nonpar.pval.df[count, ] <- Nonpar.pval.df.each[4, ]
  
  crit.LR.df[count,] <- crit.df[4,]
  colnames(estimate.df) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
  rownames(estimate.df) <- c("T=1","T=2","T=3","T=4","T=5")
  
  colnames(crit.df) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
  rownames(crit.df) <- c("T=1","T=2","T=3","T=4","T=5")
  
}



# colnames(crit.df.boot) <-  c("M=1","M=2","M=3","M=4","M=5", "M=6","M=7","M=8","M=9","M=10")
# rownames(crit.df.boot) <- c("T=1","T=2","T=3","T=4","T=5")

count <- length(ind.names)


# Combine the data frames
combined_df <- rbind(
  cbind(estimate.LR.df, original_df = "LR"),
  cbind(AIC.df, original_df = "AIC"),
  cbind(BIC.df, original_df = "BIC"),
  cbind(Nonpar.df, original_df = "Nonpar"),
  cbind(Nonpar.pval.df, original_df = "Nonpar pval")
)

write.csv(combined_df,file="results/Empirical/Chile_plain_nonpar.csv")


