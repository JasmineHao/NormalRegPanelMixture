# Get the directory of the current file
# Install this.path if not already installed
if (!requireNamespace("this.path", quietly = TRUE)) {
  install.packages("this.path")
}

# Get the current file path
current_file_path <- this.path::this.path()

current_file_dir <- dirname(current_file_path)
source(file.path(current_file_dir, "functions.R"))



# Function to find the model where AIC stops decreasing
find_model_stop <- function(aic_values) {
  # Calculate differences between consecutive AIC values
  differences <- diff(aic_values)
  
  # Find the first instance where AIC starts increasing
  model_stop <- which(differences > 0)[1]
  
  # Handle cases where AIC never increases
  if (is.na(model_stop)) {
    return(length(aic_values))
  }
  
  # Return the selected model and its AIC value
  return(model_stop)
}


cl <- makeCluster(15)

getResult <- function(Data,nrep,an,cl, M, parlist){
  registerDoParallel(cl)


  aic_table <- matrix(0,nrow=nrep,ncol=M_max)
  bic_table <- matrix(0,nrow=nrep,ncol=M_max)
  mem_seq_1 <-matrix(0,nrow=nrep,ncol=1)
  mem_seq_5 <-matrix(0,nrow=nrep,ncol=1)
  mem_seq_10 <-matrix(0,nrow=nrep,ncol=1)
  lr.estim_table <- matrix(0,nrow=nrep,ncol=M_max)
  rk.mean.estim_table <- matrix(0,nrow=nrep,ncol=M_max)
  rk.max.estim_table <- matrix(0,nrow=nrep,ncol=M_max)
  
  mem_seq_rk_mean <- matrix(0,nrow=nrep,ncol=1)
  mem_seq_rk_max <-matrix(0,nrow=nrep,ncol=1)
  
  T <- nrow(Data[,1]$Y)
  T.pair.list <- pairwise_combinations(1:T)
  for (ii in 1:nrep){
    # library(NormalRegPanelMixture)
    t <- Sys.time()
    data <- Data[,ii]
    aic <- rep(0,M_max)
    bic <- rep(0,M_max)
    lr.estim <- rep(0,M_max)
    rk.mean <- rep(0, M_max)
    rk.max <- rep(0, M_max)
    
    crit.1 <- rep(0,M_max)
    crit.5 <- rep(0,M_max)
    crit.10 <- rep(0,M_max)
    crit.rk.mean <- rep(0, M_max)
    crit.rk.max <- rep(0, M_max)
    

    test <- 1
    test.rk <- 1
    for(m in 1:M_max){
      out.h0 <- normalpanelmixPMLE(y = data$Y, x = data$X, z = data$Z, m = m, vcov.method = "none")
      an <- anFormula(phi, M, N, T)      
      aic[m] <- out.h0$aic
      bic[m] <- out.h0$bic
      out.h1 <- normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an=(an),parallel = FALSE)
      lr.estim[m] <- 2 * max(out.h1$penloglik - out.h0$loglik)
      if (test) {
        crit <- try(NormalRegPanelMixture::regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, parallel = TRUE, cl=cl, nrep=1000)$crit)
        if (class(crit) == "try-error"){
          crit <- NormalRegPanelMixture::regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, an = an, z = data$Z, parallel = TRUE, cl = cl)$crit
        }
        
      } else {
        crit <- c(Inf,Inf,Inf) 
      }
      
      result_rk <- compute_rk_statistics_pairwise_T(data, T.pair.list, N, m, n.grid = (m + 1))
      
      rk.mean[m] <- mean(result_rk$rk)
      rk.max[m] <- max(result_rk$rk)
      if (test.rk) {
        stats_KP_boot <- construct_stat_KP_P_bootstrap(result_rk$P_k_list, result_rk$Sigma_P_list, T.pair.list, N, BB, result_rk$lambda_c, r.test = m, n.grid = (m + 1), transform = "P")
        crit.rk.mean[m] <- quantile(rowMeans(stats_KP_boot$rk_b), 0.95)
        crit.rk.max[m]  <- quantile(apply(stats_KP_boot$rk_b, 1, max), 0.95)
      } else {
        crit.rk.mean[m] <- Inf
        crit.rk.max[m] <- Inf
      }
      
      
      if (lr.estim[m] < crit[1]) {
        test <- 0
      }
      if (rk.mean[m] < crit.rk.mean[m]){
        test.rk <- 0 
      }
      
      crit.1[m] <- crit[3]
      crit.5[m] <- crit[2]
      crit.10[m] <- crit[1]
      
    }
    
    lr.estim_table[ii, ] <- lr.estim
    rk.mean.estim_table[ii, ] <- rk.mean
    rk.max.estim_table[ii, ] <- rk.max

    aic_table[ii,] <- aic
    bic_table[ii,] <- bic
    mem_seq_1[ii,] <- min(which(lr.estim < crit.1))
    mem_seq_5[ii,] <- min(which(lr.estim < crit.5))
    mem_seq_10[ii, ] <- min(which(lr.estim < crit.10))
    mem_seq_rk_mean[ii, ] <- min(which(rk.mean < crit.rk.mean))
    mem_seq_rk_max[ii, ] <- min(which(rk.max < crit.rk.max))

      
    print(paste("simulation", ii, "out of", nrep))
    print(Sys.time() - t) 
 }

  # Find the model where AIC stops decreasing (local minimum)
  differences <- diff(aic_values)  # Calculate differences between consecutive AIC values
  model_stop <- which(differences > 0)[1]  # Find the first instance where AIC starts increasing

  
  aic_freq <- apply(aic_table,1,find_model_stop)
  bic_freq <- apply(bic_table,1,find_model_stop)
  
  return(list(aic=aic_freq,bic=bic_freq,mem_seq_1=mem_seq_1,mem_seq_5=mem_seq_5,mem_seq_10=mem_seq_10, rk.mean.estim_table=rk.mean.estim_table,rk.max.estim_table=rk.max.estim_table, mem_seq_rk_mean=mem_seq_rk_mean, mem_seq_rk_max = mem_seq_rk_max))
}

count_freq <- function(stats){
  stats_freq_tab <- rep(0,7)
  for(m in 1:7){
    stats_freq_tab[m] = sum(stats == m)
  }
  stats_freq_tab <- stats_freq_tab/length(stats)
  return(stats_freq_tab)
}

# Setting the parameters

# Parset = list(list(c(0.1547859, 0.4820520, 0.3631621),c(-1.0984991, -0.2158495, 0.7547124),c(1.5474731, 0.4887465, 0.5186135)))

Parset = list(list(c(0.1172174, 0.4086027, 0.4741799),c(-1.6401041, -0.3275637, 0.6876973),c(1.0428541, 0.5219205, 0.6094685)))

N <- 196
T <- 3
M <- 3 #Number of Type
p <- 0 #Number of Z
q <- 0 #Number of X
M_max <- 7

count <- 1
nrep <- 100
BB <- 199

alpha = Parset[[count]][[1]]
mu = Parset[[count]][[2]]
sigma = Parset[[count]][[3]]

t <- Sys.time()
phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = NULL,
            beta = NULL, N = N, T = T, M = M, p = p, q = q, X=NULL)

phi.data.pair <- GenerateSample(phi,nrep)
Data = phi.data.pair$Data
phi = phi.data.pair$phi

result <- getResult(Data,nrep,an,cl,M, parlist)


result_freq_table <- rbind(count_freq(result$aic),
count_freq(result$bic),
count_freq(result$mem_seq_1),
count_freq(result$mem_seq_5),
count_freq(result$mem_seq_10),
count_freq(result$mem_seq_rk_mean),
count_freq(result$mem_seq_rk_max))

rownames(result_freq_table) <- c("aic", "bic", "lr 1%", "lr 5%", "lr 10%", "rk mean 5%", "rk max 5%")
colnames(result_freq_table) <- c("M=1","M=2","M=3","M=4","M=5","M=6","M=7")

write.csv(100 * result_freq_table, "empirical_test_M3.csv", row.names = TRUE)

