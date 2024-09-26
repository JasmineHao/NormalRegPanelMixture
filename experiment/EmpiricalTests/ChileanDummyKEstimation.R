library(stargazer)
library(reshape)
# library(normalregMix)
library(foreign)
library(NormalRegPanelMixture)
library(haven)
library(dplyr)
library(splines)


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


bs_degree <- 2

count <- 1

# model <- "k_dummy_imex_ciiu"
# model_M_list <- c(7,6,6)

model <- "k_dummy_ciiu"
model_M_list <- c(7,7,9)

# model <- "k_dummy_imex"
# model_M_list <- c(7,5,5)

# model <- "k_dummy"
# model_M_list <- c(9,8,9)

# candidate models
# c("k_dummy", "k_dummy_imex", "k_dummy_ciiu", ”k_dummy_imex_ciiu“)

# for (each.code in ind.code)
for (each.code in ind.code)  {
  t <- Sys.time()
  print(each.code)
  ind.each <- subset(df,ciiu_3d==each.code)
  ind.name <- ind.each$ciiu3d_descr[1]
  ind.each$y <- log(ind.each$GO)
  ind.each$lnm <- log(ind.each$WI)
  ind.each$lnl <- log(ind.each$L)
  ind.each$lnk <- log(ind.each$K)
  ind.each$export <- ind.each$export # dummy for export 
  ind.each$import <- ind.each$import # dummy for import
  
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
  ind.each.k_hl   <- ind.each.x > median(ind.each.x) # indicator for high capital v.s. low capital
  
  
  
  # data <- list(Y = t(ind.each.y), X = cbind(ind.each.x)  ,  Z = ind.each.ciiu_dummy)
  # Find all columns in ind.each that contain '_bs_'
  bs_columns <- grep("_bs_", colnames(ind.each.t), value = TRUE)
  
  # c("k_dummy", "k_dummy_imex", "k_dummy_ciiu", ”k_dummy_imex_ciiu“)
  if (model == "k_dummy"){
    data <- list(Y = t(ind.each.y), X = cbind(ind.each.k_hl*1) ,  Z = NULL)
    
  }else if (model == "k_dummy_imex"){
    data <- list(Y = t(ind.each.y), X = cbind(ind.each.k_hl*1, ind.each.t[, c("import","export")])  ,  Z = NULL)
    
  }else if (model == "k_dummy_ciiu"){
    data <- list(Y = t(ind.each.y), X = cbind(ind.each.k_hl*1) ,  Z = ind.each.ciiu_dummy)
    
  }else if (model == "k_dummy_imex_ciiu"){
    data <- list(Y = t(ind.each.y), X = cbind(ind.each.k_hl*1, ind.each.t[, c("import","export")])  ,  Z = ind.each.ciiu_dummy)
  }
  
  if (model == "k_dummy" || model == "k_dummy_ciiu"){
    x.type <- expand.grid(c(0,1))
    
    data_X <- data.frame(data$X)
    x.type <- as.matrix(x.type)
    
    colnames(data_X) <- c("A")
    # Use dplyr to assign group numbers by matching the rows
    df_grouped <-  data_X %>%
      left_join(data.frame(x.type) %>% mutate(group_num = row_number()), 
                by = c("A" = "Var1"))
    
  } else if (model == "k_dummy_imex" || model == "k_dummy_imex_ciiu"){
    x.type <- expand.grid(c(0,1),c(0,1),c(0,1))
    
    data_X <- data.frame(data$X)
    x.type <- as.matrix(x.type)
    
    colnames(data_X) <- c("A", "B", "C")
    # Use dplyr to assign group numbers by matching the rows
    df_grouped <-  data_X %>%
      left_join(data.frame(x.type) %>% mutate(group_num = row_number()), 
                by = c("A" = "Var1", "B" = "Var2", "C" = "Var3"))
    
  }
  
  N <- dim(ind.each.y)[1]
  
  
  h1.coefficient = NULL
  
  estimate.crit <- 1
  
  # Estimate the null model
  M <- model_M_list[count]
  out.h0 <- regpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "Hessian",in.coefficient=NULL, ninits=5)
  
  q       <- ncol(data$X)
  alpha   <- out.h0$parlist$alpha  # m-vector
  mubeta  <- out.h0$parlist$mubeta  # q+1 by m
  sigma   <- out.h0$parlist$sigma  # m-vector
  gam     <- out.h0$parlist$gam
  
  
  estim_df <- array(0, dim=c(M, nrow(x.type)))
  var_df   <- array(0, dim=c(M, nrow(x.type)))
  
  # store the post probability dim(out.h0$postprobs)
  
  # Convert out.h0$postprobs to a dataframe
  postprobs_df <- as.data.frame(out.h0$postprobs)
  postprobs_df <- postprobs_df[rep(1:nrow(postprobs_df), each = 3), ]
  
  # Bind df_grouped$group_num as a new column to the dataframe
  postprobs_df$group_num <- df_grouped$group_num
  
  # Group by group_num and calculate the mean for each group
  mean_by_type <- postprobs_df %>%
    group_by(group_num) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
  
  # Print the result
  print(mean_by_type)
  
  for (m in 1:M){
    for (t in 1:nrow(x.type)){
      mubeta_ind <- array(0, dim = dim(mubeta))
      mubeta_ind[1,m] <- 1
      mubeta_ind[2:(q+1),m] <- as.vector(x.type[t,])
      coef_ind = array(0, dim=length(out.h0$coefficients))
      coef_ind[(M+1):((2+q)*M)] <- as.vector(mubeta_ind )
      if (length(gam) > 0){
        coef_ind[((3+q)*M + 1):length(coef_ind)] <- 1 / length(gam)
      }
      estim_df[m,t] <- coef_ind %*% out.h0$coefficients
      var_df[m,t] <- coef_ind %*% out.h0$vcov %*% coef_ind
    }
  }
  
  # mubeta  <- matrix(coefficients[(m+1):((2+q)*m)], nrow=q+1, ncol=m)  # q+1 by m
  # gam <- [((2+q)*m+1):((3+q)*m)]
  # out.h0$coefficients[(M+1):((2+q)*M)]
  
  count = count + 1
  write.csv(mean_by_type, paste("results/Empirical/Additional/",model, "_", ind.name ,"_type_prob.csv",sep=""))
  
  write.csv(cbind(estim_df,var_df), 
            paste("results/Empirical/Additional/",model, "_", ind.name ,"_statistics.csv",sep=""))

  stargazer(cbind(estim_df,var_df), type = "text")
  
}



