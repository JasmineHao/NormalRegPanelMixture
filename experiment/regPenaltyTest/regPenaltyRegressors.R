library(stargazer)
#Get the misclassfication
GetMisclTerm <- function(phi) {
  m <- phi$M
  
  if (m == 2)
  {
    omega.12  <- omega.12(phi)
    return (log(omega.12 /(0.5-omega.12)))
  }
  
  if (m == 3) # ln(omega_12 omega_23 / (0.5-omega_12)(0.5-omega_23))
  {
    omega.123 <- omega.123(phi)
    omega.12 <- omega.123[1]
    omega.23 <- omega.123[2]
    return (log(omega.12 * omega.23 / ((0.5-omega.12)*(0.5-omega.23))))
  }
  omega.1234 <- omega.1234(phi)
  omega.12 <- omega.1234[1]
  omega.23 <- omega.1234[2]
  omega.34 <- omega.1234[3]
  # (m == 4) # ln(omega_12 omega_23 omega_34 / (0.5-omega_12)(0.5-omega_23)(0.5-omega_34))
  return (log(omega.12 * omega.23 * omega.34 /
                ((0.5-omega.12)*(0.5-omega.23)*(0.5-omega.34))))
  
}

###########################################################

df.M2 <- read.csv(file="results/regPenaltyTest/penaltyTestM2.csv")
fit.data <- list(y = log(
  df.M2$nom.size/ (0.25 - df.M2$nom.size) ) ,
  t = 1 /  df.M2$T,
  n = 1 /  df.M2$N,
  an = log( df.M2$an/ (1 - df.M2$an )),
  omega = log(df.M2$omega / (1 - df.M2$omega)))
# fit.data$y <- replace(fit.data$y,fit.data$y == Inf,5)
# fit.data$y <- replace(fit.data$y,is.na(fit.data$y),5)


M <- 2
p <- 0
q <- 0
Nset <- c(100,500)
Tset <- c(2,5,10)
alphaset <- list(c(0.5,0.5),c(0.2,0.8))
muset <- list(c(-1,1),c(-0.5,0.5))
sigmaset <- list(c(1, 1), c(1.5, 0.75))
nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset) * length(sigmaset)
# an.data <- matrix(0,nr=nset,nc=3)

count <- 0
for (N in Nset){
  for (T in Tset){
    for (mu in muset){
      for (alpha in alphaset){
        for (sigma in sigmaset){
          count <- count + 1
          phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
                     beta = beta, N = N, T = T, M = M, p = p, q = q)
          misclterm  <- GetMisclTerm(phi)          
          # an <- anFormula(phi,M,N,T) * 0.5 
          # an.data[count, ] <-
          fit.data$omega[count] <- log(misclterm / (1 - misclterm))
          
        }
      }
    }
  }
}

an.model.M2 <- lm(y ~ t + n + an + omega , data=fit.data,na.action = na.omit)
print(summary(an.model.M2))
coef.M2 <- an.model.M2$coefficients
an.formula.M2 <- (-1.386294 - coef.M2[1] - coef.M2[2]*fit.data$t - coef.M2[3]*fit.data$n - coef.M2[4]*fit.data$an)/coef.M2[5]
an.formula.M2 <- exp(an.formula.M2)/(1 + exp(an.formula.M2))   
m_2.an <- 0.25 * mean(an.formula.M2)



###########################################################



df.M4 <- read.csv(file="results/regPenaltyTest/penaltyTestM4.csv")
fit.data <- list(y = log(
  df.M4$nom.size/ (0.25 - df.M4$nom.size) ) ,
  t = 1 /  df.M4$T,
  n = 1 /  df.M4$N,
  an = log( df.M4$an/ (1 - df.M4$an )),
  omega = log(df.M4$omega / (1 - df.M4$omega)))
# fit.data$y <- replace(fit.data$y,fit.data$y == Inf,5)
# fit.data$y <- replace(fit.data$y,is.na(fit.data$y),5)
an.model.M4 <- lm(y ~ t + n + an + omega , data=fit.data,na.action = na.omit)
print(summary(an.model.M4))


###########################################################
M <- 3
alphaset <- list(c(1/3,1/3,1/3),c(0.25,0.5,0.25))
muset 		<- list(c(-4, 0, 4), c(-4, 0, 5), c(-5, 0, 5), c(-4, 0, 6), c(-5, 0, 6), c(-6, 0, 6))
sigmaset <- list(c(1, 1, 1), c(0.75, 1.5, 0.75))

nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset) * length(sigmaset)
an.data <- matrix(0,nr=nset,nc=3)

count <- 0
for (N in Nset){
  for (T in Tset){
    for (mu in muset){
      for (alpha in alphaset){
        for (sigma in sigmaset){
          count <- count + 1
          phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
                     beta = beta, N = N, T = T, M = M, p = p, q = q)
          
          an <- anFormula(phi,M,N,T) * 0.5 
          an.data[count, ] <-
            cbind(an, phi$N,phi$T)
          
        }
      }
    }
  }
}

m_3.an <- mean(an.data[,1])


###########################################################
df.M4 <- read.csv(file="results/regPenaltyTest/penaltyTestM4.csv")
fit.data <- list(y = log(
  df.M4$nom.size/ (1 - df.M4$nom.size) ) ,
  t = 1 /  df.M4$T,
  n = 1 /  df.M4$N,
  an = log( df.M4$an/ (1 - df.M4$an )),
  omega = log(df.M4$omega / (1 - df.M4$omega)))
fit.data$y <- replace(fit.data$y,fit.data$y == Inf,5)
fit.data$y <- replace(fit.data$y,is.na(fit.data$y),5)
# an.model.M4 <- lm(y ~ t + n + an + omega , data=fit.data,na.action = na.omit)
# print(summary(an.model.M4))

M <- 4
alphaset <- list(c(0.25, 0.25, 0.25, 0.25))
muset 		<- list(c(-4,-1,1,4), c(-5,-1,1,5), c(-6,-2,2,6), c(-6,-1,2,5), c(-5,0,2,4), c(-6,0,2,4))
sigmaset <- list(c(1, 1, 1, 1), c(1, 0.75, 0.5, 0.25))

nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset) * length(sigmaset)
an.data <- matrix(0,nr=nset,nc=3)

count <- 0
for (N in Nset){
  for (T in Tset){
    for (mu in muset){
      for (alpha in alphaset){
        for (sigma in sigmaset){
          count <- count + 1
          phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
                     beta = beta, N = N, T = T, M = M, p = p, q = q)
          
          misclterm  <- GetMisclTerm(phi)          
          # an <- anFormula(phi,M,N,T) * 0.5 
          # an.data[count, ] <-
          fit.data$omega[count] <- misclterm
          
        }
      }
    }
  }
}

an.model.M4 <- lm(y ~ t + n + an + omega , data=fit.data,na.action = na.omit)
print(summary(an.model.M4))
coef.M4 <- an.model.M4$coefficients
an.formula.M4 <- (-2.9 - coef.M4[1] - coef.M4[2]*fit.data$t - coef.M4[3]*fit.data$n - coef.M4[4]*fit.data$an)/coef.M4[5]
an.formula.M4 <- exp(an.formula.M4)/(1 + exp(an.formula.M4))   
m_4.an <- 0.25 * mean(an.formula.M4)
print(m_4.an)
