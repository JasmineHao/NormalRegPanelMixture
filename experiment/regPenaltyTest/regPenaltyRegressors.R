library(stargazer)
library(NormalRegPanelMixture)
#Get the misclassfication
GetMisclTerm <- function(phi) {
  m <- phi$M
  
  if (m == 2)
  {
    omega.12  <- omega.12(phi)
    return (log(omega.12 /(1-omega.12)))
  }
  
  if (m == 3) # ln(omega_12 omega_23 / (0.5-omega_12)(0.5-omega_23))
  {
    omega.123 <- omega.123(phi)
    omega.12 <- omega.123[1]
    omega.23 <- omega.123[2]
    return (log(omega.12 * omega.23 / ((1-omega.12)*(1-omega.23))))
  }
  omega.1234 <- omega.1234(phi)
  omega.12 <- omega.1234[1]
  omega.23 <- omega.1234[2]
  omega.34 <- omega.1234[3]
  # (m == 4) # ln(omega_12 omega_23 omega_34 / (0.5-omega_12)(0.5-omega_23)(0.5-omega_34))
  return (log(omega.12 * omega.23 * omega.34 /
                ((1-omega.12)*(1-omega.23)*(1-omega.34))))
  
}

###########################################################
df.M1 <- read.csv(file="results/regPenaltyTest/penaltyTestM1.csv")
df.M2 <- read.csv(file="results/regPenaltyTest/penaltyTestM2.csv")
df.M3 <- read.csv(file="results/regPenaltyTest/penaltyTestM3.csv")
df.M4 <- read.csv(file="results/regPenaltyTest/penaltyTestM4.csv")

###########################################################
# Specification 1: 1 /T and 1 / N
###########################################################
# df.M2 <- df.M2[df.M2$an < 0.4,]

fit.data <- list(y = log(
  df.M1$nom.size/ (1 - df.M1$nom.size)) - log(5/95)  ,
  t = 1 /  df.M1$T,
  n = 1 /  df.M1$N,
  an = log( df.M1$an/ (1 - df.M1$an )))

an.model.M1 <- lm(y ~ t + n + an  , data=fit.data,na.action = na.omit)
coef.M1 <- an.model.M1$coefficients
print(coef.M1)
an.formula.M1 <- (coef.M1[1] +coef.M1[2]*fit.data$t + coef.M1[3]*fit.data$n )/coef.M1[4]
an.formula.M1 <- 1 /(1 + exp(an.formula.M1))   
m_1.an <- 1 * mean(an.formula.M1)


fit.data <- list(y = log(
  df.M2$nom.size/ (1 - df.M2$nom.size)) - log(5/95)  ,
  t = 1 /  df.M2$T,
  n = 1 /  df.M2$N,
  an = log( df.M2$an/ (1 - df.M2$an )),
  omega = df.M2$omega)

an.model.M2 <- lm(y ~ t + n + an + omega , data=fit.data,na.action = na.omit)
coef.M2 <- an.model.M2$coefficients
print(coef.M2)
an.formula.M2 <- (coef.M2[1] +coef.M2[2]*fit.data$t + coef.M2[3]*fit.data$n + coef.M2[5]*fit.data$omega)/coef.M2[4]
an.formula.M2 <- 1 /(1 + exp(an.formula.M2))   
m_2.an <- 1 * mean(an.formula.M2)




fit.data <- list(y = log(
  df.M3$nom.size/ (1 - df.M3$nom.size)) - log(5/95)  ,
  t = 1 /  df.M3$T,
  n = 1 /  df.M3$N,
  an = log( df.M3$an/ (1 - df.M3$an )),
  omega = df.M3$omega)
an.model.M3 <- lm(y ~ t + n + an + omega , data=fit.data,na.action = na.omit)
coef.M3 <- an.model.M3$coefficients
print(coef.M3)
an.formula.M3 <- (coef.M3[1] +coef.M3[2]*fit.data$t + coef.M3[3]*fit.data$n + coef.M3[5]*fit.data$omega)/coef.M3[4]
an.formula.M3 <- 1 /(1 + exp(an.formula.M3))   
m_3.an <- 1 * mean(an.formula.M3)


fit.data <- list(y = log(
  df.M4$nom.size/ (1 - df.M4$nom.size)) - log(5/95)  ,
  t = 1 /  df.M4$T,
  n = 1 /  df.M4$N,
  an = log( df.M4$an/ (1 - df.M4$an )),
  omega = df.M4$omega)
an.model.M4 <- lm(y ~ t + n + an + omega , data=fit.data,na.action = na.omit)
coef.M4 <- an.model.M4$coefficients
print(coef.M4)
an.formula.M4 <- (coef.M4[1] +coef.M4[2]*fit.data$t + coef.M4[3]*fit.data$n + coef.M4[5]*fit.data$omega)/coef.M4[4]
an.formula.M4 <-  1 /(1 + exp(an.formula.M3))    
m_4.an <- mean(an.formula.M4)


stargazer(an.model.M1, an.model.M2,an.model.M3,an.model.M4)



###########################################################
# Specification 2: log(T) and log(N)
###########################################################

fit.data <- list(y = log(
  df.M2$nom.size/ (1 - df.M2$nom.size) ) ,
  t = log(df.M2$T),
  n = log(df.M2$N),
  an = log( df.M2$an/ (1 - df.M2$an )),
  omega = df.M2$omega)
an.model.M2.2 <- lm(y ~ t + n + an + omega , data=fit.data,na.action = na.omit)
print(an.model.M2.2$coefficients)

c(-3.47934199,  0.14566884, -0.09219391, -0.13457721, -0.05850536)

###########################################################
# Specification 3: sqrt(T) and sqrt(N)
###########################################################

fit.data <- list(y = log(
  df.M2$nom.size/ (1 - df.M2$nom.size) ) ,
  t = 1 / sqrt(df.M2$T),
  n = 1/ sqrt(df.M2$N),
  an = log( df.M2$an/ (1 - df.M2$an )),
  omega = df.M2$omega)
an.model.M2.3 <- lm(y ~ t + n + an + omega , data=fit.data,na.action = na.omit)
print(an.model.M2.3$coefficients)

c(-3.69387007, -0.51991487,  2.68422631, -0.13457721, -0.05850536 )

stargazer(an.model.M2,an.model.M2.2,an.model.M2.3)

# coef.M2 <- an.model.M2$coefficients
# an.formula.M2 <- (log(5/95) - coef.M2[1] - coef.M2[2]*fit.data$t - coef.M2[3]*fit.data$n - coef.M2[4]*fit.data$an)/coef.M2[5]
# an.formula.M2 <- exp(an.formula.M2)/(1 + exp(an.formula.M2))   
# m_2.an <- 0.5 * mean(an.formula.M2)

# fit.data$y <- replace(fit.data$y,fit.data$y == Inf,5)
# fit.data$y <- replace(fit.data$y,is.na(fit.data$y),5)


# M <- 2
# p <- 0
# q <- 0
# Nset <- c(100,500)
# Tset <- c(2,5,10)
# alphaset <- list(c(0.5,0.5),c(0.2,0.8))
# muset <- list(c(-1,1),c(-0.5,0.5))
# sigmaset <- list(c(1, 1), c(1.5, 0.75))
# nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset) * length(sigmaset)
# an.data <- matrix(0,nr=nset,nc=3)

# count <- 0
# for (N in Nset){
#   for (T in Tset){
#     for (mu in muset){
#       for (alpha in alphaset){
#         for (sigma in sigmaset){
#           count <- count + 1
#           phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
#                      beta = beta, N = N, T = T, M = M, p = p, q = q)
#           # misclterm  <- GetMisclTerm(phi)
#           # print(misclterm)
#           # an <- anFormula(phi,M,N,T) * 0.5 
#           # an.data[count, ] <-
#           # fit.data$omega[count] <- misclterm
#           
#         }
#       }
#     }
#   }
# }




###########################################################

