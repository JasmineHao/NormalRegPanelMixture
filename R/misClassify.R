#' @description Computes omega_{j|i} defined in (2.1) of Maitra and Melnykov (2010)
#' @export
#' @title omega.ji
#' @name omega.ji
#' @param phi_i 3 by 1 column consisting of alpha, mu, sigma of ith component
#' @param phi_j 3 by 1 column consisting of alpha, mu, sigma of jth component
#' @return omega_{j|i}
#' @references Maitra, R., and Melnykov, V. (2010)
#' Simulating Data to Study Performance of Finite Mixture Modeling and Model-Based Clustering Algorithms,
#' \emph{Journal of Computational and Graphical Statistica},
#' \bold{19}, 354--376.
# Returns a misclassification rate omega_ji given two components i, j,
# i.e. the probability of choosing component j where
# the true model is ith component.
omega.ji <- function(phi_i, phi_j) {
  alpha_i <- phi_i[1]
  alpha_j <- phi_j[1]
  mu_i <- phi_i[2]
  mu_j <- phi_j[2]
  sigma_i <- phi_i[3]
  sigma_j <- phi_j[3]
  
  a <- (1/sigma_j^2 - 1/sigma_i^2)
  b <- mu_i / sigma_i^2 - mu_j / sigma_j^2
  c <- mu_j^2 / sigma_j^2 - mu_i^2 / sigma_i^2
  
  if (sigma_i == sigma_j){
    if (mu_i > mu_j){
      omega_ji = pnorm((2 * log(alpha_j/alpha_i) - c)/(2*b),
                       mean = mu_i, sd = sigma_i)
    }
  else{
    omega_ji = 1 - pnorm((2 * log(alpha_j/alpha_i) - c)/(2*b),
                         mean = mu_i, sd = sigma_i)}
  }else{
    d <- 2 * log(alpha_j * sigma_i / (alpha_i * sigma_j)) - c + (b^2 / a)
    da <- max(d/a, 0)
    if (sigma_i > sigma_j){
      omega_ji = pnorm(sqrt(da)-b/a, mean = mu_i, sd = sigma_i) -
      pnorm(-sqrt(da)-b/a, mean = mu_i, sd = sigma_i)}
    else{
      omega_ji = 1 +
      pnorm(-sqrt(da)-b/a, mean = mu_i, sd = sigma_i) -
      pnorm(sqrt(da)-b/a, mean = mu_i, sd = sigma_i)}
  }
  return (omega_ji)
}

#' @description Computes omega_{12} defined in Maitra and Melnykov (2010)
#' @export
#' @title omega.12
#' @name omega.12
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gam
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gam_1, ..., gam_m))
#' @return The misclassification rate omega_ij
#' @references Maitra, R., and Melnykov, V. (2010)
#' Simulating Data to Study Performance of Finite Mixture Modeling and Model-Based Clustering Algorithms,
#' \emph{Journal of Computational and Graphical Statistica},
#' \bold{19}, 354--376.
omega.12 <- function(parlist)
  # Computes omega_{12} for testing H_0:m=2 against H_1:m=3
{
  phi1 <- c(alpha = parlist$alpha[1], mu = parlist$mu[1], sigma = parlist$sigma[1])
  phi2 <- c(alpha = parlist$alpha[2], mu = parlist$mu[2], sigma = parlist$sigma[2])
  
  part1 <- omega.ji(phi1, phi2)
  part2 <- omega.ji(phi2, phi1)
  
  return((part1 + part2) / 2)
}  # end function omega.12


#' Computes omega_{12} and omega_{23} defined in Maitra and Melnykov (2010)
#' @export
#' @title omega.123
#' @name omega.123
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gamma_1, ..., gamma_m))
#' @return A 2 by 1 vector whose first element is omega_12 and second element is omega_23
#' @references Maitra, R., and Melnykov, V. (2010)
#' Simulating Data to Study Performance of Finite Mixture Modeling and Model-Based Clustering Algorithms,
#' \emph{Journal of Computational and Graphical Statistica},
#' \bold{19}, 354--376.
omega.123 <- function(parlist)
{
  phi1 <- c(alpha = parlist$alpha[1], mu = parlist$mu[1], sigma = parlist$sigma[1])
  phi2 <- c(alpha = parlist$alpha[2], mu = parlist$mu[2], sigma = parlist$sigma[2])
  phi3 <- c(alpha = parlist$alpha[3], mu = parlist$mu[3], sigma = parlist$sigma[3])
  
  part1 <- omega.ji(phi1, phi2)
  part2 <- omega.ji(phi2, phi1)
  w12 <- (part1 + part2)/2
  
  part3 <- omega.ji(phi2, phi3)
  part4 <- omega.ji(phi3, phi2)
  w23 <- (part3 + part4)/2
  
  return(c(w12, w23))
  
}  # end function omega.123

#' @description Computes omega_{12}, omega_{23}, and omega_{34} defined in Maitra and Melnykov (2010)
#' @export
#' @title omega.1234
#' @name omega.1234
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gamma_1, ..., gamma_m))
#' @return A 3 by 1 vector consisting of omega_12, omega_23, and omega_34
#' @references Maitra, R., and Melnykov, V. (2010)
#' Simulating Data to Study Performance of Finite Mixture Modeling and Model-Based Clustering Algorithms,
#' \emph{Journal of Computational and Graphical Statistica},
#' \bold{19}, 354--376.
omega.1234 <- function(parlist)
{
  phi1 <- c(alpha = parlist$alpha[1], mu = parlist$mu[1], sigma = parlist$sigma[1])
  phi2 <- c(alpha = parlist$alpha[2], mu = parlist$mu[2], sigma = parlist$sigma[2])
  phi3 <- c(alpha = parlist$alpha[3], mu = parlist$mu[3], sigma = parlist$sigma[3])
  phi4 <- c(alpha = parlist$alpha[4], mu = parlist$mu[4], sigma = parlist$sigma[4])
  
  part1 <- omega.ji(phi1, phi2)
  part2 <- omega.ji(phi2, phi1)
  w12 <- (part1 + part2)/2
  
  part3 <- omega.ji(phi2, phi3)
  part4 <- omega.ji(phi3, phi2)
  w23 <- (part3 + part4)/2
  
  part5 <- omega.ji(phi3, phi4)
  part6 <- omega.ji(phi4, phi3)
  w34 <- (part5 + part6)/2
  
  return(c(w12, w23, w34))
  
}  # end function omega.1234

coef.to.list <- function(coefficients, z = NULL) {
  # Convert coefficients to list
  len     <- length(coefficients)
  p       <- 0
  gam   <- NULL
  
  if (!is.null(z)) {
    z <- as.matrix(z)
    p <- ncol(z)
    gam <- coefficients[(len-p+1):len]
  }
  
  m <- (len-p)/3
  if (round(m) != m) {
    stop("The dimension of the coefficients is incompatible with z. Please check the data.")
  }
  
  param   <- matrix(coefficients[1:(len-p)], nrow=m, ncol=3)
  alpha   <- param[, 1]
  mu      <- param[, 2]
  sigma   <- param[, 3]
  
  a = list(alpha = alpha, mu = mu, sigma = sigma, gam = gam)
  
  a
  
}



#' @description Computes a_n based on empirical results found in Kasahara and Shimotsu (2015)
#' @export
#' @title anFormula
#' @name anFormula
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gamma_1, ..., gamma_m))
#' @param m The number of components in the mixture
#' @param n The number of observations
#' @param q The dimension of x (by default, 0)
#' @return a_n used to initialize values
#' @references Kasahara, H., and Shimotsu, K. (2015)
#' Testing the Number of Components in Normal Mixture Regression Models,
#' \emph{Journal of the American Statistical Association},
#' \bold{110}, 1632--1645.
anFormula <- function(parlist, m, n, t, q = 0)
  # Computes a_n for testing H_0 of m components
  # against H_1 of m+1 components
{
  
  if (q != 0){ # an when the dimension of X is not zero.
    return (switch(as.character(m), "1" =  0.1617017, "2" = 0.002509957, "3" =  0.05670169, "4" = 0.4858468, 0.5))
    #parlist$mu = parlist$mubeta[1,]
  }
  
    # return (0.5)
  if (m == 1) {
    b <- c(-0.61551296, 0.77642528, 28.14318354, -0.01554419)
    x <- (  b[1] + b[2]/t + b[3]/n ) / b[4] 
    an <- 1 / (1 + exp(x))
  }
  else if (m == 2) {
    omega <- omega.12(parlist)
    omega <- pmin(pmax(omega, 1e-16), 0.5-1e-16)  # an becomes NaN if omega[j]=0 or 1
    omega.term <- log(omega /(1-omega))
    
    b <-   c(-0.8112790,  -0.2882271,   4.6374028,  -0.1012959,  -0.1973225)
    x <- (  b[1] + b[2]/t + b[3]/n + b[5] * omega.term ) / b[4]   # maxa=1
    an <- 1 / (1 + exp(x))
  }
  else if (m == 3) {
    omega <- omega.123(parlist)
    omega <- pmin(pmax(omega, 1e-16), 1-1e-16)  # an becomes NaN if omega[j]=0 or 1
    omega.12 <- omega[1]
    omega.23 <- omega[2]
    omega.term <- log(omega.12 * omega.23 / ((1-omega.12)*(1-omega.23)))
    
    b <- c(-0.679611458, 0.611474005, 21.155661588, -0.110969483, 0.002174285)
    x <- (  b[1] + b[2]/t + b[3]/n + b[5] * omega.term ) / b[4]   # maxa=1
    an <- 1 / (1 + exp(x))
    
    # an <- 0.80 * x / (1 + x)
    #   x <- exp(-1.678 - 0.232 * log(t_omega) - 175.50/n)
    #   an <- 1.5 * x / (1 + x)
  } else if (m >= 4) {
    omega <- omega.1234(parlist)
    omega <- pmin(pmax(omega, 1e-16), 1-1e-16)  # an becomes NaN if omega[j]=0 or 1
    omega.12 <- omega[1]
    omega.23 <- omega[2]
    omega.34 <- omega[3]
    omega.term <- log(omega.12 * omega.23 * omega.34 / 
                        ((1-omega.12)*(1-omega.23)*(1-omega.34)))
    b <- c(-0.73506646,  0.25780605,  8.58533984, -0.12765448, -0.01269813)
    
    x <- (  b[1] + b[2]/t + b[3]/n + b[5] * omega.term ) / b[4]   # maxa=1
    an <- 1 / (1 + exp(x))
    
  }
  if (is.nan(an) | an ==0){ an = 1/n}
  return(an)
}  # end function anFormula



#' @description Computes the tuning parameter \eqn{a_n} based on
#' empirical formulas obtained by a similar method to Kasahara and Shimotsu (2015).
#' @export
#' @title anFormula
#' @name anFormula.t0
#' @param parlist parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gamma_1, ..., gamma_m)).
#' @param m number of components in the mixture.
#' @param n number of observations.
#' @param q dimension of x (default is 0).
#' @return tuning parameter \eqn{a_n}.
#' @references Kasahara, H., and Shimotsu, K. (2015)
#' Testing the Number of Components in Normal Mixture Regression Models,
#' \emph{Journal of the American Statistical Association},
#' \bold{110}, 1632--1645.
#' 
anFormula.t0 <- function(parlist, m, n, q = 0)
  # Computes a_n for testing H_0 of m components
  # against H_1 of m+1 components
{
  if (q != 0) # an when the dimension of X is not zero.
    return (switch(as.character(q), "1" = 0.3, "2" = 2.0, "3" = 2.4, "4" = 2.4, 0.3))
  
  if (m == 1) {
    an <- 0.30
  }
  else if (m == 2) {
    omega <- omega.12(parlist)
    omega <- pmin(pmax(omega, 1e-16), 1-1e-16)  # an becomes NaN if omega[j]=0 or 1
    omega.term <- log(omega /(1-omega)) 
    
    # coefficients of -(intercept, misclterm, nterm, -atermcoeff^2)/atermcoeff
    b <- c(-4.937477, -0.845460, -56.494216, -0.21091) 
    x <- exp(b[1] + b[2] * omega.term + b[3] / n - log(2) / b[4])  # maxa=1
    an <- 0.25 * x / (1 + x)
    #   x <- exp(-1.642 - 0.434 * log(omega / (1 - omega)) - 101.80/n)  # maxa=2
    #   an <- 1.8 * x / (1 + x)
  }
  else if (m == 3) {
    omega <- omega.123(parlist)
    omega <- pmin(pmax(omega, 1e-16), 1-1e-16)  # an becomes NaN if omega[j]=0 or 1
    omega.12 <- omega[1]
    omega.23 <- omega[2]
    omega.term <- log(omega.12 * omega.23 / ((1-omega.12)*(1-omega.23)))
    
    b <- c(-2.4481555, -0.2034425, -56.9171687, -0.27104) 
    x <- exp(b[1] + b[2] * omega.term + b[3] / n - log(2) / b[4])  # maxa=1
    an <- 0.25 * x / (1 + x)
    # an <- 0.80 * x / (1 + x)
    #   x <- exp(-1.678 - 0.232 * log(t_omega) - 175.50/n)
    #   an <- 1.5 * x / (1 + x)
  } else if (m == 4) {
    omega <- omega.1234(parlist)
    omega <- pmin(pmax(omega, 1e-16), 1-1e-16)  # an becomes NaN if omega[j]=0 or 1
    omega.12 <- omega[1]
    omega.23 <- omega[2]
    omega.34 <- omega[3]
    omega.term <- log(omega.12 * omega.23 * omega.34 / 
                        ((1-omega.12)*(1-omega.23)*(1-omega.34)))
    b <- c(-5.3663749, -0.2462147, -199.6375112, -0.300460) 
    x <- exp(b[1] + b[2] * omega.term + b[3] / n - log(2) / b[4])  # maxa=1
    an <- 0.25 * x / (1 + x)
  }
  else 
    an <- 1.0
  
  return (an)
}  # end function anFormula
