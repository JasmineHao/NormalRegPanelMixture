

#  get.t.from.lambda
# lambda Take in lambda, list of mu, sigma, beta
# @return the second order expansion
get.t.from.lambda <- function(lambda){
  mu_j <- lambda[1]
  sigma_j <- lambda[2]
  beta_j <- lambda[length(lambda)-2:length(lambda)]
  q_j <- length(lambda)
  q <- q_j - 2
  t <- rep(0,q_j * (q_j - 1) /2 + q_j)
  t[1] <- mu_j^2
  t[2] <- mu_j * (sigma_j^2)
  t[3] <- sigma_j^4
  if (q > 0){
    t[4:(3+q)] <- mu_j * beta_j
    t[(4+q):(3+2*q)] <- (sigma_j^2) * beta_j
    t[(4+2*q):(3+3*q)] <- beta_j^2
  }
  # :(3+3*q+(q-1)*q/2)
  if (q > 1) {
    qq <-  1
    for (j in 1:(q-1)) {
      for (i in (j+1):q) {
        t[(3+3*q+qq)] <- 2*beta_j[j]*beta_j[i]
        qq <- qq+1
      }
    }
  }
  return(t)
}

hermite <- function(Z, sigma)
# Computes the normalized Hermite polynomial for computing
# the critical values and p-values of the modified EM test
# Input
#   Z (n by m) : normalized data
#  sigma (m by 1): parameters
# Output
#   H (n by m by 4): normalized Hermite polynomials
{
n <- nrow(Z)
m    <- length(sigma)

H <- array(0,dim=c(n,m,4))
H[,,1] <- Z/sigma
H[,,2] <- t(t(Z^2-1)/2/sigma^2)
H[,,3] <- t(t(Z^3-3*Z)/6/sigma^3)
H[,,4] <- t(t(Z^4-6*Z^2+3)/24/sigma^4)
return(H)
}  # end function hermite


tKR <- function (a,b) {
  # Computes the transpose of the Khatri-Rao product of t(a) and t(b)
  n <- nrow(a)
  k <- ncol(a)
  KR <- matrix(unlist(lapply(1:k, function(i) (a[,i]*b))),nrow=n)
  KR
}


coef.to.list <- function(coefficients, z = NULL) {
# 　Convert coefficients to list
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

#' @description Generates mixed normal random variables with regressor x
#' @title rnormregpanelmix
#' @name rnormregpanelmix
#' @param n The number of observations
#' @param t Number of time periods
#' @param x n by k-1 matrix that does NOT include a constant
#' @param alpha m by 1 vector that represents proportions of components
#' @param mubeta k by m matrix that represents (mu times k regression coefficients) on x for m components
#' @param sigma m by 1 vector that represents sd of components
#' @return n by 1 vector that is formed by regressor x
rnormregpanelmix <- function (n, t, x = NULL, alpha, mubeta, sigma) {
  # Generates mixed normal random variables with regressor x
  # Input
  #  n : number of observations
  #   x : (n by k-1) matrix NOT including a constant
  #   alpha  : m-vector
  #  mubeta  : k by m matrix
  #  sigma  : m-vector
  # Output
  #  y : n by 1 vector
  # if (normalregMix.test.on) # initial values controlled by normalregMix.test.on
  #   set.seed(normalregMix.test.seed)
  
  m     <- length(alpha)
  nt    <- n*t
  mubeta   <- matrix(mubeta, ncol=m)
  y <- rep(0,nt)
  if (!is.null(x)){
    x <- as.matrix(x)
    if (nrow(x) != nt) { stop("y and x must have the same number of rows.") }
    x1   <- cbind(1,x)
    ii   <- sample(m, n, replace=TRUE, prob=alpha)
    
    for (nn in 1:n){
      y[(t*(nn-1)+1):(t*nn)]  <- rnorm(t, mean = x1[(t*(nn-1)+1):(t*nn),]%*% mubeta[, ii[nn]] , sd = sigma[ii[nn]])
    }
  } else {
    ii   <- sample(m, n, replace=TRUE, prob=alpha)
    for (nn in 1:n){
      y[(t*(nn-1)+1):(t*nn)]  <- rnorm(t, mean = mubeta[, ii[nn]] , sd = sigma[ii[nn]])
    }
  }
  
  y
  
} 



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

  if (sigma_i == sigma_j)
    if (mu_i > mu_j)
      omega_ji = pnorm((2 * log(alpha_j/alpha_i) - c)/(2*b),
                       mean = mu_i, sd = sigma_i)
  else
    omega_ji = 1 - pnorm((2 * log(alpha_j/alpha_i) - c)/(2*b),
                         mean = mu_i, sd = sigma_i)
  else {
    d <- 2 * log(alpha_j * sigma_i / (alpha_i * sigma_j)) - c + (b^2 / a)
    da <- max(d/a, 0)
    if (sigma_i > sigma_j)
      omega_ji = pnorm(sqrt(da)-b/a, mean = mu_i, sd = sigma_i) -
      pnorm(-sqrt(da)-b/a, mean = mu_i, sd = sigma_i)
    else
      omega_ji = 1 +
      pnorm(-sqrt(da)-b/a, mean = mu_i, sd = sigma_i) -
      pnorm(sqrt(da)-b/a, mean = mu_i, sd = sigma_i)
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




# Recursive function to compute Hermite polynomials
hermite_recursive <- function(n, x) {
  if (n == 0) {
    return(rep(1, length(x)))
  } else if (n == 1) {
    return(x)
  } else {
    H_n_minus_1 <- hermite_recursive(n - 1, x)
    H_n_minus_2 <- hermite_recursive(n - 2, x)
    return( x * H_n_minus_1 -  (n - 1) * H_n_minus_2)
  }
}


calculate_P_matrix <- function(data_c, T.even, T.odd, n.grid=2, type="indicator"){
  
  T <- nrow(data_c)
  N <- ncol(data_c)
  # Create a list to store the indicator matrices
  indicator_list <- list()
  
  if (type=="indicator"){
    for (t in 1:T){
      # Calculate the quantiles
      quantiles <- quantile(data_c[t,], probs = seq(0, 1, length.out = n.grid + 1))
      quantiles[1] <- -Inf
      quantiles[n.grid+1] <- Inf
      # Create the indicator matrix
      indicator_matrix <- matrix(0, nrow = N, ncol = n.grid)
      for (n in 1:length(data_c[t,])) {
        for (m in 1:n.grid) {
          if (data_c[t,n] > quantiles[m] & data_c[t,n] <= quantiles[m + 1]) {
            indicator_matrix[n, m] <- 1
          }
        }
      }
      # Add the indicator matrix to the list
      indicator_list[[t]] <- indicator_matrix
    }
  } else{
    for (t in 1:T){
      # Create the indicator matrix
      hermite_matrix <- sapply(0:(n.grid-1), function(k) hermite_recursive(k, data_c[t,]))
      indicator_list[[t]] <- hermite_matrix
    }
  }
  
  
  # Initialize the result matrix
  
  P_dim = length(Reduce(kronecker, lapply(T.even, function(t) indicator_list[[t]][1, ])))
  result_matrix_even <- matrix(0, nrow = N, ncol = P_dim)
  result_matrix_odd <- matrix(0, nrow = N, ncol = P_dim)
  
  # Iterate over the rows and compute the Kronecker product for each row
  for (n in 1:N) {
    kronecker_even_t <- Reduce(kronecker, lapply(T.even, function(t) indicator_list[[t]][n, ]))
    kronecker_odd_t <- Reduce(kronecker, lapply(T.odd, function(t) indicator_list[[t]][n, ]))
    
    # Add the Kronecker product to the result matrix
    result_matrix_even[n, ] <- kronecker_even_t
    result_matrix_odd[n, ] <- kronecker_odd_t
  }
  P_raw <- t(result_matrix_even) %*% result_matrix_odd
  # P_raw
}



calculate_W_P <- function(data,T.even, T.odd, n.grid=2, BB=199, type="indicator"){
  data_c <- data$Y
  n_size <- ncol(data_c)
  
  P_c <- calculate_P_matrix(data_c, T.even, T.odd, n.grid, type=type)
  
  
  ru <- matrix(runif(n_size * BB), nrow = n_size, ncol = BB)
  n_element <- length(as.vector(P_c))
  # Initialize the vec_P_b matrix
  vec_P_b <- matrix(0, nrow = BB, ncol = n_element)
  
  # Loop to generate vec_P_b
  for (i in 1:BB) {
    index <- ceiling(ru[, i] * n_size)
    data_b <- data$Y[,index]
    P_b <- calculate_P_matrix(data_b, T.even, T.odd, n.grid, type=type)
    vec_P_b[i, ] <- as.vector(P_b)
  }
  
  
  # Calculate mean_vec_P_b
  mean_vec_P_b <- as.vector(P_c) # colMeans(vec_P_b)
  
  # Initialize the W_b matrix
  W_b <- matrix(0, nrow =  n_element, ncol = n_element)
  
  # Fill the W_b matrix
  for (i in 1:n_element) {
    for (j in i:n_element) {
      if (i == j) {
        W_b[i, i] <- (1 / (BB - 1)) * sum((vec_P_b[, i] - mean_vec_P_b[i])^2)
      } else {
        W_b[i, j] <- (1 / (BB - 1)) * sum((vec_P_b[, i] - mean_vec_P_b[i]) * (vec_P_b[, j] - mean_vec_P_b[j]))
        W_b[j, i] <- W_b[i, j]
      }
    }
  }
  W_c = n_size*W_b 
  
  return(list(P_c = P_c, W_c = W_c))
}

construct_stat_KP <- function(P_c, W_c, r.test, n_size){
  # Perform SVD decomposition of the matrix
  P_svd <- svd(P_c)
  tol_s <- 1e-10
  # Extract the U, D, and V components of the SVD decomposition
  D <- P_svd$d
  U <- P_svd$u
  V <- P_svd$v
  
  U_12 <- U[1:(r.test),(r.test+1):ncol(U)]
  V_12 <- V[1:(r.test),(r.test+1):ncol(U)]
  
  if (r.test == 1){
    U_12 <- t(as.matrix(U_12))
    V_12 <- t(as.matrix(V_12))
  }
  
  U_22 <- U[(r.test+1):ncol(U),(r.test+1):ncol(U)]
  V_22 <- V[(r.test+1):ncol(V),(r.test+1):ncol(V)]
  

  svd_U22 <- svd(U_22)
  sqrtm_u_22 <- svd_U22$u %*% diag(svd_U22$d) %*% t(svd_U22$u)
  inv_u_22 <- svd_U22$u %*% diag(1/svd_U22$d) %*% t(svd_U22$v)
  
  svd_V22 <- svd(V_22)
  sqrtm_v_22 <- svd_V22$u %*% diag(svd_V22$d) %*% t(svd_V22$u)
  inv_v_22 <- svd_V22$u %*% diag(1/svd_V22$d) %*% t(svd_V22$v)
  
  
  # Construct the A_q_o and B_q_o matrix. 
  # A_q_o <- t( sqrtm_u_22 %*% solve(t(U_22)) %*% cbind(t(U_12), t(U_22)))
  A_q_o <- t( sqrtm_u_22 %*% inv_u_22 %*% cbind(t(U_12), t(U_22)))
  
  B_q_o <- sqrtm_v_22 %*% inv_v_22 %*% cbind(t(V_12), t(V_22))
  lambda_q <- t(A_q_o) %*% P_c %*% t(B_q_o)
  kron_BA_o <- kronecker(B_q_o, t(A_q_o))
  Omega_q <- kron_BA_o %*% W_c %*% t(kron_BA_o)
  
  if (qr(Omega_q)$rank == nrow(Omega_q)) {
    
    r <- nrow(Omega_q)
    rk_c <- n_size * sum(as.vector(lambda_q) * solve(Omega_q) %*% as.vector(lambda_q))
  } else {
    svd_result <- svd(Omega_q)
    s <- svd_result$d
    r <- sum(s > tol_s * max(s))
    inv_Omega_q <- svd_result$v[, 1:r] %*% diag(1 / s[1:r]) %*% t(svd_result$u[, 1:r])
    rk_c <- n_size * sum( t(as.vector(lambda_q)) %*% inv_Omega_q %*% as.vector(lambda_q))
  } 
  
  AIC_c = rk_c - 2*r
  BIC_c = rk_c - log(n_size)*r  
  HQ_c  = rk_c - 2*log(log(n_size))*r   
  return(rk_c)
}


return_p_val <- function(data_P_W, m, N){
  
  P_c <- data_P_W$P_c 
  W_c <- data_P_W$W_c
  s_1 <- dim(P_c)[1]
  s_2 <- dim(P_c)[2]
  
  rk_c <- construct_stat_KP(P_c, W_c, m, N)
  df <- (s_1 - m) * (s_2 - m)
  
  p.val <- 1 - pchisq(rk_c, df, lower.tail=TRUE)
  
  return(p.val) 
}
