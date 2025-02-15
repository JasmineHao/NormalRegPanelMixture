#' @description Computes the critical values of the modified EM test.
#' @export
#' @title normalpanelmixCrit
#' @name normalpanelmixCrit
#' @param y n by t vector of data
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gamma_1, ..., gamma_m))
#' @param z n by p matrix of regressor associated with gamma
#' @param values 3 by 1 Vector of length 3 (k = 1, 2, 3) at which the p-values are computed
#' @param nrep The number of replications used to compute p-values
#' @return A list with the following items:
#' \item{crit}{3 by 3 matrix of (0.1, 0.05, 0.01 critical values), jth row corresponding to k=j}
#' \item{pvals}{A vector of p-values at k = 1, 2, 3}
normalpanelmixCrit <- function(y, parlist, z = NULL, values = NULL, nrep = 10000)
{
  # if (normalregMix.test.on) # initial values controlled by normalregMix.test.on
  #   set.seed(normalregMix.test.seed)

  t <- nrow(y)
  n <- ncol(y)
  nt <- n*t
  y <- as.vector(y)
  p <- 0

  alpha <- parlist$alpha
  mu    <- parlist$mu
  sigma <- parlist$sigma
  gam <- parlist$gam
  m     <- length(alpha)

  if (!is.null(z)) {
    z <- as.matrix(z)
    p <- ncol(z)
    y   <- y - z %*% gam
  }

  pvals <- NULL

  if (m==1){
    crit <- qchisq(c(0.1,0.05,0.01), 2, lower.tail=F)
    if (!is.null(values))
    {
      k <- length(values)
      pvals <- pchisq(values, 2, lower.tail=F)
    }

  } else { # m>=2

    set.seed(123456)

    Z0 <- t((t(matrix(rep.int(y,m), ncol=m))-mu)/sigma)  # normalized data, nt by m
    f <- t(t(exp(-Z0^2/2)/sqrt(2*pi))/sigma)      # pdf, nt by m
    f0 <- colSums(t(f)*alpha)              # FIXIT(JASMINE), this is wrong for panel data data pdf, nt by 1

    S_alpha <- (f[, 1:(m-1)] - f[ ,m])/f0
    H <- hermite(Z0,sigma)

    S_mu    <- t(t(H[,,1]*f)*alpha)/f0    # n by m
    S_sigma <- t(t(H[,,2]*f)*alpha)/f0    # n by m
    S_gam <- rowSums(S_mu)*z

    S_lambda21 <- t(t(H[,,3]*f)*alpha)/f0  # n by m
    S_lambda22 <- t(t(H[,,4]*f)*alpha)/f0  # n by m

    S_lambda <- matrix(rbind(S_lambda21,S_lambda22),ncol=2*m)  # score wrt lambda (n by 2m)

    S_eta <- cbind(S_alpha,S_gam,S_mu,S_sigma)
    I11 <- t(S_eta) %*% S_eta/nt      # I_eta
    I21 <- t(S_lambda) %*% S_eta/nt    # I_{lambda eta}
    I22 <- t(S_lambda) %*% S_lambda/nt  # I_lambda, 2*m by 2*m

    Iall <- rbind(cbind(I11,t(I21)),cbind(I21,I22))
    if (rcond(Iall) < .Machine$double.eps)
    {
      eig <- eigen(Iall, symmetric=TRUE)
      tol2 <- (1e-14)*eig$values[1]
      vals <- eig$values
      vals[vals < tol2] <- tol2
      Iall.mod <- eig$vectors %*% (vals * t(eig$vectors))
      I11 <- Iall.mod[1:(p+3*m-1),1:(p+3*m-1)]
      I21 <- Iall.mod[(p+3*m):(p+5*m-1),1:(p+3*m-1)]
      I22 <- Iall.mod[(p+3*m):(p+5*m-1),(p+3*m):(p+5*m-1)]
    }

    I221 <- I22 - I21%*%solve(I11,t(I21))  # I_{lambda.eta}, which is also (var(W_lambda.eta))^{-1}

    e <- eigen(I221, symmetric=TRUE)    # eigenvalue decomposition is slower than chol but more stable
    u <- t(e$vec %*% (t(e$vec) * sqrt(e$val)) %*% matrix(rnorm(nrep*2*m), 2*m, nrep))

    EM <- matrix(0, nrow=nrep, ncol=m)

    for (jj in 1:m) {    # computes m quadratic forms
      uu <- u[,(2*jj-1):(2*jj)]
      Ivinv <- solve(I221[(2*jj-1):(2*jj),(2*jj-1):(2*jj)])
      EM[,jj] <- rowSums((uu%*%Ivinv)*uu)
    }

    max_EM <- apply(EM,1,max)  # max of m local modified EM stats
    max_EM_sort <- sort(max_EM)
    q <- ceiling(nrep*c(0.90,0.95,0.99))
    crit <- max_EM_sort[q]

    if (!is.null(values))
    {
      k <- length(values)
      pvals <- rowMeans(t(matrix(rep.int(max_EM_sort,k),ncol=k)) >= values)
    }

  }   # end else (m>=2)
  return(list(crit = crit, pvals = pvals))
}  # end function normalpanelmixCrit

#' @description Computes the bootstrap critical values of the modified EM test.
#' @export
#' @title normalpanelmixCritBoot
#' @name normalpanelmixCritBoot
#' @param y n by 1 vector of data
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gamma_1, ..., gamma_m))
#' @param z n by p matrix of regressor associated with gamma
#' @param values 3 by 1 Vector of length 3 (k = 1, 2, 3) at which the p-values are computed
#' @param ninits The number of initial candidates to be generated
#' @param nbtsp The number of bootstrap observations; by default, it is set to be 199.
#' @param parallel Determines whether package \code{doParallel} is used for calculation
#' @param cl Cluster used for parallelization (optional)
#' @return A list with the following items:
#' \item{crit}{3 by 3 matrix of (0.1, 0.05, 0.01 critical values), jth row corresponding to k=j}
#' \item{pvals}{A vector of p-values at k = 1, 2, 3}
normalpanelmixCritBoot <- function (y, parlist, z = NULL, values = NULL, ninits = 10,
                               nbtsp = 199, parallel = FALSE, cl = NULL, an = 0.5) {
  # if (normalregMix.test.on) # initial values controlled by normalregMix.test.on
  #   set.seed(normalregMix.test.seed)

  
  t <- nrow(y)
  n <- ncol(y)
  nt <- n*t
  y   <- as.vector(y) # TODO(Jasmine) Is this really needed?
  
  alpha <- parlist$alpha
  mu    <- parlist$mu
  sigma <- parlist$sigma
  gam <- parlist$gam
  m     <- length(alpha)
  an    <- anFormula(parlist = parlist, m = m, n = n)

  pvals <- NULL

  # Generate bootstrap observations
  ii    <- sample(m, nbtsp*nt, replace = TRUE, prob = alpha)
  ybset <- rnorm(nbtsp*nt, mean = mu[ii], sd = sigma[ii])
  ybset <- matrix(ybset, nrow = nt, ncol = nbtsp)

  if (!is.null(z)) {
    zgam <- as.matrix(z) %*% gam
    ybset <- ybset + replicate(nbtsp, as.vector(zgam))
  }
  if (parallel) { 
    if (is.null(cl))
      cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    out <- foreach (i.btsp = 1:nbtsp) %dopar% {
      
      normalpanelmixMEMtest (ybset[,i.btsp], t=t, m = m,
                        z = z, an = an, ninits = ninits, crit.method = "none") }
    on.exit(cl)
  }
  else
    out <- apply(ybset, 2, normalpanelmixMEMtest, t = t, m = m, z = z, an = an, ninits = ninits)

  emstat.b <- sapply(out, "[[", "emstat")  # 3 by nbstp matrix

  emstat.b <- t(apply(emstat.b, 1, sort))

  q <- ceiling(nbtsp*c(0.90,0.95,0.99))
  crit <- emstat.b[, q]

  if (!is.null(values)) { pvals <- rowMeans(emstat.b > values) }

  return(list(crit = crit, pvals = pvals))
}  # end function normalpanelmixCritBoot
