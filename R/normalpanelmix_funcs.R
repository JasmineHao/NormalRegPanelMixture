#' @description Given a pair of h and tau and data, compute ordinary &
#' penalized log-likelihood ratio resulting from MEM algorithm at k=1,2,3,
#' tailored for parallelization.
#' @export
#' @title normalpanelmixMaxPhiStep
#' @name normalpanelmixMaxPhiStep
#' @param htaupair A set of h and tau
#' @param y n by 1 vector of data
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gamma_1, ..., gamma_m))
#' @param z n by p matrix of regressor associated with gamma
#' @param p Dimension of z
#' @param an a term used for penalty function
#' @param ninits The number of randomly drawn initial values.
#' @param ninits.short The number of candidates used to generate an initial phi, in short MEM
#' @param epsilon.short The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param epsilon The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}.
#' @param maxit.short The maximum number of iterations in short EM.
#' @param maxit The maximum number of iterations.
#' @param verb Determines whether to print a message if an error occurs.
#' @return A list of phi, log-likelihood, and penalized log-likelihood resulting from MEM algorithm.
normalpanelmixMaxPhiStep <- function (htaupair, y, parlist, z = NULL, p,
                                 an,
                                 ninits, ninits.short,
                                 epsilon.short, epsilon,
                                 maxit.short, maxit,
                                 verb)
{
  alpha0 <- parlist$alpha

  m      <- length(alpha0)
  m1     <- m+1
  k      <- 1
  t  <- nrow(y)
  n  <- ncol(y)
  h      <- as.numeric(htaupair[1])
  tau    <- as.numeric(htaupair[2])

  mu0    <- parlist$mu
  mu0h   <- c(-1e+10,mu0,1e+10)        # m+2 by 1
  sigma0 <- parlist$sigma
  sigma0h<- c(sigma0[1:h],sigma0[h:m]) # m+1 by 1
  gam0 <- parlist$gam
  if (is.null(z)) {
    ztilde <- matrix(0) # dummy
    gam <- NULL
  }else{
    ztilde <- as.matrix(z)
  }
  # generate initial values
  tmp <- normalpanelmixPhiInit(y = y, parlist = parlist, z = z, h=h, tau = tau, ninits = ninits.short)

  # short EM
  b0 <- as.matrix(rbind(tmp$alpha, tmp$mu, tmp$sigma, tmp$gam))
  out.short <- cppnormalpanelmixPMLE(b0, as.vector(y), ztilde, mu0h, sigma0h, m1, p, t, an, maxit.short, ninits.short,
                                epsilon.short, tau, h, k)
  components <- order(out.short$penloglikset, decreasing = TRUE)[1:ninits]
  if (verb && any(out.short$notcg)) {
    cat(sprintf("non-convergence rate at short-EM = %.3f\n",mean(out.short$notcg)))
  }
  # long EM
  b1 <- as.matrix(b0[ ,components])
  
  out <- cppnormalpanelmixPMLE(b1, as.vector(y), ztilde, mu0h, sigma0h, m1, p, t, an, maxit, ninits, epsilon, tau, h, k)
  
  index     <- which.max(out$penloglikset)
  alpha <- b1[1:m1,index]
  mu <- b1[(1+m1):(2*m1),index]
  sigma <- b1[(1+2*m1):(3*m1),index]
  if (!is.null(z)) {
    gam     <- b1[(3*m1+1):(3*m1+p),index]
  }
  mu.order  <- order(mu)
  alpha     <- alpha[mu.order]
  mu        <- mu[mu.order]
  sigma     <- sigma[mu.order]
  sigma0h <- sigma0h[mu.order]
  b <- as.matrix( c(alpha, mu, sigma, gam) )

  # initilization
  loglik <-  vector("double", 3)
  penloglik <-  vector("double", 3)
  coefficient <- vector("double", length(b))

  penloglik[1] <- out$penloglikset[[index]]
  loglik[1]    <- out$loglikset[[index]]
  for (k in 2:3) {
    ninits <- 1
    maxit <- 1 # One iteration ahead is enough
    # Two EM steps
    out <- cppnormalpanelmixPMLE(b, as.vector(y), ztilde, mu0h, sigma0h, m1, p,t, an, maxit, ninits, epsilon, tau, h, k)
    alpha <- b[1:m1,1] # b has been updated
    
    tau <- alpha[h] / (alpha[h] + alpha[h+1]) 
    
    mu <- b[(1+m1):(2*m1),1]
    sigma <- b[(1+2*m1):(3*m1),1]
    if (!is.null(z)) {
      gam     <- b[(3*m1+1):(3*m1+p),1]
    }
    loglik[k]    <- out$loglikset[[1]]
    penloglik[k]   <- out$penloglikset[[1]]

    # Check singularity: if singular, break from the loop
    if ( any(sigma < 1e-06) || any(alpha < 1e-06) || is.na(sum(alpha)) ) {
      loglik[k]    <- -Inf
      penloglik[k]   <- -Inf
      break
    }

    mu.order  <- order(mu)
    alpha     <- alpha[mu.order]
    mu        <- mu[mu.order]
    sigma     <- sigma[mu.order]
    sigma0h <- sigma0h[mu.order]
  }
  coefficient <- as.matrix( c(alpha, mu, sigma, gam) ) # at k=3
  parlist <- list(alpha = alpha, mubeta = mu, sigma = sigma, gam = gam)
  return (list(coefficient = coefficient, loglik = loglik, penloglik = penloglik, parlist = parlist))
}




#' @description Compute ordinary & penalized log-likelihood ratio resulting from
#' MEM algorithm at k=1,2,3.
#' @export
#' @title normalpanelmixMaxPhi
#' @name normalpanelmixMaxPhi
#' @param y n by 1 vector of data
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gamma_1, ..., gamma_m))
#' @param z n by p matrix of regressor associated with gamma
#' @param an a term used for penalty function
#' @param tauset A set of initial tau value candidates
#' @param ninits The number of randomly drawn initial values.
#' @param epsilon.short The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param epsilon The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}.
#' @param maxit.short The maximum number of iterations in short EM.
#' @param maxit The maximum number of iterations.
#' @param verb Determines whether to print a message if an error occurs.
#' @param parallel Determines whether package \code{doParallel} is used for calculation
#' @param cl Cluster used for parallelization; if it is \code{NULL}, the system will automatically
#' create a new one for computation accordingly.
#' @return A list with items:
#' \item{loglik}{Log-likelihood resulting from MEM algorithm at k=1,2,3.}
#' \item{penloglik}{Penalized log-likelihood resulting from MEM algorithm at k=1,2,3.}
normalpanelmixMaxPhi <- function (y, parlist, z = NULL, an, tauset = c(0.1,0.3,0.5), ninits = 10, epsilon.short = 1e-02, epsilon = 1e-08,
                             maxit.short = 500, maxit = 2000,
                             verb = FALSE,
                             parallel = FALSE,
                             cl = NULL, k_max = 3) {
  # Given a parameter estimate of an m component model and tuning paramter an,
  # maximize the objective function for computing the modified EM test statistic
  # for testing H_0 of m components against H_1 of m+1 for a univariate normal finite mixture

  warn  <- options(warn=-1) # Turn off warnings

  if (!is.null(z)) {
    z     <- as.matrix(z)
    p     <- ncol(z)
  } else {
    p <- 0
  }
  m <- length(parlist$alpha)

  ninits.short <- ninits*10*(1+p)*m

  loglik.all <- matrix(0,nrow=m*length(tauset),ncol=k_max)
  penloglik.all <- matrix(0,nrow=m*length(tauset),ncol=k_max)
  coefficient.all <- matrix(0,nrow=m*length(tauset),ncol=(3*(m+1)+p))

  if (parallel) {
    if (is.null(cl)){
      cl <- makeCluster(detectCores())
      # print(detectCores())
    }
    registerDoParallel(cl)
    results <- foreach (t = 1:length(tauset),
                        .export = 'normalpanelmixMaxPhiStep', .combine = c)  %:%
      foreach (h = 1:m) %dopar% {
        normalpanelmixMaxPhiStep (c(h, tauset[t]), y, parlist, z, p,
                             an,
                             ninits, ninits.short,
                             epsilon.short, epsilon,
                             maxit.short, maxit,
                             verb, k_max = k_max) }
    on.exit(cl)
    loglik.all <- t(sapply(results, "[[", "loglik"))
    penloglik.all <- t(sapply(results, "[[", "penloglik"))
    coefficient.all <- t(sapply(results, "[[", "coefficient"))
  }
  else{
    for (h in 1:m){
      for (t in 1:length(tauset)) {
        rowindex <- (t-1)*m + h
        tau <- tauset[t]
        result <- normalpanelmixMaxPhiStep (c(h, tauset[t]), y, parlist, z, p,
                                            an,
                                            ninits, ninits.short,
                                            epsilon.short, epsilon,
                                            maxit.short, maxit,
                                            verb, k_max = k_max)
        loglik.all[rowindex,] <- result$loglik
        penloglik.all[rowindex,] <- result$penloglik
        coefficient.all[rowindex,] <- result$coefficient
      }
    }
  }
  # loglik <- apply(loglik.all, 2, max)  # 3 by 1 vector
  # penloglik <- apply(penloglik.all, 2, max)  # 3 by 1 vector
  index <- which.max(loglik.all[ ,k_max]) # a par (h,m) that gives the highest likelihood at k=3
  coefficient <- as.vector(coefficient.all[index,])
  loglik <- loglik.all[index,k_max]
  penloglik <- penloglik.all[index,k_max]
  out <- list(coefficient = coefficient, loglik = loglik, penloglik = penloglik, parlist=result$parlist)

  out


}  # end normalpanelmixMaxPhi


#' @description Given a pair of h and tau and data, compute ordinary &
#' penalized log-likelihood ratio resulting from MEM algorithm at k=1,2,3,
#' tailored for parallelization.
#' @export
#' @title normalpanelmixMaxPhiStep
#' @name normalpanelmixMaxPhiStep
#' @param htaupair A set of h and tau
#' @param y n by 1 vector of data
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gamma_1, ..., gamma_m))
#' @param z n by p matrix of regressor associated with gamma
#' @param p Dimension of z
#' @param an a term used for penalty function
#' @param ninits The number of randomly drawn initial values.
#' @param ninits.short The number of candidates used to generate an initial phi, in short MEM
#' @param epsilon.short The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param epsilon The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}.
#' @param maxit.short The maximum number of iterations in short EM.
#' @param maxit The maximum number of iterations.
#' @param verb Determines whether to print a message if an error occurs.
#' @return A list of phi, log-likelihood, and penalized log-likelihood resulting from MEM algorithm.
normalpanelmixMaxPhiStep <- function (htaupair, y, parlist, z = NULL, p,
                                 an,
                                 ninits, ninits.short,
                                 epsilon.short, epsilon,
                                 maxit.short, maxit,
                                 verb, k_max = 3)
{
  alpha0 <- parlist$alpha

  m      <- length(alpha0)
  m1     <- m+1
  k      <- 1
  t  <- nrow(y)
  n  <- ncol(y)
  h      <- as.numeric(htaupair[1])
  tau    <- as.numeric(htaupair[2])

  mu0    <- parlist$mu
  mu0h   <- c(-1e+10,mu0,1e+10)        # m+2 by 1
  sigma0 <- parlist$sigma
  sigma0h<- c(sigma0[1:h],sigma0[h:m]) # m+1 by 1
  gam0 <- parlist$gam
  if (is.null(z)) {
    ztilde <- matrix(0) # dummy
    gam <- NULL
  }else{
    ztilde <- as.matrix(z)
  }
  # generate initial values
  tmp <- normalpanelmixPhiInit(y = y, parlist = parlist, z = z, h=h, tau = tau, ninits = ninits.short)

  # short EM
  b0 <- as.matrix(rbind(tmp$alpha, tmp$mu, tmp$sigma, tmp$gam))
  out.short <- cppnormalpanelmixPMLE(b0, as.vector(y), ztilde, mu0h, sigma0h, m1, p, t, an, maxit.short, ninits.short,
                                epsilon.short, tau, h, k)
  components <- order(out.short$penloglikset, decreasing = TRUE)[1:ninits]
  if (verb && any(out.short$notcg)) {
    cat(sprintf("non-convergence rate at short-EM = %.3f\n",mean(out.short$notcg)))
  }
  # long EM
  b1 <- as.matrix(b0[ ,components])
  
  out <- cppnormalpanelmixPMLE(b1, as.vector(y), ztilde, mu0h, sigma0h, m1, p, t, an, maxit, ninits, epsilon, tau, h, k)
  
  index     <- which.max(out$penloglikset)
  alpha <- b1[1:m1,index]
  mu <- b1[(1+m1):(2*m1),index]
  sigma <- b1[(1+2*m1):(3*m1),index]
  if (!is.null(z)) {
    gam     <- b1[(3*m1+1):(3*m1+p),index]
  }
  mu.order  <- order(mu)
  alpha     <- alpha[mu.order]
  mu        <- mu[mu.order]
  sigma     <- sigma[mu.order]
  sigma0h <- sigma0h[mu.order]
  b <- as.matrix( c(alpha, mu, sigma, gam) )

  # initilization
  loglik <-  vector("double", k_max)
  penloglik <-  vector("double", k_max)
  coefficient <- vector("double", length(b))

  penloglik[1] <- out$penloglikset[[index]]
  loglik[1]    <- out$loglikset[[index]]
  for (k in 2:k_max) {
    ninits <- 1
    maxit <- 1
    # Two EM steps
    out <- cppnormalpanelmixPMLE(b, as.vector(y), ztilde, mu0h, sigma0h, m1, p,t, an, maxit, ninits, epsilon, tau, h, k)
    alpha <- b[1:m1,1] # b has been updated
    tau <- alpha[h] / (alpha[h] + alpha[h + 1])
    mu <- b[(1+m1):(2*m1),1]
    sigma <- b[(1+2*m1):(3*m1),1]
    if (!is.null(z)) {
      gam     <- b[(3*m1+1):(3*m1+p),1]
    }
    loglik[k]    <- out$loglikset[[1]]
    penloglik[k]   <- out$penloglikset[[1]]

    # Check singularity: if singular, break from the loop
    if ( any(sigma < 1e-06) || any(alpha < 1e-06) || is.na(sum(alpha)) ) {
      loglik[k]    <- -Inf
      penloglik[k]   <- -Inf
      break
    }

    mu.order  <- order(mu)
    alpha     <- alpha[mu.order]
    mu        <- mu[mu.order]
    sigma     <- sigma[mu.order]
    sigma0h <- sigma0h[mu.order]
  }
  coefficient <- as.matrix( c(alpha, mu, sigma, gam) ) # at k=3
  parlist <- list(alpha = alpha, mubeta = mu, sigma = sigma, gam = gam)
  return (list(coefficient = coefficient, loglik = loglik, penloglik = penloglik, parlist = parlist))
}

#' @description  Sequentially performs MEM test given the data for y and x
#' on the null hypothesis H_0: m = m_0 where m_0 is in {1, 2, ..., maxm}
#' @export
#' @title normalpanelmixMEMtestSeq
#' @name normalpanelmixMEMtestSeq
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x (if exists)
#' @param z n by p matrix of regressor associated with gamma
#' @param maxm The maximum number of components set as null hypothesis in the mixture
#' @param ninits The number of randomly drawn initial values.
#' @param maxit The maximum number of iterations.
#' @param nbtsp The number of bootstrap observations; by default, it is set to be 199
#' @param parallel Determines whether package \code{doParallel} is used for calculation
#' @param cl Cluster used for parallelization; if it is \code{NULL}, the system will automatically
#' create a new one for computation accordingly.
#' @param crit.bootstrap.from The minimum m in null hypothesis to have critical values
#' calculated from bootstrap for the test statistics
#' @return A list of with the following items:
#' \item{alpha}{maxm by maxm matrix, whose i-th column is a vector of alphas estimated given the null hypothesis m_0 = i}
#' \item{mu}{maxm by maxm matrix, whose i-th column is a vector of mus estimated given the null hypothesis m_0 = i}
#' \item{sigma}{maxm by maxm matrix, whose i-th column is a vector of sigmas estimated given the null hypothesis m_0 = i}
#' \item{beta}{A list of length maxm, whose i-th element is a q times i matrix of betas estimated given the null hypothesis m_0 = i}
#' \item{gam}{maxm by maxm matrix, whose i-th column is a vector of gammas estimated given the null hypothesis m_0 = i}
#' \item{emstat}{A maxm vector of values of modified EM statistics of the model at m_0 = 1, 2, ..., maxm}
#' \item{pvals}{A maxm by 3 matrix whose i-th row indicates a vector of p-values at k = 1, 2, 3}
#' \item{aic}{A maxm vector of Akaike Information Criterion of the fitted model at m_0 = 1, 2, ..., maxm}
#' \item{bic}{A maxm vector of Bayesian Information Criterion of the fitted model at m_0 = 1, 2, ..., maxm}
#' \item{loglik}{A maxm vector of log-likelihood values of the model at m_0 = 1, 2, ..., maxm}
#' \item{penloglik}{A maxm vector of penalized log-likelihood values of the model at m_0 = 1, 2, ..., maxm}
#' \item{pmle.result}{A list of output from normalpanelmixPMLE under the number of components selected by sequantial hypothesis testing}
#' @examples
#' data(faithful)
#' attach(faithful)
#' normalpanelmixMEMtestSeq(y = eruptions)
normalpanelmixMEMtestSeq <- function (y, x = NULL, z = NULL,  maxm = 3, ninits = 10, maxit = 2000,
                                 nbtsp = 199, parallel = FALSE, cl = NULL,
                                 crit.bootstrap.from = 3) {
  # Compute the modified EM test statistic for testing H_0 of m components
  # against H_1 of m+1 components for a univariate finite mixture of normals


  if (!is.null(x))
    return (regmixMEMtestSeq(y = y, x = x, z = z, maxm = maxm, ninits = ninits, maxit = maxit,
                             nbtsp = nbtsp, parallel = parallel, cl = cl,
                             crit.bootstrap.from = crit.bootstrap.from))

  t <- dim(y)[1] #Number of year
  n <- dim(y)[2] #Number of firms
  # y   <- as.vector(y)
  p   <- 0
  if (!is.null(z)) {
    z <- as.matrix(z)
    p <- ncol(z)
    gam <- matrix(0, nrow = p, ncol = maxm)
  }
  else
    gam <- NULL

  out   <- vector('list', length = maxm)
  aic    <- bic <- double(maxm)
  pvals   <- emstat <- matrix(0, nrow = maxm, ncol = 3)
  loglik  <- penloglik <- double(maxm)

  alpha   <- mu <- sigma <- matrix(0, nrow = maxm, ncol = maxm)

  # Test H_0:m=1, H_0:m=2, ...
  binit <- NULL
  for (m in 1:maxm){
    pmle.result   <- normalpanelmixPMLE(y = y, m = m, z = z, vcov.method = "none",
                                   ninits = ninits, maxit = maxit, binit = binit)
    loglik[m] <- loglik0 <- pmle.result$loglik
    penloglik[m] <- penloglik0 <- pmle.result$penloglik
    aic[m]  <- pmle.result$aic
    bic[m]  <- pmle.result$bic

    parlist <- pmle.result$parlist
    alpha0  <- parlist$alpha
    mu0     <- parlist$mu
    sigma0  <- parlist$sigma
    gam0  <- parlist$gam

    alpha[,m] <- c(alpha0, double(maxm - m))
    mu[,m]     <- c(mu0, double(maxm - m))
    sigma[,m] <- c(sigma0, double(maxm - m))

    cat(sprintf("%d-component model estimate:\n",m))
    tab = as.matrix(rbind(alpha0, mu0, sigma0))
    rownames(tab) <- c("alpha", "mu", "sigma")
    colnames(tab) <- c(paste("comp", ".", 1:m, sep = ""))
    print(tab, digits = 4)

    if (!is.null(z)){
      gam[, m] <- gam0
      cat("gam =", gam0,"\n")
    }
    cat(sprintf("\nAIC, BIC, and log-likelihood of 1 to %.i", m), "component models \n")
    cat(c("AIC    =", sprintf(' %.2f', aic[1:m])), "\n")
    cat(c("BIC    =", sprintf(' %.2f', bic[1:m])), "\n")
    cat(c("loglik =", sprintf('%.2f', loglik[1:m])), "\n\n")

    if (m <= maxm){

      cat(sprintf("Testing the null hypothesis of %d components\n", m))

      an    <- anFormula(parlist = parlist, m = m, n = n, t = t)
      par1  <- normalpanelmixMaxPhi(y = y, parlist = parlist, z = z, an = an,
                               ninits = ninits, maxit = maxit, parallel = parallel)
      emstat.m  <- 2*(par1$penloglik - loglik0)

      # use the estimate of b as one of the initial values
      binit <- par1$coefficient

      cat(c("modified EM-test statitic ", sprintf('%.3f',emstat.m)),"\n")
      if (m <= crit.bootstrap.from) {
        em.out <- normalpanelmixCrit(y=y, parlist=parlist, z=z, values = emstat.m)
        cat(c("asymptotic p-value       ", sprintf('%.3f',em.out$pvals)),"\n \n")
      } else {
        em.out <- normalpanelmixCritBoot(y=y, parlist=parlist, z=z, values = emstat.m,
                                    ninits = ninits, nbtsp = nbtsp, parallel = parallel, cl = cl)
        cat(c("bootstrap p-value        ", sprintf('%.3f',em.out$pvals)),"\n \n")
      }
      pvals[m,]     <- em.out$pvals
      emstat[m,]    <- emstat.m
    }
  }

  for (m in 1:maxm)
    if ( pvals[m,2] >= 0.05 ) {
      cat(sprintf("\nThe number of components selected by Sequential Hypothesis Testing (alpha=0.05) = %.i", m), " \n")
      cat(sprintf("The number of components selected by AIC = %.i", which.min(aic)), " \n")
      cat(sprintf("The number of components selected by BIC = %.i", which.min(bic)), " \n")
      binit <- as.vector(c(alpha[1:m,m], mu[1:m,m], sigma[1:m,m],  gam[,m]))
      pmle.result   <- normalpanelmixPMLE(y = y, m = m, z = z, vcov.method = "none", ninits = 2, maxit = maxit, binit = binit)
      cat(sprintf("\nThe summary of the estimated %.i", m), "component model: \n")
      print(summary(pmle.result))
      break
    }

  a = list(alpha = alpha, mu = mu, sigma = sigma, gam = gam, emstat = emstat, pvals = pvals, aic = aic, bic = bic, loglik = loglik, penloglik = penloglik, pmle.result = pmle.result)

  a
}  # end normalpanelmixMEMtestSeq


#' @description Generates lists of parameters for initial candidates used by
#' the modified EM test for mixture of normals.
#' @title normalpanelmixPhiInit
#' @name normalpanelmixPhiInit
#' @param y n by 1 vector of data
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma,
#' and gamma in the form of (alpha = (alpha_1, ..., alpha_m),
#' mu = (mu_1, ..., mu_m), sigma = (sigma_1, ..., sigma_m),
#' gam = (gamma_1, ..., gamma_m))
#' @param z n by p matrix of regressor associated with gamma
#' @param h h used as index for pivoting
#' @param tau Tau used to split the h-th component
#' @param ninits number of initial values to be generated
#' @return A list with the following items:
#' \item{alpha}{m+1 by ninits matrix for alpha}
#' \item{mu}{m+1 by ninits matrix for mu}
#' \item{sigma}{m+1 by ninits matrix for sigma}
#' \item{gam}{m+1 by ninits matrix for gamma}
normalpanelmixPhiInit <- function(y, parlist, z = NULL, h, tau, ninits = 1) {
  # if (normalregMix.test.on) # initial values controlled by normalregMix.test.on
  #   set.seed(normalregMix.test.seed)

  t <- nrow(y)
  n <- ncol(y)
  p <- ncol(z)

  mu0 <- parlist$mu
  sigma0 <- parlist$sigma
  alpha0 <- parlist$alpha
  m <- length(alpha0)

  if (is.null(z)) {
    gam <- NULL
  } else {
    gam0 <- parlist$gam
    y <- y - z %*% gam0
    gam <- matrix(runif(p * ninits, min = 0.5, max = 1.5), nrow = p) * gam0
  }

  if (m >= 2) {
    mid <- (mu0[1:(m - 1)] + mu0[2:m]) / 2 # m-1 by 1
    lb0 <- c(min(y), mid) # m by 1
    lb <- c(lb0[1:h], lb0[h:m]) # m+1 by 1
    ub0 <- c(mid, max(y)) # m by 1
    ub <- c(ub0[1:h], ub0[h:m]) # m+1 by 1
  } else {
    lb <- c(min(y), min(y))
    ub <- c(max(y), max(y))
  }

  mu <- matrix(runif((m + 1) * ninits, min = lb, max = ub), nrow = m + 1)

  sigma.hyp <- c(sigma0[1:h], sigma0[h:m]) # m+1 by 1
  sigma <- matrix(runif((m + 1) * ninits, min = sigma.hyp * 0.25, max = sigma.hyp * 2), nrow = m + 1)

  alpha.hyp <- c(alpha0[1:h], alpha0[h:m]) # m+1 by 1
  alpha.hyp[h:(h + 1)] <- c(alpha.hyp[h] * tau, alpha.hyp[h + 1] * (1 - tau))
  alpha <- matrix(rep.int(alpha.hyp, ninits), nrow = m + 1)

  list(alpha = alpha, mu = mu, sigma = sigma, gam = gam)
} # end function normalpanelmixPhiInit


normalpanelmixPMLE.M1 <- function(y, parlist, x = NULL, z = NULL,an=NULL, vcov.method = c("Hessian", "OPG", "none"), ninits = 25, epsilon = 1e-08, maxit = 2000, epsilon.short = 1e-02, maxit.short = 500, binit = NULL){
  t <- dim(y)[1] # Number of year
  n <- dim(y)[2] # Number of firms
  nt <- n * t
  p <- 0
  ninits.short <- ninits * 10 * (1 + p) 

  m <- length(parlist$alpha)
  m1 <- m + 1
  if (!is.null(x)) {
    return(regpanelmixPMLE(
      y = y, x = x, m = m, z = z, vcov.method = vcov.method,
      ninits = ninits, epsilon = epsilon, maxit = maxit,
      epsilon.short = epsilon.short, maxit.short = maxit.short,
      binit = binit
    ))
  }
  if (!is.null(z)) {
    z <- as.matrix(z)
    p <- ncol(z)
    if (nrow(z) != nt) {
      stop("y and z must have the same number of rows.")
    }
    ls.out <- lsfit(z, y)
    sd0 <- sqrt(mean(ls.out$residuals^2))
  } else {
    sd0 <- sd(y) * sqrt((n - 1) / n)
  }
  
  loglik.all <- matrix(0, nrow = m, ncol = 1)
  penloglik.all <- matrix(0, nrow = m, ncol = 1)
  coefficient.all <- matrix(0, nrow = m, ncol = (3 * (m + 1) + p))
  if (is.null(an)) {
    an <- 1 / n
  }
  for (h in 1:m) {
    mu0 <- parlist$mu
    mu0h <- c(-1e+10, mu0, 1e+10) # m+2 by 1
    sigma0 <- parlist$sigma
    sigma0h <- c(sigma0[1:h], sigma0[h:m]) # m+1 by 1
    #sigma0h <- rep(sd0, m)
    gam0 <- parlist$gam
    tau <- 0.5
    if (is.null(z)) {
      ztilde <- matrix(0) # dummy
      gam <- NULL
    } else {
      ztilde <- as.matrix(z)
    }
    # generate initial values
    k <- 0
    psih <- 1 # use hard bounds for estimation
    
    tmp <- normalpanelmixPhiInit(y = as.vector(y), parlist = parlist, z = z, h = h, tau = tau, ninits = ninits.short)
    
    b0 <- as.matrix(rbind(tmp$alpha, tmp$mu, tmp$sigma, tmp$gam))
    out.short <- cppnormalpanelmixPMLE( b0, as.vector(y), ztilde, mu0h, sigma0h, m1, p, t, an, maxit.short, ninits.short, epsilon.short, 0.5, h, k, psih)
    components <- order(out.short$penloglikset, decreasing = TRUE)[1:ninits]
    b1 <- b0[, components] # b0 has been updated
    
    out <- cppnormalpanelmixPMLE( b1, as.vector(y), ztilde, mu0h, sigma0h, m1, p, t, an, maxit, ninits, epsilon, 0.5, h, k, psih)
    index <- which.max(out$penloglikset)
    alpha <- b1[1:m1, index] # b0 has been updated
    mu <- b1[(1 + m1):(2 * m1), index]
    sigma <- b1[(1 + 2 * m1):(3 * m1), index]
    if (!is.null(z)) {
      gam <- b1[(3 * m1 + 1):(3 * m1 + p), index]
    }
    penloglik <- out$penloglikset[index]
    loglik <- out$loglikset[index]
    parlist.M1 <- list(alpha = alpha, mubeta = mu, sigma = sigma, gam = gam)
    coefficients <- unlist(parlist.M1)

    loglik.all[h, ] <- loglik
    penloglik.all[h, ] <- penloglik
    coefficient.all[h, ] <- coefficients
  }
  # cppnormalpanelmixPMLE(b1, as.vector(y), ztilde, mu0h, sigma0h, m1, p, t, an, maxit, ninits, epsilon, tau, h, k)
  index <- which.max(penloglik.all[, 1])
  coefficient <- as.vector(coefficient.all[index, ])
  loglik <- loglik.all[index, 1]
  penloglik <- penloglik.all[index, 1]
  out <- list(coefficient = coefficient, loglik = loglik, penloglik = penloglik)
  out
}

#' @description Estimates parameters of a finite panel mixture of univariate normals by
#' penalized maximum log-likelhood functions.
#' @export
#' @title normalpanelmixPMLE
#' @name normalpanelmixPMLE
#' @param y n by 1 vector of data
#' @param x n by q matrix of data for x (if exists)
#' @param m The number of components in the mixture
#' @param z n by p matrix of regressor associated with gamma
#' @param vcov.method Method used to compute the variance-covariance matrix, one of \code{"Hessian"} and \code{"OPG"}.
#' The default option is \code{"Hessian"}. When \code{method = "Hessian"}, the variance-covarince matrix is
#' estimated by the Hessian using the formula given in Boldea and Magnus (2009).
#' When \code{method = "OPG"}, the outer product of gradients is used.
#' @param ninits The number of randomly drawn initial values.
#' @param epsilon The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}.
#' @param maxit The maximum number of iterations.
#' @param epsilon.short The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param maxit.short The maximum number of iterations in short EM.
#' @param binit The initial value of parameter vector that is included as a candidate parameter vector
#' @return  A list of class \code{normalpanelmix} with items:
#' \item{coefficients}{A vector of parameter estimates. Ordered as \eqn{\alpha_1,\ldots,\alpha_m,\mu_1,\ldots,\mu_m,\sigma_1,\ldots,\sigma_m,\gam}.}
#' \item{parlist}{The parameter estimates as a list containing alpha, mu, and sigma (and gam if z is included in the model).}
#' \item{vcov}{The estimated variance-covariance matrix.}
#' \item{loglik}{The maximized value of the log-likelihood.}
#' \item{penloglik}{The maximized value of the penalized log-likelihood.}
#' \item{aic}{Akaike Information Criterion of the fitted model.}
#' \item{bic}{Bayesian Information Criterion of the fitted model.}
#' \item{postprobs}{n by m matrix of posterior probabilities for observations}
#' \item{components}{n by 1 vector of integers that indicates the indices of components
#' each observation belongs to based on computed posterior probabilities}
#' \item{call}{The matched call.}
#' \item{m}{The number of components in the mixture.}
#' @note \code{normalpanelmixPMLE} maximizes the penalized log-likelihood function
#' using the EM algorithm with combining short and long runs of EM steps as in Biernacki et al. (2003).
#' \code{normalpanelmixPMLE} first runs the EM algorithm from \code{ninits}\eqn{* 4m(1 + p)} initial values
#' with the convertence criterion \code{epsilon.short} and \code{maxit.short}.
#' Then, \code{normalpanelmixPMLE} uses \code{ninits} best initial values to run the EM algorithm
#' with the convertence criterion \code{epsilon} and \code{maxit}.
#' @references     Biernacki, C., Celeux, G. and Govaert, G. (2003)
#' Choosing Starting Values for the EM Algorithm for Getting the
#' Highest Likelihood in Multivariate Gaussian Mixture Models,
#' \emph{Computational Statistics and Data Analysis}, \bold{41}, 561--575.
#'
#' Boldea, O. and Magnus, J. R. (2009)
#' Maximum Likelihood Estimation of the Multivariate Normal Mixture Model,
#' \emph{Journal of the American Statistical Association},
#' \bold{104}, 1539--1549.
#'
#' Chen, J., Tan, X. and Zhang, R. (2008)
#' Inference for Normal Mixtures in Mean and Variance,
#' \emph{Statistica Sinica}, \bold{18}, 443--465.
#'
#' McLachlan, G. J. and Peel, D. (2000) \emph{Finite Mixture Models}, John Wiley \& Sons, Inc.
#' @examples
#' data(faithful)
#' attach(faithful)
#'
#' normalpanelmixPMLE(y = eruptions, m = 1)
#' normalpanelmixPMLE(y = eruptions, m = 2)
#'
#' out <- normalpanelmixPMLE(y = eruptions, m = 2)
#' summary(out)
normalpanelmixPMLE <- function(y, x = NULL, m = 2, z = NULL, vcov.method = c("Hessian", "OPG", "none"),
                               ninits = 25, epsilon = 1e-08, maxit = 2000,
                               epsilon.short = 1e-02, maxit.short = 500, binit = NULL, in.coefficient = NULL) {
  t <- dim(y)[1] # Number of year
  n <- dim(y)[2] # Number of firms
  nt <- n * t
  p <- 0
  y <- as.vector(y)
  if (!is.null(x)) {
    return(regpanelmixPMLE(
      y = y, x = x, m = m, z = z, vcov.method = vcov.method,
      ninits = ninits, epsilon = epsilon, maxit = maxit,
      epsilon.short = epsilon.short, maxit.short = maxit.short,
      binit = binit
    ))
  }

  gam <- NULL
  ninits.short <- ninits * 10 * (1 + p) * m
  vcov.method <- match.arg(vcov.method)

  if (!is.null(z)) {
    z <- as.matrix(z)
    p <- ncol(z)
    if (nrow(z) != nt) {
      stop("y and z must have the same number of rows.")
    }
    ls.out <- lsfit(z, y)
    sd0 <- sqrt(mean(ls.out$residuals^2))
  } else {
    sd0 <- sd(y) * sqrt((n - 1) / n)
  }

  # if (m == 1) {
  #   if (!is.null(z)) {
  #     mu     <- unname(ls.out$coeff[1])
  #     gam <- unname(ls.out$coeff[2:(1 + p)])
  #     res    <- ls.out$residuals
  #     sigma <- sqrt(mean(res*res))
  #   } else {
  #     mu     <- mean(y)
  #     sigma  <- sd0
  #   }
  #
  #   loglik   <- - (n/2) *(1 + log(2*pi) + 2*log(sigma))
  #   aic      <- -2*loglik + 2*(m-1+2*m+p)
  #   bic      <- -2*loglik + log(n)*(m-1+2*m+p)
  #   penloglik <- loglik
  #
  #   parlist <- list(alpha = 1, mu = mu, sigma = sigma, gam = gam)
  #   coefficients <- c(alpha = 1, mu = mu, sigma = sigma, gam = gam)
  #   postprobs <- rep(1, n)
  #
  # } else {  # m >= 2

  # generate initial values
  tmp <- normalpanelmixPMLEinit(y = y, z = z, ninits = ninits.short, m = m)

  # the following values for (h, k, tau, an) are given by default
  # h       <- 0  # setting h=0 gives PMLE
  # k       <- 0  # k is set to 0 because this is PMLE
  # tau     <- 0.5  # tau is set to 0.5 because this is PMLE
  an <- 1 / n # penalty term for variance
  sigma0 <- rep(sd0, m)
  mu0 <- double(m + 1) # dummy
  
  if (is.null(z)) {
    ztilde <- matrix(0) # dummy
  } else {
    ztilde <- z
  }
  # short EM
  b0 <- as.matrix(rbind(tmp$alpha, tmp$mu, tmp$sigma, tmp$gam))
  if (!is.null(binit)) {
    b0[, 1] <- binit
  }

  out.short <- cppnormalpanelmixPMLE(
    b0, y, ztilde, mu0, sigma0, m, p, t, an, maxit.short,
    ninits.short, epsilon.short
  )

  # long EM
  components <- order(out.short$penloglikset, decreasing = TRUE)[1:ninits]
  b1 <- b0[, components] # b0 has been updated
  if ((!is.null(in.coefficient)) & (dim(b1)[1] == length(in.coefficient))) {
    b1[, components] <- in.coefficient
  }
  out <- cppnormalpanelmixPMLE(b1, y, ztilde, mu0, sigma0, m, p, t, an, maxit, ninits, epsilon)

  index <- which.max(out$penloglikset)
  alpha <- b1[1:m, index] # b0 has been updated
  mu <- b1[(1 + m):(2 * m), index]
  sigma <- b1[(1 + 2 * m):(3 * m), index]
  if (!is.null(z)) {
    gam <- b1[(3 * m + 1):(3 * m + p), index]
  }
  penloglik <- out$penloglikset[index]
  loglik <- out$loglikset[index]
  postprobs <- matrix(out$post[, index], nrow = n)

  aic <- -2 * loglik + 2 * (m - 1 + 2 * m + p)
  bic <- -2 * loglik + log(n) * (m - 1 + 2 * m + p)

  mu.order <- order(mu)
  alpha <- alpha[mu.order]
  mu <- mu[mu.order]
  sigma <- sigma[mu.order]

  postprobs <- postprobs[, mu.order]
  if (m > 1) {
    colnames(postprobs) <- c(paste("comp", ".", 1:m, sep = ""))
  }
  parlist <- list(alpha = alpha, mubeta = mu, sigma = sigma, gam = gam)
  coefficients <- unlist(parlist)

  # } # end m >= 2

  if (vcov.method == "none") {
    vcov <- NULL
  } else {
    vcov <- normalpanelmixVcov(y = y, coefficients = coefficients, z = z, vcov.method = vcov.method)
  }

  a <- list(
    coefficients = coefficients, parlist = parlist, vcov = vcov, loglik = loglik,
    penloglik = penloglik, aic = aic, bic = bic, postprobs = postprobs,
    components = getComponentcomponents(postprobs),
    call = match.call(), m = m, label = "PMLE"
  )

  class(a) <- "normalregMix"

  a
} # end function normalpanelmixPMLE



#' @description Generate initial values used by the PMLE of mixture of normals
#' @title normalpanelmixPMLEinit
#' @name normalpanelmixPMLEinit
#' @param y n by 1 vector of data
#' @param z n by p matrix of regressor associated with gamma
#' @param ninits number of initial values to be generated
#' @param m The number of components in the mixture
#' @return A list with the following items:
#' \item{alpha}{m by ninits matrix for alpha}
#' \item{mu}{m by ninits matrix for mu}
#' \item{sigma}{m by ninits matrix for sigma}
#' \item{gam}{m by ninits matrix for gam}
normalpanelmixPMLEinit <- function(y, z = NULL, ninits = 1, m = 2) {
  # if (normalregMix.test.on) # initial values controlled by normalregMix.test.on
  #   set.seed(normalregMix.test.seed)

  t <- nrow(y)
  n <- ncol(y)
  y <- as.vector(y)

  p <- ncol(z)
  gam <- NULL
  if (!is.null(z)) {
    out <- lsfit(z, y)
    gam0 <- out$coef[-1]
    gam <- matrix(runif(p * ninits, min = 0.5, max = 1.5), nrow = p) * gam0
    y <- out$residuals + out$coef[1]
  }

  alpha <- matrix(runif(m * ninits), nrow = m)
  alpha <- t(t(alpha) / colSums(alpha))
  mu <- matrix(runif(m * ninits, min = min(y), max = max(y)), nrow = m)
  sigma <- matrix(runif(m * ninits, min = 0.01, max = 2) * sd(y), nrow = m)

  list(alpha = alpha, mu = mu, sigma = sigma, gam = gam)
} # end function normalpanelmixPMLEinit


#' @description Computes the variance-covariance matrix of the MLE of
#' m-component normal mixture.
#' @title normalpanelmixVcov
#' @name normalpanelmixVcov
#' @param y n by 1 vector of data
#' @param coefficients (alpha_1, ..., alpha_m, mu_1, ..., mu_m, sigma_1, ..., sigma_m, gam)
#' @param z n by p matrix of regressor associated with gamma
#' @param vcov.method Method used to compute the variance-covariance matrix,
#' one of \code{"Hessian"} and \code{"OPG"}. #' The default option is \code{"Hessian"}.
#' When \code{method = "Hessian"}, the variance-covarince matrix is
#' estimated by the Hessian using the formula given in Boldea and Magnus (2009).
#' When \code{method = "OPG"}, the outer product of gradients is used.
#' @return The variance-covariance matrix of the MLE of
#' m-component normal mixture given the data and coefficients.
#' @references   Boldea, O. and Magnus, J. R. (2009)
#' Maximum Likelihood Estimation of the Multivariate Normal Mixture Model,
#' \emph{Journal of the American Statistical Association},
#' \bold{104}, 1539--1549.
normalpanelmixVcov <- function(y, coefficients, z = NULL, vcov.method = c("Hessian", "OPG"))
{
t  <- nrow(y)
n  <- ncol(y)
y <- as.vector(y)
len <- length(coefficients)
p <- 0
gam  <- NULL
vcov.method <- match.arg(vcov.method)

if (!is.null(z)) {
  z <- as.matrix(z)
  p <- ncol(z)
  gam <- coefficients[(len-p+1):len]
}

m <- (len-p)/3
if (round(m) != m) {
  stop("The dimension of the coefficients is incompatible with z. Please check the data.")
}

alpha   <- coefficients[1:m]
mu      <- coefficients[(m+1):(2*m)]
sigma   <- coefficients[(2*m+1):(3*m)]

if (m == 1) {
  if (is.null(z)) {
    I <- n*diag(c(1/sigma^2, 2/sigma^2))  # information matrix
  } else {
    z1 <- cbind(1, z)
    I <- matrix(0, nrow=p+2, ncol=p+2)
    I[1:(p+1), 1:(p+1)] <- t(z1) %*% z1/sigma^2
    I[(p+2), (p+2)] <- n*2/sigma^2
    s.1 <- c(1,(p+2), 2:(p+1))
    I <- I[s.1, s.1]
  }
  vcov  <- solve(I)
  # Because the variance is parameterized as sigma^2, we convert it to sigma
  c.mat.vec <- c(1, (1/sigma^(1/2))/2, rep(1, p))
  vcov <- diag(c.mat.vec) %*% vcov %*% diag(c.mat.vec)

} else { # end of if (m == 1)
  # m >= 2

  # Compute posterior probabilities, and adjust y if z is present
  sigma0  <- rep(1, m)  # dummy
  mu0     <- double(m)  # dummy
  an      <- 1/n  # penalty term for variance
  h       <- 0
  tau     <- 0.5
  k       <- 0
  epsilon <- 1e-08
  maxit = 2
  ninits = 1

  b <- matrix( rep( coefficients, ninits), ncol = ninits)
  if (is.null(z)) {
    out.p <- cppnormalpanelmixPMLE(b, y, matrix(0),  mu0, sigma0, m, p ,t , an, maxit, ninits, epsilon, tau, h, k)

  } else {
    out.p <- cppnormalpanelmixPMLE(b, y, z,  mu0, sigma0, m, p, t, an, maxit, ninits, epsilon, tau, h, k)
    # Adjust y
    y <- y - z %*% gam
  }
  post <- matrix(out.p$post, nrow=n)

  p1 <- seq(1, (2*m-1), by=2) # sequence of odd numbers, 1,3,...,2*m-1
  p2 <- seq(2, (2*m), by=2)    # sequence of even numbers, 2,4,...,2*m

  # Matrices used in computing vcov

  a <- diag(1/alpha[-m], nrow=m-1, ncol=m-1)
  a <- cbind(a, -1/alpha[m])  # m-1 by m matrix of a_i's
  abar <- a %*% t(post) # m-1 by n

  Z0 <- t((t(matrix(rep.int(y, m), ncol=m))-mu)/sigma)  # normalized data, n by m
  f <- t(t(exp(-Z0^2/2)/sqrt(2*pi))/sigma)      # pdf, n by m
  phi <- t(t(f)*alpha)                # n by m
  f0 <- rowSums(phi)                  # data pdf, n by 1

  vinv <- 1/(sigma*sigma)

  b <- t(t(Z0)/sigma)  # n by m
  B <- t(vinv - t(b*b))  # n by m

  c0 <- array(0, dim=c(n, m, 2))
  c0[, , 1] <- b
  c0[, , 2] <- -B/2

  # Computes Hessian-based I
  if (vcov.method == "Hessian")  {
    other.method = "OPG"
    C0 <- array(0, dim=c(n, m, 2, 2))
    C0[, , 1, 1] <- t(matrix(vinv, nrow=m, ncol=n))  # n by m
    C0[, , 2, 1] <- C0[, , 1, 2] <- t(t(b)*vinv)    # n by m
    C0[, , 2, 2] <- t((vinv -2*t(B))*vinv)/2       # n by m

    Q.pi <- - abar %*% t(abar)  # m-1 by m-1
    Q.pi.theta <- matrix(0, nrow=m-1, ncol=2*m)  # m-1 by 2m
    for (i in 1:m){
      zi <- a[, i] - abar  # m-1 by n
      wi <- c0[, i, ]*post[, i]  # n by 2
      Q.i <- colSums(tKR(wi, t(zi)))  # 2*(m-1) vector
      # first m-1 elements correspond to mu x pi
      # second m-1 elements correspond to sigma x pi
      Q.pi.theta[, i] <- Q.i[1:(m-1)]
      Q.pi.theta[, m+i] <- Q.i[m:(2*(m-1))]
    }

    Q.theta <- matrix(0, nrow=2*m, ncol=2*m)
    for (i in 2:m){ # off-diagonal blocks
      for (j in 1:(i-1)){
        wi  <- c0[, i, ]*post[, i]
        wj  <- c0[, j, ]*post[, j]
        Q.ij <- - colSums(tKR(wi, wj))
        Q.theta[(2*i-1):(2*i), (2*j-1):(2*j)] = t(matrix(Q.ij, nrow=2, ncol=2))
      }
    }

    Q.theta <- Q.theta + t(Q.theta)
    for (i in 1:m){ # diagonal blocks
      C.ii <- array(C0[, i, , ], dim=c(n, 2, 2))
      Q.ii.1 <- apply(C.ii*post[, i], c(2, 3), sum)
      w.ii <- tKR(c0[, i, ], c0[, i, ])*post[, i]*(1-post[, i])
      Q.ii.2 <- matrix(colSums(w.ii), nrow=2, ncol=2)
      Q.theta[(2*i-1):(2*i), (2*i-1):(2*i)] <- -Q.ii.1 + Q.ii.2
    }
    # odd rows and columns of Q.theta correspond to mu
    # even rows and columns of Q.theta correspond to sigma
    Q.theta <- Q.theta[c(p1,p2), c(p1,p2)] # first block = wrt mu, second blosk = wrt sigma

    dimI <- m-1+2*m
    I <- matrix(0, nrow=dimI, ncol=dimI)
    I[1:(m-1), 1:(m-1)] <- - Q.pi
    I[1:(m-1), m:dimI]  <- - Q.pi.theta
    I[m:dimI, 1:(m-1)]  <- - t(Q.pi.theta)
    I[m:dimI, m:dimI]   <- - Q.theta

    if (!is.null(z)) {
      dbar  <-  z*rowSums(post*b) # n by p
      Q.gam.theta <- matrix(0, nrow=p, ncol=2*m)  # p by 2*m matrix
      for (i in 1:m) {
        C.i <- array(C0[, i, 1, ], dim=c(n, 2))  # n by 2
        Q.i.1 <- colSums(tKR(-C.i+b[, i]*c0[, i, ], z*post[, i]))
        Q.i.2 <- colSums(tKR(c0[, i, ]*post[, i], dbar))  # p*2 vector
        Q.gam.theta[, (2*i-1):(2*i)] <- matrix(Q.i.1-Q.i.2, nrow=p, ncol=2)
      }

      Q.gam.theta <- Q.gam.theta[, c(p1, p2),  drop=FALSE]  # p by 2*m
      w1 <- (post*b)%*%t(a) - rowSums(post*b)*t(abar)  # n by m-1
      Q.pi.gam.0 <- colSums(tKR(w1, z))  # (m-1)*p vector
      Q.pi.gam  <- matrix(Q.pi.gam.0, nrow=m-1, ncol=p)
      Q.gam     <- - t(z) %*% (z*rowSums(post*B)) -
        matrix(colSums(tKR(dbar, dbar)), nrow=p, ncol=p)
      I <- cbind(I, -rbind(Q.pi.gam, t(Q.gam.theta)))
      I <- rbind(I, -cbind(t(Q.pi.gam), Q.gam.theta, Q.gam))
    }  # end if (!is.null(z))

  } else  { # compute I with (vcov.method == "OPG")
    other.method = "Hessian"
    score <- t(abar)
    for (j in 1:m) { score <- cbind(score, c0[, j, ]*post[, j]) }

    ind <- c(c(1:(m-1)), p1+m-1, p2+m-1)
    score <- score[, ind]
    I <- t(score) %*% score

    if (!is.null(z)) {
      dbar  <-  z*rowSums(post*b) # n by p
      score <- cbind(score, dbar)
      I <- t(score) %*% score
    }

  } # end if (vcov.method == "OPG")

  vcov <- try(solve(I))
  if (class(vcov) == "try-error" || any(diag(vcov) <0) ) {
    vcov <- matrix(NaN, nrow = 3*m-1+p, ncol = 3*m-1+p)
    warning("Fisher information matrix is singular and/or the
            variance is estimated to be negative. Consider using vcov.method=\"",other.method,"\".")
  }

  # Because the variance is parameterized as sigma^2, we convert is to sigma

  c.mat.vec <- c(rep(1, m-1+m), (1/sigma^(1/2))/2, rep(1, p))
  vcov <- diag(c.mat.vec) %*% vcov %*% diag(c.mat.vec)

  # Add the variance of alpha_m
  len   <- length(coefficients)
  M.mat <- diag(len-1)
  M.mat <- rbind(M.mat[1:(m-1), ], c(rep(-1, m-1), rep(0, len-m)),
                 M.mat[m:(len-1), ])

  vcov     <- M.mat %*% vcov %*% t(M.mat)

}   # end else (i.e., m >= 2)

vcov

}  # end function normalpanelmixVcov
