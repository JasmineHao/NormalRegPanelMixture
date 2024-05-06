#' Add together two numbers.
#'
#' @export
#' @param alpha m vector.
#' @param mu m vector.
#' @param sigma m vector.
#' @param gamma p vector.
#' @param beta q vector
#' @param N int
#' @param T int
#' @param M number of types
#' @param p columns of z
#' @param q columns of x
#' @return list of items
#' \item{Y}
#' \item{x}
#' \item{z}
#' @examples
#'generateData(alpha,mu,sigma,gamma,beta,N,T,M,p,q)
# The example of x, z

# x = matrix(rnorm(N*T*q),nc=q)
# z = matrix(rnorm(N*T*p),nc=p)
generateData <- function(alpha,mu,sigma,gamma,beta,N,T,M,p,q,x=NULL,z=NULL){
  sprintf("N = %d",N)
  sprintf("T = %d",T)

  R <- matrix(0,nrow=N,ncol=M)
  if (sum(alpha) != 1){
    alpha <- alpha / sum(alpha)
  }
  
  if (length(alpha) != M | length(mu) != M){
    stop('M must be the size of alpha')
  }
  
  prior <- runif(N)
  alpha.cum <- c(0,cumsum(alpha))
   
  if (M > 1){
    for (m in 1:M){
      lb <-  alpha.cum[m]
      ub <- alpha.cum[m+1] 
      
      R[,m] <- 1 * ((prior <= ub) & (prior > lb))
      }
  }else{
    R <- matrix(1,nrow=N,ncol=M)
  }
  Y = matrix(0,T,N)

  if ((q != 0) && (is.null(x))){
    x = matrix(rnorm(N*T*q),nc=q)
  }
  if ((p != 0) && (is.null(z))){
    x = matrix(rnorm(N*T*p),nc=p)
  }
  

  mu_R <- R %*% mu
  sigma_R <- R %*% sigma
  u <- matrix(rnorm(N*T),nrow=T,ncol=N)

  for (nn in 1:N) {
    y_nn <- rep(0,T)
    y_nn <- mu_R[nn] + sigma_R[nn] * u[,nn]
    if (q > 1 ){
      beta_R <- R %*% beta
      y_nn = y_nn + x[(T*(nn-1)+1) : (T*nn),] %*% beta_R[nn,]
     }else if (q == 1){
       beta_R <- R %*% as.vector(beta)
       y_nn = y_nn + x[(T*(nn-1)+1) : (T*nn),] * beta_R[nn]
     }else{}

    if (p > 1){
      y_nn = y_nn + z[(T*(nn-1)+1) : (T*nn),] %*% gamma
    }else if (p ==1){
      y_nn = y_nn + z[(T*(nn-1)+1) : (T*nn),] * gamma
    }else{ }
    Y[,nn] <- y_nn
  }
  if (p == 0){z = NULL}
  if (q == 0){x = NULL}
  
  return(list(Y=Y,Z=z,X=x))
  }



#' Add together two numbers.
#'
#' @export
#' @param alpha m vector.
#' @param mu m vector.
#' @param sigma m vector.
#' @param gamma p vector.
#' @param mu0 m vector.
#' @param sigma0 m vector.
#' @param gamma p vector.
#' @param beta q vector
#' @param N int
#' @param T int
#' @param M number of types
#' @param p columns of z
#' @param q columns of x
#' @return list of items
#' \item{Y}
#' \item{x}
#' \item{z}
#' @examples
#'generateDataAR1(alpha,mu,sigma,gamma,beta,N,T,M,p,q)
# The example of x, z

# x = matrix(rnorm(N*T*q),nc=q)
# z = matrix(rnorm(N*T*p),nc=p)
generateDataAR1 <- function(alpha,mu,sigma,gamma,beta,mu0,sigma0,gamma0,beta0,N,T,M,p,q,p.0,q.0,x=NULL,z=NULL, x0 = NULL, z0=NULL){
  sprintf("N = %d",N)
  sprintf("T = %d",T)

  # Draw the initial type
  R <- matrix(0,nrow=N,ncol=M)
  if (sum(alpha) != 1){
    alpha <- alpha / sum(alpha)
  }
  
  if (length(alpha) != M | length(mu) != M){
    stop('M must be the size of alpha')
  }
  
  prior <- runif(N)
  alpha.cum <- c(0,cumsum(alpha))
   
  if (M > 1){
    for (m in 1:M){
      lb <-  alpha.cum[m]
      ub <- alpha.cum[m+1] 
      
      R[,m] <- 1 * ((prior <= ub) & (prior > lb))
      }
  }else{
    R <- matrix(1,nrow=N,ncol=M)
  }
  Y = matrix(0,T,N)

  if ((q != 0) && (is.null(x))){
    x = matrix(rnorm(N*T*q),nc=q)
  }
  if ((p != 0) && (is.null(z))){
    z = matrix(rnorm(N*T*p),nc=p)
  }
  
  if ((q.0 != 0) && (is.null(x0))) {
    x0 = matrix(rnorm(N * q.0), nc = q.0)
  }
  if ((p != 0) && (is.null(z))) {
    z0 = matrix(rnorm(N * p), nc = p)
  }

  mu_R <- R %*% mu
  sigma_R <- R %*% sigma

  mu_R0 <- R %*% mu0
  sigma_R0 <- R %*% sigma0
  u <- matrix(rnorm(N*T),nrow=T,ncol=N)
  
  
  y0 <- mu_R0 + sigma_R0 * rnorm(N)
  
  if (q.0 > 1) {
    beta_R0 <- R %*% beta0
    y0 = y0 + rowSums(x0 * beta_R0)
  } else if (q.0 == 1) {
    beta_R0 <- R %*% as.vector(beta0)
    
    y0 = y0 + x0 * beta_R0
  } else {

  }
  
  if (q > 1) {
    beta_R <- R %*% beta
  } else if (q == 1) {
    beta_R <- R %*% as.vector(beta)
  }
     
  for (nn in 1:N) {
    y_nn <- rep(0, T)
    y_nn <- mu_R[nn] + sigma_R[nn] * u[, nn]
    y_l1 <- y0[nn]
    for (tt in 1:T){
      x[(T * (nn - 1) + tt), 1] <- y_l1 # update the lagged value
      x_t <- x[(T * (nn - 1) + tt), ]
      y_nn[tt] <- y_nn[tt] + x_t %*% beta_R[nn, ]
      
      if (p > 1){
        y_nn = y_nn + z[(T*(nn-1)+1) : (T*nn),] %*% gamma
      }else if (p ==1){
        y_nn = y_nn + z[(T*(nn-1)+1) : (T*nn),] * gamma
      }else{ }
      y_l1 <- y_nn[tt]
    }
    Y[,nn] <- y_nn
  }
  
  
  return(list(Y=Y,Z=z,X=x, Y0=y0, X0=x0, Z0=z0))
  }
