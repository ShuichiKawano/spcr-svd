### spcrSVD-LADMM

### INPUT
# x: design matrix
# y: a vector of response variable
# k: the number of principal components
# lambda_V: lambda_V
# lambda_beta: lambda_beta
# rho_V, rho_beta: penalty parameter in LADMM - default is 1.
# tol: tolerance of the algorithm - default is 1e-5.
# beta_int_i, beta_0_i, beta_i, V0_i, V1_i, V_i, Lambda_1_i, Lambda_2_i, lambda_3_i
#   : initial values of the parameters - default is NULL. This is optional. 
#     Initial values are automatically generated in this algorithm. 

### OUTPUT
# V: estimate of the parameter V
# V0: estimate of the parameter V0
# beta: estimate of the parameter beta
# beta0: estimate of the parameter beta0
# beta_int: estimate of the intercept
# Lambda_0: estimate of the parameter Lambda_0
# lambda_0: estimate of the parameter lambda_0

### Soft-threshhold
softsh <- function(z, eta)
{
  if( eta < abs(z) ){
    if( z > 0 ) {
      z - eta
    } else {
      z + eta
    }
  } else {
    0
  }
}

### Linear ADMM for spcrSVD
spcrSVD_linearized <- function(x, y, k, w, lambda_V, lambda_beta, rho_V=1, rho_beta=1, tol=1e-5,
                     beta_int_i=NULL, beta_0_i=NULL, beta_i=NULL, V0_i=NULL,V_i=NULL, Lambda_0_i=NULL, lambda_0_i=NULL)
{
  n <- nrow(x)
  p <- ncol(x)
  
  ### START initialization
  if( is.null(beta_int_i) ){
    # Not prepared initial values
    beta_int <- mean(y)
    beta_0 <- beta <- rep(0, k)
    V0 <- V <- as.matrix( svd(x)$v[ ,1:k] )
    Lambda_0 <- matrix( 0, nrow(V), ncol(V) )
    lambda_0 <- rep(0, k)
  }else{
    # Prepared initial values
    beta_int <- beta_int_i
    beta_0 <- beta_0_i
    beta <- beta_i
    V0 <- V0_i
    V <- V_i
    Lambda_0 <- Lambda_0_i
    lambda_0 <- lambda_0_i
  }
  ### END initialization
  
  loss_old <- 2
  loss_new <- 1
  itr <- 1
  while( ( abs(loss_old - loss_new) / abs(loss_old) ) >= tol )
  {
    ### Spectral radius : nu (the first eigenvalue of the matrix beta t(beta) kronecker t(x) x, which is the second derivative of t(beta) t(V) t(x) x V beta with respect to V)
    NU <- kronecker(beta%*%t(beta), t(x)%*%x)
    nu <- max(eigen(NU)$values)

    ### Estimation of Z
    Z <- x %*% V0
    
    ### Estimation of V0
    SVD_V <- svd( w*t(x)%*%Z/n + 0.5*rho_V*(V+Lambda_0) )
    V0 <- as.matrix(SVD_V$u %*% t(SVD_V$v))

    ### Estimation of V
    Temp <- (2*n/(n*rho_V+2*nu))*( (t(x)%*%(y - beta_int)%*%t(beta) - t(x)%*%x%*%V%*%beta%*%t(beta))/n + nu*V/n - 0.5*rho_V*(Lambda_0 - V0) )
    for( i in 1:p ) for( j in 1:k ) V[ i, j ] <- softsh( Temp[i, j], lambda_V/( (n*rho_V+2*nu)/n ) )
    V <- sweep(V, 2, sqrt(apply((V+1e-10)^2, 2 ,sum)), FUN ="/")
    
    ### Estimation of beta
    beta <- solve( t(V)%*%t(x)%*%x%*%V/n + 0.5*rho_beta*diag(k) ) %*% ( t(V) %*% t(x) %*% (y - beta_int)/n + 0.5*rho_beta*(beta_0 - lambda_0) )
    
    ### Estimation of beta_0
    for( i in 1:k ) beta_0[ i ] <- softsh( beta[ i ] + lambda_0[ i ], lambda_beta/rho_beta )
    
    ### Estimation of beta_int
    beta_int <- mean( y - x %*% V %*% beta )
    
    ### Estimation of Lambda_0
    Lambda_0 <- Lambda_0 + V - V0
    
    ### Estimation of lambda_1
    lambda_0 <- lambda_0 + beta - beta_0
    
    loss_old <- loss_new
    loss_new <- mean( (y - beta_int - x%*%V%*%beta)^2 ) + w * sum(diag( t(x - Z%*%t(V0)) %*% (x - Z%*%t(V0))))/(n*p) + 
    lambda_V*sum(abs(V))/(p*k) + lambda_beta*sum(abs(beta_0))/k
    itr <- itr + 1
    if(itr == 10000) break
  }
  return(list(V=V, V0=V0, beta=beta, beta0=beta_0, beta_int=beta_int, Lambda_0=Lambda_0, lambda_0=lambda_0))
}
