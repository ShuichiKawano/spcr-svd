### spcrSVD

### INPUT
# x: design matrix
# y: a vector of response variable
# k: the number of principal components
# lambda_V: lambda_V
# lambda_beta: lambda_beta
# rho1, rho2, rho3: penalty parameter in ADMM - default is 1.
# tol: tolerance of the algorithm - default is 1e-5.
# beta_int_i, beta_0_i, beta_i, V0_i, V1_i, V_i, Lambda_1_i, Lambda_2_i, lambda_3_i
#   : initial values of the parameters - default is NULL. This is optional. 
#     Initial values are automatically generated in this algorithm. 

### OUTPUT
# V: estimate of the parameter V
# V0: estimate of the parameter V0
# V1: estimate of the parameter V1
# beta: estimate of the parameter beta
# beta0: estimate of the parameter beta0
# beta_int: estimate of the intercept
# Lambda_1: estimate of the parameter Lambda_1
# Lambda_2: estimate of the parameter Lambda_2
# lambda_3: estimate of the parameter lambda_3

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

### ADMM for spcrSVD
spcrSVD <- function(x, y, k, w, lambda_V, lambda_beta, rho1=1, rho2=1, rho3=1, tol=1e-5, beta_int_i=NULL, beta_0_i=NULL, 
                     beta_i=NULL, V0_i=NULL, V1_i=NULL, V_i=NULL, Lambda_1_i=NULL, Lambda_2_i=NULL, lambda_3_i=NULL)
{
  n <- nrow(x)
  p <- ncol(x)
  
  ### START initialization
  if( is.null(beta_int_i) ){
    ### Not prepared initial values
    beta_int <- mean(y)
    beta_0 <- beta <- rep(0, k)
    V0 <- V1 <- V <- as.matrix( svd(x)$v[ ,1:k] )
    Lambda_2 <- Lambda_1 <- matrix( 0, nrow(V), ncol(V) )
    lambda_3 <- rep(0, k)
  }else{
    ### Prepared initial values
    beta_int <- beta_int_i
    beta_0 <- beta_0_i
    beta <- beta_i
    V0 <- V0_i
    V1 <- V1_i
    V <- V_i
    Lambda_2 <- Lambda_2_i
    Lambda_1 <- Lambda_1_i
    lambda_3 <- lambda_3_i
  }
  ### END initialization
  
  loss_old <- 2
  loss_new <- 1

  itr <- 1
  while( ( abs(loss_old - loss_new) / abs(loss_old) ) >= tol )
  {	
    ### Estimation of Z
    Z <- x %*% V
    
    ### Estimation of V
    SVD_V <- svd( (w/n) * t(x) %*% Z + 0.5 * rho1 * ( V0 - Lambda_1 ) )
    V <- as.matrix(SVD_V$u %*% t(SVD_V$v))
    
    ### Estimation of V_0
    for( i in 1:p ) for( j in 1:k ) V0[ i, j ] <- softsh( (rho1*(V[ i, j ] + Lambda_2[ i, j ]) + rho2*(V[ i, j ] + Lambda_2[ i, j ]))/(rho1+rho2), 
                                                          lambda_V/(rho1 + rho2))
    V0 <- sweep(V0, 2, sqrt(apply((V0+1e-10)^2, 2 ,sum)), FUN ="/")
    
    ### Estimation of V_1
    vecV1 <- solve( kronecker(beta%*%t(beta), t(x)%*%x)/n + 0.5*rho2*kronecker(diag(k), diag(p)) )%*%c( t(x)%*%(y - beta_int)%*%t(beta)/n + 0.5*rho2*(V0 - Lambda_2) )
    V1 <- matrix(vecV1, p, k)
    V1 <- sweep(V1, 2, sqrt(apply((V1+1e-10)^2, 2 ,sum)), FUN ="/")
    
    ### Estimation of beta
    beta <- solve( t(V1)%*%t(x)%*%x%*%V1/n + 0.5*rho3*diag(k) ) %*% ( t(V1) %*% t(x) %*% (y - beta_int)/n + 0.5*rho3*(beta_0 - lambda_3) )
    
    ### Estimation of beta_0
    for( i in 1:k ) beta_0[ i ] <- softsh( beta[ i ] + lambda_3[ i ], lambda_beta/rho3 )
    
    ### Estimation of beta_int
    beta_int <- mean( y - x %*% V1 %*% beta )
    
    ### Estimation of Lambda_1
    Lambda_1 <- Lambda_1 + ( V - V0 )
    
    ### Estimation of Lambda_2
    Lambda_2 <- Lambda_2 + ( V1 - V0 )
    
    ### Estimation of lambda_3
    lambda_3 <- lambda_3 + ( beta - beta_0 )
    
    loss_old <- loss_new
    loss_new <- mean( (y - beta_int - x%*%V1%*%beta)^2 ) + w * sum(diag( t(x - Z%*%t(V)) %*% (x - Z%*%t(V))))/(n*p) + 
    lambda_V*sum(abs(V0))/(p*k) + lambda_beta*sum(abs(beta_0))/k
    itr <- itr + 1
    if(itr == 10000) break
  }
  return(list(V=V, V0=V0, V1=V1, beta=beta, beta0=beta_0, beta_int=beta_int, Lambda_1=Lambda_1, Lambda_2=Lambda_2, lambda_3=lambda_3))
}
