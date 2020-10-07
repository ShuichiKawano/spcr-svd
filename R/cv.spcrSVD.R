### Cross-validation for spcrSVD

### INPUT
# x: design matrix
# y: a vector of response variable
# k: the number of principal components
# w: tuning parameter on the PCA loss function
# fold: the number of fold - default is 5.
# rho1, rho2, rho3: penalty parameter in ADMM - default is 1.
# tol: tolerance of the algorithm - default is 1e-5.

### OUTPUT
# lambda_V_CV: lambda_V that minimizes CV
# lambda_beta_CV: lambda_beta that minimizes CV
# CV_mat: matrix of CV 
# lambda_V_seq: sequence of lambda_V
# lambda_beta_seq: sequence of lambda_beta

source("R/spcrSVD.R")

cv.spcrSVD <- function(x, y, k, w, fold, rho1=1, rho2=1, rho3=1, tol=1e-5){
  
  n <- nrow(x)
  p <- ncol(x)
  
  ### determination of lambda_V_max
  V0 <- V1 <- V <- as.matrix( svd(x)$v[ ,1:k] )
  Z <- x %*% V
  SVD_V <- svd( (w/n) * t(x) %*% Z + 0.5 * rho1 * V0 )
  V <- as.matrix(SVD_V$u %*% t(SVD_V$v))
  lambda_V_max <- (rho1+rho2)*max(abs(V)) - 1e-5
  
  ### determination of lambda_beta_max
  beta <- solve( t(V1)%*%t(x)%*%x%*%V1/n + 0.5*rho3*diag(k) ) %*% ( t(V1) %*% t(x) %*% (y - mean(y))/n )
  lambda_beta_max <- rho3*max(abs(beta)) - 1e-5

  ### Candidate values of tuning parameters (lamdba_beta, lambda_gamma)
  lambda_V_candidate <- rev( seq( 0.005, lambda_V_max, length=10 ) )
  lambda_beta_candidate <- rev( seq( 0.005, lambda_beta_max, length=10 ) )

  ### Preparation of CV
  temp_num <- 1:n
  extract_num <- as.list(NULL)
  for(i in 1:fold)
  {
    temp_num2 <- sample(temp_num, n/fold, replace=FALSE)
    extract_num[i][[1]] <- temp_num2
    temp_num <- setdiff( temp_num, temp_num2 )
  }
  
  ### CV_mat : estimated CV errors
  CV_mat <- matrix( 0, length(lambda_V_candidate), length(lambda_beta_candidate) )
  
  xall <- x
  yall <- y
  
  ### "fold"-fold CV
  for(itr_CV in 1:fold)
  {	
    x_ori <- xall[ -extract_num[[itr_CV]], ]
    y_ori <- yall[ -extract_num[[itr_CV]] ]
    x_test_cv <- xall[ extract_num[[itr_CV]], ]
    y_test_cv <- yall[ extract_num[[itr_CV]] ]
    
    y <- y_ori - mean(y_ori)
    x <- sweep(x_ori, 2, apply(x_ori,2,mean))
    y_test_cv <- y_test_cv - mean(y_ori)
    x_test_cv = sweep(x_test_cv, 2, apply(x_ori,2,mean))
    
    ### Initialization when itr_CV=1
    if( itr_CV == 1 ){
      beta_int <- mean(y)
      beta_0 <- beta <- rep(0, k)
      V0 <- V1 <- V <- as.matrix( svd(x)$v[ ,1:k] )
      Lambda_2 <- Lambda_1 <- matrix( 0, nrow(V), ncol(V) )
      lambda_3 <- rep(0, k)
    }

    ####### START Estimate parameters
    for( itr_lambda_V in 1:length(lambda_V_candidate) )
    {
      lambda_V <- lambda_V_candidate[itr_lambda_V]
      for( itr_lambda_beta in 1:length(lambda_beta_candidate) )
      {
        lambda_beta <- lambda_beta_candidate[itr_lambda_beta]
        result <- spcrSVD(x=x, y=y, k=k, w=w, lambda_V=lambda_V, lambda_beta=lambda_beta, rho1=rho1, rho2=rho2, rho3=rho3, tol=tol,
                           beta_int_i=beta_int, beta_0_i=beta_0, beta_i=beta, V0_i=V0, V1_i=V1, V_i=V, Lambda_1_i=Lambda_1, Lambda_2_i=Lambda_2, 
                           lambda_3_i=lambda_3)
        ### CV-error
        s_cv <- mean( ( y_test_cv - result$beta_int - x_test_cv %*% result$V1 %*% result$beta )^2 )
        ### Strock of CV-error
        CV_mat[ itr_lambda_V, itr_lambda_beta ] <- CV_mat[ itr_lambda_V, itr_lambda_beta ] + s_cv
        
        ### Initial values of next steps
        beta_int <- result$beta_int
        beta_0 <- result$beta0
        beta <- result$beta
        V0 <- result$V0
        V1 <- result$V1
        V <- result$V
        Lambda_1 <- result$Lambda_1
        Lambda_2 <- result$Lambda_2
        lambda_3 <-result$lambda_3
        
        ### stock of the results for the first itr_lambda_beta
        if( itr_lambda_beta == 1 ) INIT_itr <- result
        ### stock of the results for the first step among all setting
        if( (itr_CV == 1) && (itr_lambda_beta == 1) && (itr_lambda_V == 1) ) INIT_ITR <- result 
        
      }
      ### initial values for the next itr_lambda_V
      beta_int <- INIT_itr$beta_int
      beta_0 <- INIT_itr$beta0
      beta <- INIT_itr$beta
      V0 <- INIT_itr$V0
      V1 <- INIT_itr$V1
      V <- INIT_itr$V
      Lambda_1 <- INIT_itr$Lambda_1
      Lambda_2 <- INIT_itr$Lambda_2
      lambda_3 <-INIT_itr$lambda_3
      
    }
    ####### END Estimate parameters
    CV_mat <- CV_mat/fold
    
    ### initial values for the next CV_itr
    beta_int <- INIT_ITR$beta_int
    beta_0 <- INIT_ITR$beta0
    beta <- INIT_ITR$beta
    V0 <- INIT_ITR$V0
    V1 <- INIT_ITR$V1
    V <- INIT_ITR$V
    Lambda_1 <- INIT_ITR$Lambda_1
    Lambda_2 <- INIT_ITR$Lambda_2
    lambda_3 <-INIT_ITR$lambda_3
  }
  
  ### START Search of min CV
  whichiminCandi.col <- c( which.min(CV_mat[1, ]), which.min(CV_mat[2, ]), which.min(CV_mat[3, ]), which.min(CV_mat[4, ]), which.min(CV_mat[5, ]), 
                           which.min(CV_mat[6, ]), which.min(CV_mat[7, ]), which.min(CV_mat[8, ]), which.min(CV_mat[9, ]), which.min(CV_mat[10, ]) )
  minCandi.col <- c( min(CV_mat[1, ]), min(CV_mat[2, ]), min(CV_mat[3, ]), min(CV_mat[4, ]), min(CV_mat[5, ]), min(CV_mat[6, ]), min(CV_mat[7, ]), 
                     min(CV_mat[8, ]), min(CV_mat[9, ]), min(CV_mat[10, ]) )
  
  whichiminCandi.row <- c( which.min(CV_mat[ ,1]), which.min(CV_mat[ ,2]), which.min(CV_mat[ ,3]), which.min(CV_mat[ ,4]), which.min(CV_mat[ ,5]), 
                           which.min(CV_mat[ ,6]), which.min(CV_mat[ ,7]), which.min(CV_mat[ ,8]), which.min(CV_mat[ ,9]), which.min(CV_mat[ ,10]) )
  minCandi.row <- c( min(CV_mat[ ,1]), min(CV_mat[ ,2]), min(CV_mat[ ,3]), min(CV_mat[ ,4]), min(CV_mat[ ,5]), min(CV_mat[ ,6]), min(CV_mat[ ,7]), 
                     min(CV_mat[ ,8]), min(CV_mat[ ,9]), min(CV_mat[ ,10]) )
  ### END Search of min CV
  
  ### Selected tuning parameters by CV
  lambda_V <- lambda_V_candidate[ whichiminCandi.row[ which.min(minCandi.row) ] ]
  lambda_beta <- lambda_beta_candidate[ whichiminCandi.col[ which.min(minCandi.col) ] ]
  return(list(lambda_V_CV=lambda_V, lambda_beta_CV=lambda_beta, CV_mat=CV_mat, lambda_V_seq=lambda_V_candidate, lambda_beta_seq=lambda_beta_candidate))
}
