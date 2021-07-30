# Sparse principal component regression via singular value decomposition approach. 
This is an R source code for performing sparse principal component regression via singular value decomposition approach (spcrSVD). The directory *R* consists of four files as follows. 
- **spcrSVD.R** provides spcrSVD via an ADMM algorithm. 
- **cv.spcrSVD.R** computes cross-validation of spcrSVD. 
- **spcrSVD_linearized.R** provides spcrSVD via a linearized ADMM algorithm. 
- **cv.spcrSVD_linearized.R** computes cross-validation of spcrSVD_linearized. 

spcrSVD is introduced in the paper:
Kawano, S. (2021) Sparse principal component regression via singular value decomposition approach. *Advances in Data Analysis and Classification (in press)* (doi:[ 10.1007/s11634-020-00435-2](https://doi.org/10.1007/s11634-020-00435-2)).

## Usage example
Read source files.
```
source("R/spcrSVD.R")
source("R/cv.spcrSVD.R")
source("R/spcrSVD_linearized.R")
source("R/cv.spcrSVD_linearized.R")
```

Setting of simulation.
```
library(MASS)
n <- 50; ErrorVariance <- 1; k <- 1; np <- 2; dummy_np <- 8
w <- 1e-1; rho1 <- 1; rho2 <- 1; rho3 <- 1; rho_V <- 1; rho_beta <- 1; tol <- 1e-5
Sigma <- diag( rep(1,np) )
for(i in 1:nrow(Sigma)) for(j in 1:ncol(Sigma)) Sigma[i ,j] = 0; Sigma[1,1] = 1; Sigma[2,2] = 1
nu0 <- c(2,1)
x_ori <- mvrnorm(n, rep(0, nrow(Sigma)), Sigma)
x_dummy <- matrix( rnorm(dummy_np*n), n, dummy_np )
x_ori <- cbind(x_ori, x_dummy)
y_ori <- nu0[1]*x_ori[ ,1] + nu0[2]*x_ori[ ,2] + rnorm(n, 0, ErrorVariance)
y <- y_ori - mean(y_ori)
x <- sweep(x_ori, 2, apply(x_ori,2,mean))
```

Perform spcrSVD
```
# Perform spcrSVD in the file spcrSVD.R
spcrSVD(x=x, y=y, k=k, w=w, lambda_V=1e-1, lambda_beta=1e-1, rho1=rho1, rho2=rho2, rho3=rho3, tol=tol)

# Perform cv.spcrSVD with five-fold in the file cv.spcrSVD.R
cv.spcrSVD(x=x, y=y, k=k, w=w, fold=5, rho1=rho1, rho2=rho2, rho3=rho3, tol=tol)

# Perform spcrSVD_linearized in the file spcrSVD_linearized.R
spcrSVD_linearized(x=x, y=y, k=k, w=w, lambda_V=1e-1, lambda_beta=1e-1, rho_V=rho_V, rho_beta=rho_beta, tol=tol)

# Perform cv.spcrSVD_linearized with five-fold in the file cv.spcrSVD_linearized.R
cv.spcrSVD_linearized(x=x, y=y, k=k, w=w, fold=5, rho_V=rho_V, rho_beta=rho_beta, tol=tol)
```
