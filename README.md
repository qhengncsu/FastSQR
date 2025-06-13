## Installation
To install the package, the right development tools must be in place. Specifically,
Windows users need [RTools](https://cran.r-project.org/bin/windows/Rtools/), while
Mac users need [Xcode and a Fortran compiler](https://cran.r-project.org/bin/macosx/tools/), 
and type the following commands in a terminal to download and extract an OpenMP runtime:
```
curl -O https://mac.r-project.org/openmp/openmp-17.0.6-darwin20-Release.tar.gz
sudo tar fvxz openmp-17.0.6-darwin20-Release.tar.gz -C /
```
Linux users need a GCC compiler. Most Linux systems in a server environment already have 
a GCC compiler with OpenMP support available. If that is not the case, you can try
```
sudo apt install build-essential
```
To install the pacakge in R, use the following R code:
```R
library(devtools)
install_github("qhengncsu/FastSQR")
```

## Quick Start
To start using the package, we provide the following simulated example:
```R
# generate data
set.seed(1234)
# design matrix is 1000*10000
X = matrix(rnorm(1000*10000),nrow=1000)
# coefficient is 10 1s followed by 9900 zeros
beta = c(rep(1,10),rep(0,ncol(X)-10))
# noise is t distribution with 2 degrees of freedom
y = X%*%beta+rt(1000, df=2)+5 #intercept is 5

# fit a quantile regression model at the 0.5 quantile,
# using smoothing parameter h=0.25, and along 50 values
# of penalty parameters lambda, the algorithm chooses
# its own sequence of lambda
res1 = myglmnet(X,y,family="quantile", nlambda=50, h=0.25, tau = 0.5)
# type res1 to view the degree of freedoms, the proportion
# of deviance explained, and the sequence of lambda

# fit a quantile regression model at the 0.75 quantile
res2 = myglmnet(X,y,family="quantile", nlambda=50, h=0.25, tau = 0.5)

# fit a quantile regression model using user specified sequence of lambda
lambdas = rev(seq(0.02,0.1,length=50))
res2 = myglmnet(X,y,family="quantile", lambda=lambdas, nlambda=50, h=0.25, tau = 0.5)
```
