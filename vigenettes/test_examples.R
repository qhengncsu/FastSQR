set.seed(1234)
X = matrix(rnorm(10000*10000),nrow=10000)
beta = c(rep(1,10),rep(0,ncol(X)-10))
y = X%*%beta+rt(10000, df=2)+5
#y = as.double(rbinom(1000,size=1,prob=0.5))
library(myglmnet)
library(conquer)
start.time <- Sys.time()
res1 = myglmnet(X,y,family="quantile",nlambda=50, h = 0.25, tau=0.9)
end.time <- Sys.time()
print(end.time-start.time)
lambdas = res1$lambda
start.time <- Sys.time()
res2 = conquer.reg(X,y,lambda=lambdas,h=0.25,para.elastic=1,epsilon=1e-4,tau=0.9)
end.time <- Sys.time()
print(end.time-start.time)
tau = 0.5
quantile_loss1 <- double(50)
quantile_loss2 <- double(50)
for(j in 1:50){
  pred1 <-  X%*%res1$beta[,j] + res1$a0[j]
  quantile_loss1[j] <- mean(0.5*abs(y-pred1)+(tau-0.5)*(y-pred1))
  pred2 <-  X%*%res2$coeff[2:50001,j] + res2$coeff[1,j]
  quantile_loss2[j] <- mean(0.5*abs(y-pred2)+(tau-0.5)*(y-pred2))
}
quantile_loss1
rev(quantile_loss2)

best_fit <- fit_with_tuning_real(X,y)


genotype.pfile = "sample"
n <- 2000
p = 8000
vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd=paste0('zstdcat', ' ', paste0(genotype.pfile, '.pvar.zst'))), 'CHROM'='#CHROM'), VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID
pvar <- pgenlibr::NewPvar(paste0(genotype.pfile, '.pvar.zst'))
pgen <- pgenlibr::NewPgen(paste0(genotype.pfile, '.pgen'))
snp_matrix <- prepareFeatures(pgen, vars, vars)
float_matrix = ReadList(pgen, 1:p, meanimpute =TRUE)
psam.ids <- readIDsFromPsam(paste0(genotype.pfile, '.psam'))

matmul_result <- double(n*2)
intercepts = c(1.0,3.0)
x = snp_matrix
dense_beta = matrix(rnorm(2*p),nrow=p,ncol=2)
.Call("PlinkMultvC", x@ptr, x@xim, x@covs, n, p , x@ncov, dense_beta, intercepts, matmul_result)
matmul_result = matrix(matmul_result, nrow=n)

matmul_result2 = float_matrix %*% dense_beta
for(i in 1:(n*(1-val_ratio))){
  matmul_result2[i, ] <- matmul_result2[i, ] + intercepts
}

matmul_result - matmul_result2

beta = c(rep(1,10), rep(0, p-10))
y = 5+ float_matrix %*% beta + rnorm(n)#as.matrix(phe.master[,..covariates])%*% beta_c
start.time <- Sys.time()
res1 = myglmnet(float_matrix, y, family='quantile', standardize=F, intercept=T,nlambda=30, tau=0.5)
end.time <- Sys.time()
print(end.time-start.time)

start.time <- Sys.time()
res2 = myglmnet(snp_matrix, y, family='quantile', lambda = res1$lambda, standardize=F, intercept=T,nlambda=30, tau=0.5)
end.time <- Sys.time()
print(end.time-start.time)

best_fit <- fit_with_tuning_real(float_matrix,y)
