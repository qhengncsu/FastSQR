
#install.packages('/Users/ruilinli/myglmnet', repo=NULL,type='source')
library(pgenlibr)
library(myglmnet)
library(microbenchmark)

n = 1000
p = 2000
t0 = Sys.time()
m = PlinkMatrix("/Users/ruilinli/plink-ng/toy_data.pgen", 1:n, 1:p)
m = actualize(m)
t1 = Sys.time()
t1 - t0
#plinktest(m$pgen, rnorm(100))
pgen <- pgenlibr::NewPgen("/Users/ruilinli/plink-ng/toy_data.pgen", pvar = NULL, sample_subset =1:n)
m2 = ReadList(pgen, 1:p, meanimpute =TRUE)
v = rep(1.0, p)
eta = double(n)
mytest(m@ptr, m@xim, n, p, v,eta)
beta = rnorm(p) * rbinom(p,1, 0.3)
y = m2 %*% beta
y = y/sd(y)
microbenchmark(
  myglmnet(m, y, family='gaussian', standardize=F, intercept=F),
  myglmnet(m2, y, family='gaussian', standardize=F, intercept=F),
  glmnet(m2, y, family='gaussian', standardize=F, intercept=F),
  times=1L
)
a1 = myglmnet(m, y, family='gaussian', standardize=F, intercept=F)
a1 = myglmnet(m, y, family='gaussian', standardize=F, intercept=F)
a1 = myglmnet(m, y, family='gaussian', standardize=F, intercept=F)

y = m2 %*% beta
a2 = myglmnet(m2, y, family='gaussian', standardize=F, intercept=F)
a1 = glmnet(m2, y, family='gaussian', standardize=F, intercept=F)
y2 = rep(1.0, p)
y2[y<median(y)] = 0.0

a1 = myglmnet(m, y2, family='logistic', standardize=F, intercept=F)
a1 = myglmnet(m, y2, family='logistic', standardize=F, intercept=F)
a1 = myglmnet(m, y2, family='logistic', standardize=F, intercept=F)

a2 = myglmnet(m2, y2, family='logistic', standardize=F, intercept=F)



















stop('sufficient')
library(myglmnet)
library(glmnet)
library(microbenchmark)
n = 20
p = 5
 set.seed(1)

X = matrix(rnorm(n*p),n,p)
beta = rnorm(p) * rbinom(p, 1, 0.6)

y = X %*% beta
sdy = sd(y)
y = y /sdy

y2 = rep(1.0, n)
y2[y<median(y)] = 0.0
lam = c(0.3, 0.2, 0.15, 0.1, 0.05)

ref = glmnet(X, y2, family='binomial', lambda = lam, standardize = F, intercept = T)

test = myglmnet(X, y2, family='logistic', lambda = lam,standardize = F, intercept =T)

max(ref$beta[,10] - test$beta[,10])



ref = glmnet(X, y, family='gaussian', exclude=c(2L), standardize = F, intercept = T)

test = myglmnet(X, y, family='gaussian', exclude=c(2L), standardize = F, intercept = T)

## Performance profiling
library(microbenchmark)
library(myglmnet)
library(glmnet)
n = 1500
p = 3000
set.seed(1)
X = matrix(rnorm(n*p),n,p)
beta = rnorm(p) * rbinom(p, 1, 0.3)
y = X%*% beta
y = y/sd(y)
w = rep(1/n, n)
y2 = rep(1.0, n)
y2[y < median(y)] = 0.0
l1 = glmnet(X, y, family='gaussian',  weights = w)
l2 = myglmnet(X, y, family='gaussian')

a1 = myglmnet(X, y2, family='logistic', exclude=c(2L))
a2 = glmnet(X, y2, family='binomial', exclude=c(2L))

microbenchmark(
  myglmnet(X, y, family='gaussian', exclude=c(2L)),
  glmnet(X, y, family='gaussian', exclude=c(2L)),
  times=3L
)

# warm start test
beta0 = l1$beta[,50]
lam = l1$lambda[51:length(l1$lambda)]
a3 = myglmnet(X, y, family='gaussian', beta0 =beta0, lambda=lam)
a4 = glmnet(X, y, family='gaussian', lambda = lam)

beta0 = a2$beta[,50]
lam = a2$lambda[51:length(a2$lambda)]
a3 = myglmnet(X, y2, family='logistic', beta0 =beta0, lambda=lam)
a4 = glmnet(X, y2, family='binomial', lambda = lam)

microbenchmark(
  myglmnet(X, y2, family='logistic', exclude=c(2L)),
  glmnet(X, y2, family='binomial', exclude=c(2L)),
  times=3L
)






#
# ref = glmnet(X, y, family='gaussian', intercept = F, standardize = F, weights = w)
# lambdar = ref$lambda
# # ref = glmnet(X, y, family='gaussian', lambda=lambdar, intercept = F, standardize = F)
# # tset = wrapper(X, y,lambda=lambdar)
#
# microbenchmark(glmnet(X, y, family='gaussian', lambda=lambdar, intercept = F, standardize = F,weights = w),
#                wrapper(X, y,lambda=lambdar), times=1L)



# X = matrix(rnorm(n*p),n,p)
# beta = rnorm(p) * rbinom(p, 1, 0.6)

# y = X %*% beta
# sdy = sd(y)
# y = y /sdy


# ref = glmnet(X, y, family='gaussian', intercept = T, standardize = F)
# lambdar = ref$lambda
# wrapper(X, y,lambda=lambdar)



# microbenchmark(
#   glmnet(X, y, family='gaussian',lambda = ref$lambda, intercept = T, standardize = F),
#   wrapper(X, y,lambda=lambdar),
#   times = 1L
# )

# ref2 = glmnet(X, y, family='gaussian',lambda = ref$lambda, intercept = T, standardize = F)
#
# #
# # profvis({
# #   ref2 = glmnet(X, y, family='gaussian',lambda = ref$lambda, intercept = F, standardize = F)
# # })
#
# wrapper(X, y,lambda=lambdar)
#
#
#
# ref =glmnet(X,y, family='gaussian', lambda=lambdar,intercept = T, standardize = F)
# ref$beta
#
# # ref = glmnet()
#
# # wrapper()

## Cox test


rev_cumsum = function(v){
  rev(cumsum(rev(v)))
}

compute_lw = function(f, time, c){
  f = scale(f,TRUE,FALSE)
  theta = exp(f)[,1]
  n = length(theta)
  r = rank(time, ties.method="min")
  rm = rank(time, ties.method="max")
  rskden=rev_cumsum(theta)[r]
  lsecond = cumsum(c/rskden)[rm]
  lsecond = lsecond*theta
  wsecond = cumsum(c/rskden^2)[rm]
  wsecond = (exp(2*f)[,1])*wsecond
  l = (lsecond - c)/n
  w = (lsecond - wsecond)/n
  dev = mean((log(rskden) - f) * c)
  return(list(l = l, w=w, dev = dev*2))
}

library(myglmnet)
library(survival)
library(glmnet)

n = 5000
p = 5
X = matrix(rnorm(n*p),n,p)
beta = rnorm(p)
y = rexp(n) * exp(- X%*%beta)
y = round(y) + 1.5
status = rbinom(n,1,0.5)
s = Surv(y, status)
colnames(s) = c('time', 'status')
res = myglmnet(X, s, family="cox", standardize=F, nlambda=3)
res2 = glmnet(X, s, family="cox", standardize=F, nlambda=3)





res2 = coxph(Surv(y, status)~ X)
o = order(y)
v = runif(n)
v = v/sum(v)
eta = X %*% beta
lw = compute_lw(eta[o], y[o], v[o]*n)
time = y[o]
r = rank(time, ties.method="min")
rm = rank(time, ties.method="max")
rskden=rev_cumsum(exp(eta[o]))[r]

library(myglmnet)
mytest(eta, y, v)


cox_residual = function(f, time, c){
  f=scale(f,TRUE,FALSE)
  o = order(time)
  ordered_time = time[o]
  r = rank(ordered_time, ties.method="min")
  rm = rank(ordered_time, ties.method="max")
  ef=exp(f)[o]
  rskden=rev(cumsum(rev(ef)))[r]
  grad =  cumsum(c[o]/rskden)[rm]
  grad[o] = c[o] -  ef * grad
  grad
}


resid = cox_residual(rep(0, n), y, rep(1/n,n))
resid %*% X



cdeviance <-
    function (pred = NULL, y, x = 0, offset = NULL, weights = NULL,
              beta = NULL)
{
    storage.mode(x) = "double"
    ty = y[,1]
    tevent = y[,2]
    ty = ty + (1 - tevent) * 100 * .Machine$double.eps
    nobs = as.integer(length(ty))
    nvars = as.integer(ncol(x))
    nvec=1
    if (is.null(weights))
        weights = rep(1, nobs)
    else {
        weights=nobs*weights/sum(weights)
        weights = as.double(weights)
    }
### Compute saturated loglikelihood
    wd=weights[tevent==1]
    tyd=ty[tevent==1]
    if(any(duplicated(tyd))){
        wd=tapply(wd,tyd,sum)
    }
    wd=wd[wd>0]
    lsat=-sum(wd*log(wd))
####
    if (is.null(offset))
        offset = rep(0, nobs)
    else offset=as.double(offset)
    if (is.null(beta)) {
        beta = double(0)
        nvars = as.integer(0)
    }
    else {
        beta = as.matrix(beta)
        nvec = ncol(beta)
    }
    if(!is.null(pred)){
        # trick to get a set of deviances based on predictions"
        x=as.matrix(pred)
        nvec=ncol(x)
        storage.mode(x)="double"
        beta=diag(nvec)
        nvars=as.integer(nvec)
        storage.mode(beta)="double"
    }
    nvec=as.integer(nvec)


    lsat

}


# Code to compute saturated log partial likelihood
lsat = function(y, w=NULL)
{
  yt = y[,1]
  yevent = y[,2]
  n = length(yevent)
  has_event = yevent > 0
  yt = yt[has_event]
  nevent = length(yt)
  if(is.null(w)){
    w = rep(1,n)
  } else {
    w = w[has_event]
  }

  o = order(yt)
  yt = yt[o]
  w = w[o]

  r =  rank(yt, ties.method='min')
  rm = rank(yt, ties.method='max')
  result = 0
  i = 1
  while(i < nevent + 1){
    diff = rm[i] - r[i]
    if(diff > 0){
      wsum = 0
      for(j in 1:(diff+1)){
        wj = w[i+j-1]
        if(wj != 0){
          result = result + wj * log(wj)
          wsum = wsum + wj
        }
      }
      if(wsum != 0){
        result = result - wsum * log(wsum)
      }
    }
    i = i + diff + 1
  }
  result
}