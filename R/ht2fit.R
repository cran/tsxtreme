## Copyright (C) 2017 Thomas Lugrin
## Gathering functions from old files into one reference file
## (Scope:) Heffernan-Tawn model fit in 2 stages
## List of functions: - ht2fit
##                    - nllh.ht
##                    - grad.ht
##                    - ht2step.2d
##                    - ht2step
##################################################

## >   ts: vector of reals, time series
## >   u.mar,u.dep: scalars, marginal/dependence thresholds given as a probabilities
## >   nlag: integer, number of lags to consider
## >   conditions: logical, TRUE means alpha and beta must satisfy AD/AI constraints
## <   ret: list, MLE for alpha - beta - z hat - standard deviations of alpha and beta
## .   called by user
dep2fit <- function(ts, u.mar=0, u.dep, lapl=FALSE, method.mar=c("mle","mom","pwm"), nlag=1, conditions=TRUE){
  data.up <- format.ts(ts=ts, u.mar=u.mar, u.dep=u.dep, method=method.mar,
                       lapl=lapl, nlag=nlag)
  ht2step(data=data.up, conditions=conditions)
}

##################################################
## LIKELIHOODS

## >   par: vector of reals, alpha-beta-mu-sigma^2
## >   data: 2-column matrix, first column is X>u - second column is Y|X>u
## >   conditions: boolean, TRUE means alpha and beta must satisfy AD/AI constraints
## <   ret: scalar, value of the negative log-likelihood of H+T
## .   called by ht2step.2d
## ... find estimates of Heffernan-Tawn's alpha and beta

nllh.ht <- function(par,data,conditions){
  if(par[4] <= 0){ return(Inf) }
  if(conditions){
    if(!conditions.verify(par[1], par[2], 0, data=data) ||
         !conditions.verify(par[1], par[2], 1, data=data)){ return(Inf) }
  }else{
    if(par[1]< -1 || par[1]>1 || par[2]<0 || par[2]>1){ return(Inf) }
  }
  sigma <- par[4]*data[,1]^(2*par[2])
  mu    <- par[1]*data[,1] + par[3]*data[,1]^par[2]
  return(sum(log(sigma) + (data[,2]-mu)^2/sigma))
}


## cf. above - > conditions not used
grad.ht <- function(par,data,conditions){
  sig <- par[4]*data[,1]^(2*par[2])
  mu  <- par[1]*data[,1] + par[3]*data[,1]^par[2]
  cen <- data[,2]-mu
  
  d.a <- -2*sum(cen*data[,1]/sig)
  d.b <- 2*sum(log(data[,1]) * (1 - cen*(par[3]*data[,1]^par[2]+cen)/sig))
  d.m <- -2*sum(cen/(data[,1]^(par[2])*par[4]))
  d.s <- sum((1-cen^2/sig)/par[4])
  
  return(c(d.a,d.b,d.m,d.s))
}


##################################################
## FITS

## >   data: 2-column matrix, first column is X>u (Laplace) - second column is Y|X>u
## >   conditions: boolean, TRUE means alpha and beta must satisfy AD/AI constraints
## <   ret: list, MLE for alpha - beta - z hat - standard deviations of alpha and beta
## .   called by ht2step or user
ht2step.2d <- function(data, conditions=TRUE){
  ret <- stepfit()
  if(conditions){
    bds <- conditions.bounds(0,FALSE,data)
    a   <- runif(1,bds[1],bds[2])
    bds <- conditions.bounds(a,TRUE,data)
    b   <- runif(1,bds[1],bds[2])
  }else{
    a <- runif(1, -1, 1)
    b <- runif(1, 0, 1)
  }
  res <- optim(par=c(a,b,0,1), fn=nllh.ht, gr=grad.ht, data=data, conditions=conditions,
                hessian=TRUE, method="BFGS")
  ret$a <- res$par[1]
  ret$b <- res$par[2]
  ret$res <- (data[,2] - data[,1]*ret$a)/data[,1]^ret$b
  ret$pars.se <- sqrt(diag(solve(res$hessian))[1:2])
  return(ret)
}



## >   data: matrix of reals, first column is X>u (Laplace) - other columns are |X>u
## >   conditions: boolean, TRUE means alpha and beta must satisfy AD/AI constraints
## <   ret: list, MLE for alpha - beta - z hat - standard deviations of alpha and beta
## .   called by th2fit
ht2step <- function(data, conditions=TRUE){
  n       <- dim(data)[1]
  d       <- dim(data)[2]
  ret   <- stepfit()
  ret$a <- ret$b <- numeric(d-1)
  ret$res     <- matrix(0, nrow=n, ncol=d-1, dimnames=list(NULL, paste("Z",1:(d-1), sep="")))
  ret$pars.se <- matrix(NA, nrow=2, ncol=d-1)
  ret$nlag    <- d-1
  for(i in 2:d){
    fit           <- ht2step.2d(data[,c(1,i)], conditions)
    ret$a[i-1]    <- fit$a
    ret$b[i-1]    <- fit$b
    ret$res[,i-1]     <- fit$res
    ret$pars.se[,i-1] <- fit$pars.se
  }
  return(ret)
}