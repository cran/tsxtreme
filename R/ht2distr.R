## Copyright (C) 2017 Thomas Lugrin
## Functions related to H+T model
## (Scope:) Heffernan-Tawn model fit in 2 stages
## List of functions: - p.res2
##                    - q.res2
##################################################

## >   res: vector of reals, quantiles on which to compute the distribution function of the residuals
## >   sorted.res: vector, EDF
## <   ret: vector of the same length as z, distribution of the residuals evaluated in z
## .   called by et.2.step

p.res2 <- function(res, sorted.res){
  n <- length(sorted.res)
  a <- approx(sorted.res, y=(1:n)/(n+1), xout=res, method="linear", ties="ordered", rule=2)$y
  return(a)
}



## >   p: scalar, probability (to compute quantiles of residual distribution), in [0,1]
## >   a: scalar, alpha parameter, in [-1,1]
## >   b: scalar, beta parameter, in [0,1]
## >   data: bivariate vector, (X,Y) with Y | X>u
## <   boolean, does the couple (a,b) satisfy the conditions?
## .   called by verifies.conditions to get quantiles of the EDF of the residual distribution given (a,b)
q.res2 <- function(p, a, b, data){
  if(dim(data)[2] != 2) stop("data are in the wrong format. Provide a 2-column matrix for (X,y), with Y|X>u")
  Z <- (data[,2] - a*data[,1])/data[,1]^b
  return(quantile(Z, p))
}