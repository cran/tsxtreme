## Copyright (C) 2017 Thomas Lugrin
## Implementation of current density/distribution/etc. functions
## (Scope:) Any
## List of functions: - dlapl
##                    - plapl
##                    - qlapl
##                    - rlapl
##################################################

##################################################
## LAPLACE DISTRIBUTION
dlapl <- function(x, loc=0, scale=1, log=FALSE){
  if(any(scale < 0)) stop("scale parameter must be non-negative")
  
  if(log){
    ret <- -abs(x-loc)/scale-log(2)-log(scale)
    ret[is.nan(ret)] <- 1
    return(ret)
  }else{
    ret <- exp(-abs(x-loc)/scale)/2/scale
    ret[is.nan(ret)] <- 1
    return(ret)
  }
}

plapl <- function(q, loc=0, scale=1, lower.tail=TRUE, log.p=FALSE){
  if(any(scale < 0)) stop("scale parameter must be non-negative")
  ret <- sign(q-loc)*(1-exp(-abs(q-loc)/scale))/2
  ret[is.nan(ret)] <- 0.5
  if(lower.tail){
    if(log.p)
      return(log(0.5+ret))
    else
      return(0.5+ret)
  }else{
    if(log.p)
      return(log(0.5-ret))
    else
      return(0.5-ret)
  }
}

qlapl <- function(p, loc=0, scale=1, lower.tail=TRUE, log.p=FALSE){
  if(any(p < 0) || any(p > 1)) stop("p must lie between 0 and 1")
  if(any(scale < 0)) stop("scale must be non-negative")
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
  ret <- loc-scale*sign(p-0.5)*log(1-2*abs(p-0.5))
  nas <- is.nan(ret)
  if(sum(nas)){
    p <- rep(p, length.out=length(ret))
    p <- p[nas]
    ret[nas] <- ifelse(p==1,Inf,-Inf)
  }
  return(ret)
}

rlapl <- function(n, loc=0, scale=1){
  if(any(scale < 0)) stop("scale must be non-negative")
  u <- runif(n)
  return(qlapl(u, loc=loc, scale=scale))
}