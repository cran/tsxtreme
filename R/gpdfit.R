## Copyright (C) 2017 Thomas Lugrin
## Generalised Pareto distribution fit
## (Scope:) GPD fit, marginal transform
## List of functions: - nllh.gpd
##                    - grad.gpd
##                    - nllh.gpd.0
##                    - gpd.mle
##                    - gpd.mom
##                    - gpd.pwm
##################################################

##################################################
## LIKELIHOOD-BASED GPD FIT
## >   par: vector of 2 reals, scale-shape
## >   u: scalar, threshold given as a quantile
## >   x: vector of reals, data on which to fit the GPD
## <   ret: scalar, value of the negative log-likelihood of the GPD
## .   called by gpd.mle
nllh.gpd <- function(par, u, x){
  exced <- x[x>u]-u
  n     <- length(exced)
  if(par[1] <= 0) return(Inf)
  if(abs(par[2]) < 0.00001){
    nllh     <- n*log(par[1])+sum(exced)/par[1]
  }
  else{
    pos.part <- pmax(1+par[2]/par[1]*exced, 0)
    if(min(pos.part) <= 0) return(Inf)
    nllh     <- n*log(par[1])+(1/par[2]+1)*sum(log(pos.part))
  }
  return(nllh)
}

grad.gpd <- function(par, u, x){
  exced <- x[x>u]-u
  n     <- length(exced)
  if(abs(par[2]) < 0.00001){
    d.s <- n/par[1] - sum(exced)/par[1]^2
    d.x <- 0
  }
  else{
    d.s <- ( n-(1/par[2]+1)*sum(1/(1+par[1]/(par[2]*exced))) )/par[1]
    d.x <- -1/par[2]^2*sum(log(1+par[2]*exced/par[1]) + (1/par[2]+1)*sum(1/(par[1]/exced+par[2])))
  }
  return(c(d.s,d.x))
}


## >   par: scalar, scale (shape is assumed 0)
## >   u: scalar, threshold given as a quantile
## >   x: vector of reals, data on which to fit the GPD
## <   ret: scalar, value of the negative log-likelihood of the GPD
## .   called by gpd.mle when shape.null.hyp=TRUE
nllh.gpd.0 <- function(par, u, x){
  nllh.gpd(c(par,0), u, x)
}



## >   ts: vector of reals, stationary time series
## >   u: scalar, threshold given as a quantile
## >   hessian: boolean, TRUE means computing the hessian matrix of the likelihood
## <   ret: list, threshold - MLEs - log-likelihood - standard deviations - covariance matrix
## .   called by scale.ts
## ... find MLEs of scale and shape of a GPD
gpd.mle <- function(ts, u, hessian=FALSE){
  opt <- optim(c(sd(ts),0.1), fn=nllh.gpd, gr=grad.gpd, u=u, x=ts, hessian=hessian)
  ret <- list(u=u, pars=opt$par, llh=-opt$value)
  if(hessian){
    ret$var     <- solve(opt$hessian)
    ret$pars.sd <- sqrt(diag(var))
  }
  return(ret)
}

##################################################
## MOMENT-BASED GPD FIT
## >   ts: vector of reals, stationary time series
## >   u: scalar, threshold given as a quantile
## <   ret: list, threshold - MOM estimates
## .   called by scale.ts
## ... find MOM estimates of scale and shape of a GPD
gpd.mom <- function(ts, u){
  ts <- ts-u
  m <- mean(ts)
  s <- var(ts)
  scale <- m*(m^2/s+1)/2
  shape <- (m^2/s-1)/2
  return(list(u=u, pars=c(scale,shape)))
}

gpd.pwm <- function(ts, u){
  ts <- ts-u
  ts <- sort(ts, decreasing=FALSE)
  n  <- length(ts)
  a0 <- mean(ts)
  a1 <- sum(ts[-n]*((n-1):1))/(n-1)/n
  scale <- 2*a0*a1/(a0-2*a1)
  shape <- a0/(a0-2*a1)-2
  return(list(u=u, pars=c(scale,shape)))
}


##################################################
## GENERIC FUNCTION
gpd <- function(ts, u, method){
  fct.name <- paste("gpd.",method[1], sep="")
  eval(call(fct.name, ts=ts, u=u))
}