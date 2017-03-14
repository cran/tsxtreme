## Copyright (C) 2017 Thomas Lugrin
## Marginal transforms
## (Scope:) H+T/E+T fits, marginal transformations, thresholding
## List of functions: - format.ts
##                    - scale.ts
##                    - scale.to.margin
##                    - threshold.margin
##                    - lags.matrix
##################################################


##################################################
## FORMAT TIME SERIES FOR H+T MODEL

## >   ts: vector of reals, (stationary!) time series
## >   u.mar: scalar, marginal (GPD) threshold given as a probability
## >   u.dep: scalar, conditional threshold given as a probability
## >   method: string, "mle" for max likelihood or "mom" for the method of moments or "pwm" for proba weighted moments
## >   nlag: integer, number of lags to be considered
## >   lapl: logical, is [ts] on Laplace scale already?
## <   ret: matrix with [nlag]+1 columns on the Laplace scale, with 1st column above threshold
## .   called by [etfit()], [et2fit()]
format.ts <- function(ts, u.mar, u.dep, method=c("mle","mom","pwm"), nlag, lapl=FALSE){
  if(!lapl)
    ts <- scale.ts(ts=ts, u=u.mar, method=method)$ts.L
  data <- lags.matrix(ts=ts, nlag=nlag)
  data <- threshold.margin(data=data, u=u.dep, quantile=FALSE)
  return(data)
}


##################################################
## MARGINAL TRANSFORMATIONS

## >   ts: vector of reals, (stationary!) time series
## >   u: scalar, threshold given as a probability
## >   method: string, "mle" for max likelihood or "mom" for the method of moments or "pwm" for proba weighted moments
## <   ret: list, original time series with different marginal distributions + GPD parameters
## .   called by [format.ts()]

scale.ts <- function(ts, u, method=c("mle","mom","pwm")){
  u    <- quantile(ts, u)
  n    <- length(ts)
  ts.U <- numeric(n)
  ei   <- (ts>u)
  p.u  <- sum(ei)/n
  pars <- gpd(ts, u, method)$pars
  if(pars[2] == 0){
    ts.U[ei] <- 1 - p.u*exp(-(ts[ei]-u)/pars[1])
  }
  else{
    ts.U[ei]  <- 1 - p.u*(1+pars[2]/pars[1]*(ts[ei]-u))^(-1/pars[2])
  }
  ts.U[!ei] <- rank(ts[!ei])/(n+1)
  ts.L <- numeric(n)
  ts.L[ts.U<0.5]  <- log(2*ts.U[ts.U<0.5])
  ts.L[ts.U>=0.5] <- -log(2*(1-ts.U[ts.U>=0.5]))
  return(list(ts.U=ts.U, ts.G=-log(-log(ts.U)), ts.F=-1/log(ts.U), ts.L=ts.L, pars=c(pars,p.u)))
}



## >   x: (vector of) reals, (a range of) observations for which change of scale is needed
## >   ts: vector of reals, time series in original scale
## >   u: scalar, threshold given as a probability
## >   gpd.pars: vector of 3 reals, scale-shape-probability of exceeding u
## <   ret: list, original observation(s) scaled to different marginal distributions
## .   called by ?
scale.to.margin <- function(x, ts, u, gpd.pars){
  u   <- quantile(ts,u)
  n   <- length(ts)
  x.U <- x
  if(gpd.pars[2] == 0){
    x.U[x>u] <- 1 - gpd.pars[3]*exp(-(x[x>u]-u)/gpd.pars[1])
  }
  else{
    x.U[x>u] <- 1 - gpd.pars[3]*(1+gpd.pars[2]/gpd.pars[1]*(x[x>u]-u))^(-1/gpd.pars[2])
  }
  x.U[x<=u] <- approx(sort(ts), y=(1:n)/(n+1), xout=x.U[x<=u], method="linear", ties="ordered", rule=2)$y
  if(x.U<0.5) x.L <- log(2*x.U)
  else x.L <- -log(2*(1-x.U))
  return(list(x.U=x.U, x.G=-log(-log(x.U)), x.F=-1/log(x.U), x.L=x.L))
}



## >   p: (vector of) reals, (a range of) probabilities which to transform to original [ts] scale
## >   ts: vector of reals, time series in original scale
## >   u: scalar, threshold given as a probability
## >   gpd.pars: vector of 2 reals, scale-shape parameters of the GPD fitted to [ts]
## <   a vector of the same length as [p] in [ts] scale
## .   called by thetafit, chifit, theta2fit
scale.to.original <- function(p, ts, u, gpd.pars){
  u.O <- quantile(ts, u)
  exc <- p > u
  p[exc]  <- qgpd((p[exc]-u)/(1-u), loc=u.O, scale=gpd.pars[1], shape=gpd.pars[2])
  p[!exc] <- quantile(ts, p[!exc])
  return(p)
}



##################################################
## PRE-PROCESSING FOR HT & ET

## >   data: matrix of reals, replicates in rows - variables in columns
## >   u: scalar, threshold given as a probability or quantile (cf. quantile)
## >   quantile: boolean, TRUE means u is given as a quantile - FALSE means as a probability
## <   ret: matrix of reals, only first column observations exceeding u are kept
## .   called by format.ts
threshold.margin <- function(data, u=0.95, quantile=FALSE)
{
  if(!quantile) u <- quantile(data[,1], u)
  ei      <- (data[,1]>u)
  return(data[ei,])
}



## >   ts: vector of reals, (stationary!) time series
## >   nlag: scalar, number of lags of the output lagged matrix (=> nbr of columns == nlag+1)
## <   ret: matrix of lagged versions of ts in its columns
## .   called by format.ts and thetaruns
lags.matrix <- function(ts, nlag){
  n      <- length(ts)
  dim    <- nlag+1
  ts.mat <- matrix(0, nrow=n-dim+1, ncol=dim)
  for(i in 1:dim){
    ts.mat[,i] <- ts[i:(n-dim+i)]
  }
  return(ts.mat)
}
