## Copyright (C) 2017 Thomas Lugrin
## runs estimator for theta(x,m)
## (Scope:) empirical estimate to be compared to E+T typically
## List of functions: - runs
##                    - compute.runs
##                    - block.bootstrap.sample
##################################################

##################################################
## RUNS ESTIMATOR

## >   ts: vector of reals, time series which to estimate theta(x,m) for
## >   nlag: integer, nbr of lags (run-length - 1)
## >   u.mar: probability, threshold used for marginal transform
## >   probs: vector of probabilities, x in theta(x,m)
## >   block.length: integer>0, used for the block-bootstrap CI
## >   R.boot: integer>0, nbr of repetitions used for the block-bootstrap CI - 0 means CI not computed
## >   levels: vector of probabilities, specifies which posterior quantiles to compute (+estimate)
## <   ret: list, theta with CI - nbr of exc. across [mesh] - [mesh] in Laplace/original scale
## .   called by user
thetaruns <- function(ts, lapl=FALSE, nlag=1, u.mar=0,
                      probs=seq(u.mar,0.995,length.out=30),
                      method.mar=c("mle","mom","pwm"),
                      R.boot=0, block.length=(nlag+1)*5, levels=c(.025,.975)){
  ## marginal re-scale to Laplace if needed
  if(!lapl)
    ts.L <- scale.ts(ts=ts, u=u.mar, method=method.mar)$ts.L
  data <- lags.matrix(ts=ts.L, nlag=nlag)
  mesh <- probs
  mesh.O <- scale.to.original(p=mesh, ts=ts, u=u.mar, gpd.pars=gpd(ts, u=u.mar, method=method.mar)$pars)
  mesh.L <- qexp(1-(1-mesh)*2)
  ret <- depmeasure("runs")
  ret$probs <- probs
  ret$levels <- mesh.O
  ret$nlag <- nlag
  ## set containers
  nbr.vert <- length(mesh)
  nbr.quant<- 1+length(levels)
  theta    <- matrix(0, nrow=nbr.vert, ncol=nbr.quant)
  colnames(theta) <- c("estimate",paste(levels*100,"%", sep=""))
  ## compute and store
  th.est             <- compute.runs(data,mesh.L)
  theta[,"estimate"] <- th.est[[1]]
  nbr.exc            <- th.est[[2]]
  ret$nbr.exc <- nbr.exc
  ## compute bootstrapped quantiles
  if(R.boot > 0 & nbr.quant > 1){
    ts.R <- block.bootstrap.sample(ts, block.length, R.boot)
    th.R <- matrix(0, nrow=nbr.vert, ncol=R.boot)
    for(r in 1:R.boot){
      ts.boot  <- scale.ts(ts.R[,r], u=u.mar, method=method.mar)$ts.L
      ts.mat   <- lags.matrix(ts.boot, nlag=nlag)
      th.R[,r] <- compute.runs(ts.mat,mesh.L)[[1]]
    }
    theta[,-1] <- t(apply(th.R, 1, quantile, levels))
  }else{
    theta <- theta[,"estimate",drop=FALSE]
  }
  ret$theta <- theta
  return(ret)
}



## >   data: matrix, lags matrix of a time series
## >   mesh: vector, x in theta(x,m) on Laplace scale
## <   ret: list, estimates of theta across [mesh] - nbr of exceedances across [mesh]
## .   called by thetaruns
compute.runs <- function(data, mesh){
  nbr.vert <- length(mesh)
  theta <- nbr.exc <- numeric(nbr.vert)
  for(i in 1:nbr.vert){
    exc        <- data[data[,1]>=mesh[i],,drop=FALSE]
    nbr.exc[i] <- dim(exc)[1]
    theta[i]   <- sum(apply(exc[,-1,drop=FALSE], 1, max)<mesh[i])/nbr.exc[i]
  }
  return(list(theta,nbr.exc))
}


## >   ts: vector of reals, time series which to estimate theta(x,m) for
## >   block.length: integer>0, used for the block-bootstrap CI
## >   R: integer>0, nbr of repetitions used for the block-bootstrap CI
## <   matrix, each column is a block-bootstrapped time series
## .   called by thetaruns()
block.bootstrap.sample <- function(ts, block.length, R){
  ## handle cases when block length does not divide the time series
  rest <- length(ts)%%block.length
  if(rest > 0)
    ts <- ts[-((length(ts)-rest+1):length(ts))]
  ts.mat <- matrix(ts, nrow=block.length, byrow=FALSE)
  nbr.bl <- floor(length(ts)/block.length)# == dim(ts.mat)[2]
  ts.R   <- matrix(0, nrow=length(ts), ncol=R)
  for(r in 1:R){
    sample   <- sample.int(nbr.bl, size=nbr.bl, replace=TRUE)
    ts.R[,r] <- ts.mat[,sample]
  }
  return(ts.R)
}