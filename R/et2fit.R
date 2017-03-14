## Copyright (C) 2017 Thomas Lugrin
## user interface for E+T fit
## (Scope:) Model fit in 2 stages
## List of functions: - theta2fit
##                    - th2est
##                    - thboot
##                    - blockboot
##################################################

##################################################
## FIT

## >   ts: vector of reals, time series which to estimate theta(x,m) from
## >   lapl: boolean, is [ts] already with Laplace margins?
## >   m: integer>2, run-length
## >   R: integer>0, nbr of samples used for the estimation procedure
## >   u.mar: probability, threshold used for marginal threshold - 0 indicates data are already in Laplace margins
## >   u.dep: probability, threshold used for Heffernan-Tawn model
## >   mesh: vector of probabilities, x in theta(x,m)
## >   method.mar: string, "mle" for max likelihood or "mom" for the method of moments or "pwm" for proba weighted moments
## >   method: (vector of) string(s), either "prop" or "MCi"
## >   silent: boolean, verbosity
## >   R.boot: integer>0, nbr of bootstrap samples to compute CI - default is 0 (do not compute)
## >   block.length: integer>0, length of blocks in block-bootstrap
## <   ret: list, ht fit - mesh in Laplace scale - theta - ratio of thetas
## .   called by user

#ONLY UPPER TAIL!
theta2fit <- function(ts, lapl=FALSE, nlag=1, R=1000,
                         u.mar=0, u.dep, probs=seq(u.dep,0.9999,length.out=30),
                         method.mar=c("mle","mom","pwm"), method=c("prop","MCi"),
                         silent=FALSE,
                         R.boot=0, block.length=m*5, levels=c(.025,.975)){
  data.up <- format.ts(ts=ts, u.mar=u.mar, u.dep=u.dep, method=method.mar,
                       nlag=nlag, lapl=lapl)
  m       <- nlag+1
  n       <- dim(data.up)[1]
  mesh    <- probs
  mesh.O  <- scale.to.original(p=mesh, ts=ts, u=u.mar, gpd.pars=gpd(ts, u=u.mar, method=method.mar)$pars)
  
  fit <- ht2step(data.up)
  if(!silent)
    print(paste("Likelihood fit just ended, with estimates for alpha, beta resp.:",
                signif(fit$a,3),"and",signif(fit$b,3)))
  n.vert <- length(mesh)
  sim.U  <- runif(R)
  ind    <- sample.int(n, size=R, replace=TRUE)
  sim.Z  <- fit$res[ind,, drop=FALSE]
  #compute theta for each m at each probability in mesh
  nbr.quant <- length(levels)+1
  th        <- matrix(NaN, nrow=n.vert, ncol=nbr.quant,
                    dimnames=list(NULL,c("estimate",paste(levels*100,"%",sep=""))))
  th[,"estimate"] <- th2est(fit, sim.U=sim.U, sim.Z=sim.Z, mesh=mesh, nlag=nlag, method=method)
  ret <- depmeasure("steptheta")
  ret$fit <- fit; ret$probs <- probs; ret$levels <- mesh.O; ret$nlag <- nlag
  ## compute bootstrap intervals (recursive)
  if(R.boot > 0 & nbr.quant > 1){
    boot <- thboot(par=list(ts=ts,u.mar=u.mar,u.dep=u.dep,nlag=nlag,method.mar=method.mar,R=R,mesh=mesh),
                      block.length=block.length, R=R.boot, method=method, levels=levels)
    n.vert <- dim(boot)[2]
    th[1:n.vert,2:nbr.quant] <- boot
  }else{
    th <- th[,"estimate",drop=FALSE]
  }
  ret$theta <- th
  return(ret)
}


## >   fit: list, output from ht.2.step
## >   sim.U: vector of reals, R uniform samples
## >   sim.Z: matrix of reals, R multivariate samples from residual distribution
## >   mesh: vector of reals, quantile values on which to evaluate theta (uniform scale)
## >   method: (vector of) string(s), either "prop" or "MCi"
## <   ret: estimate of theta using the method of proportions or MC integration
## .   called by theta2fit
## ... estimates which are not calculated are returned as arrays of NaNs or zeroes
th2est <- function(fit, sim.U, sim.Z, mesh, nlag, method){
  n.vert <- length(mesh)
  R      <- length(sim.U)
  mesh.L <- qexp(1-(1-mesh)*2)
  nbr.quant <- length(levels)+1
  theta <- numeric(n.vert)
  # prepare MC integration...
  HZ.u   <- numeric(R)
  sort.Z <- apply(fit$res, 2, sort)
  for(i in 1:n.vert){
    sim.L <- -log(2*(1-mesh[i]-(1-mesh[i])*sim.U))
    # method of proportion
    if(sum(grepl("prop",method))){
      sim.Y <- as.matrix(sim.L%*%t(fit$a) + exp(log(sim.L)%*%t(fit$b))*sim.Z)# [Rx(m-1)]
      theta[i] <- sum(apply(as.matrix(sim.Y), 1, max)<mesh.L[i])/R
    }# Monte Carlo integration
    else if(sum(grepl("MCi",method))){
      sim.Z.MC <- (mesh.L[i]-sim.L%*%t(fit$a))/exp(log(sim.L)%*%t(fit$b))# [Rx(m-1)]
      HZ <- 1
      for(lag in 1:nlag){
        HZ <- HZ * p.res2(sim.Z.MC[,lag], sort.Z[,lag])
        if(i == 1){ HZ.u[lag,] <- HZ }
      }
      theta[i] <- mean(HZ)
    }else{
      stop("unrecognised method. Should be either 'prop' or 'MCi'")
    }
  }
  return(theta)
}


## >   par: list, arguments to feed theta2fit
## >   block.length: integer>0, length of blocks in block-bootstrap
## >   R: integer>0, nbr of bootstrap samples to compute CI
## >   method: string, either "prop" or "MCi"
## >   levels: vector of probabilities, quantiles of bootstrap distribution to be output
## <   matrix where columns are quantiles and rows are x-levels of theta(x,m)
## .   bootstrap CI for quantile estimates from 2-step E+T
thboot <- function(par, block.length, R, method, levels){
  nbr.vert <- length(par$mesh)
  ts.boot  <- blockboot(par$ts, block.length, R)
  n        <- block.length * floor(length(par$ts)/block.length)
  theta.R  <- apply(ts.boot, 2, function(ts,u.mar,u.dep,nlag,R,mesh,method.mar,method,shape.null.hyp,silent){
                                fit=theta2fit(ts=ts,u.mar=u.mar,u.dep=u.dep,
                                              lapl=TRUE,nlag=nlag,R=R,probs=mesh,
                                              method.mar=method.mar,method=method,silent=silent)
                                return(fit$theta[,"estimate"])}, 
                  nlag=par$nlag, R=par$R, u.mar=par$u.mar, u.dep=par$u.dep, mesh=par$mesh,
                  method.mar=par$method.mar, method=method, silent=TRUE)
  theta.R <- matrix(theta.R, nrow=nbr.vert, ncol=R)
  # remove NA values
  v <- 1
  while(v <= nbr.vert && sum(is.na(theta.R[v,])) == 0){ v <- v+1 }
  if(v == 1){
    warning("unable to compute uncertainty of theta")
    v <- 2
  }
  nbr.vert <- v-1
  theta.R  <- theta.R[1:nbr.vert,, drop=FALSE]
  # compute CI
  ci <- t(apply(theta.R, 1, quantile, levels))
  return(ci)
}


## >   ts: vector of reals, original time series
## >   block.length: integer>0, length of blocks for block bootstrap
## >   R: integer>0, number of bootstrap samples
## <   ret: matrix, R columns containing the block bootstrap samples of ts
## .   called by thboot
blockboot <- function(ts, block.length, R){
  rest <- length(ts)%%block.length
  if(rest > 0){ ts <- ts[-((length(ts)-rest+1):length(ts))] }
  ts.mat <- matrix(ts, nrow=block.length)
  nbr.bl <- floor(length(ts)/block.length)# == dim(ts.mat)[2]
  ts.R   <- matrix(0, nrow=nbr.bl*block.length, ncol=R)
  for(i in 1:R){
    sample   <- sample.int(nbr.bl, size=nbr.bl, replace=TRUE)
    ts.R[,i] <- c(ts.mat[,sample])
  }
  return(ts.R)
}