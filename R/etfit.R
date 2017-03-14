## Copyright (C) 2017 Thomas Lugrin
## user inteface for etfit class
## (Scope:) Fit Eastawn model using Bayesian semiparametrics
## List of functions: - thetafit
##                    - chifit
##################################################

##################################################
## E+T R interface

## >   ts: vector of reals, time series which to estimate theta(x,m) for
## >   lapl: boolean, is the time series already transformed to Laplace margins?
## >   m: integer>2, run-length
## >   R: integer>0, nbr of samples used in a single sweep for the estimation
## >   S: integer>0, nbr of posterior sweeps used for the estimation
## >   u.mar: probability, threshold used for marginal threshold
## >   u.dep: probability, threshold used for Heffernan-Tawn model
## >   probs: vector of probabilities, x in theta(x,m)
## >   method.mar: string, "mle" for max likelihood or "mom" for the method of moments or "pwm" for proba weighted moments
## >   method: (vector of) string(s), either "prop" or "MCi"
## >   silent: boolean, verbosity
## >   fit: boolean or list returned from htfit, TRUE means the htfit must be called
## >   par: list, contains all arguments needed for the call to htfit - if empty, default values are assumed
## >   submodel: string, structure imposed on (a,b) - defaults to "fom"
## >   levels: vector of probabilities, specifies which posterior quantiles to compute (+mean+median)
## <   ret: list, ht fit - probs in Laplace scale - theta - theta posterior samples
## .   called by user

thetafit <- function(ts, lapl=FALSE, nlag=1, R=1000, S=500,
                  u.mar=0, u.dep, probs=seq(u.dep,0.9999,length.out=30),
                  method.mar=c("mle","mom","pwm"), method=c("prop","MCi"),
                  silent=FALSE,
                  fit=TRUE, prev.fit=bayesfit(), par=bayesparams(),
                  submodel=c("fom","none"), levels=c(.025,.975)){
  data.up <- format.ts(ts=ts, u.mar=u.mar, u.dep=u.dep, method=method.mar,
                       nlag=nlag, lapl=lapl)
  ret <- etfit(data=data.up, R=R, S=S, probs=probs, method=method, silent=silent,
              fit=fit, prev.fit=prev.fit, par=par, submodel=submodel, levels=levels)
  mesh.O     <- scale.to.original(p=probs, ts=ts, u=u.mar, gpd.pars=gpd(ts, u=u.mar, method=method.mar)$pars)
  ret$levels <- mesh.O
  return(ret)
}


etfit <- function(data, R, S, probs, method,
                  silent,
                  fit, prev.fit, par, submodel, levels){
  data.up <- data
  n       <- dim(data.up)[1]
  mesh    <- probs
  
  if(fit){
    if(!is.bayesparams(par)) stop("par must be of class 'bayesparams'.")
    fit <- htfit(data.up, prop.a=par$prop.a, prop.b=par$prop.b,
                prior.mu=par$prior.mu, prior.nu=par$prior.nu, prior.eta=par$prior.eta,
                trunc=par$trunc, comp.saved=par$comp.saved, maxit=par$maxit,
                burn=par$burn, thin=par$thin,
                adapt=par$adapt, batch.size=par$batch.size, mode=par$mode, submodel=submodel)
    if(!silent)
      print(paste("Fit of DDP done. Mean value for alpha & beta:",
                  signif(apply(fit$a, 2, mean),3),"and",
                  signif(apply(fit$b, 2, mean),3)))
  }else{
    if(!is.bayesfit(prev.fit)) stop("prev.fit must be of class 'bayesfit'.")
    if(prev.fit$len == 0) stop("prev.fit must contain a proper instance of class 'bayesfit'.")
    fit <- prev.fit
  }
  # set up parameters for theta
  tic()
  nit      <- dim(fit$a)[1]
  nlag     <- dim(fit$a)[2]
  mesh.L   <- qexp(1-(1-mesh)*2)
  nbr.vert <- length(mesh.L)
  sim.U    <- runif(R)
  it <- sample.int(nit, S, replace=TRUE)
  a  <- as.double(fit$a[it,])
  b  <- as.double(fit$b[it,])
  m  <- fit$mean[it,,, drop=FALSE]
  s  <- fit$sd[it,,, drop=FALSE]
  w  <- fit$w[it,, drop=FALSE]
  # prepare returned value
  ret       <- depmeasure("theta")
  ret$fit   <- fit
  ret$probs <- mesh
  ret$nlag  <- nlag
  # method of proportions: generate residuals
  if(grepl("prop",method[1])){
    sim.Z <- r.res(R=R, S=S, nlag=nlag, w=w, m=m, s=s)
  }
  nbr.quant <- length(levels)+2
  th    <- matrix(0, nrow=nbr.vert, ncol=nbr.quant,
                 dimnames=list(NULL,c("mean","median",paste(levels*100,"%",sep=""))))
  distr <- matrix(0, nrow=nbr.vert, ncol=S) # compute coverage/RMSE when true theta is known
  for(i in 1:nbr.vert){
    sim.L <- -log(2*(1-mesh[i]-(1-mesh[i])*sim.U))
    ## method of proportions
    if(grepl("prop", method[1])){
      sim.Y <- sim.L%*%t(a) + exp(log(sim.L)%*%t(b))*sim.Z# [RxS*(m-1)]
      sim.Y <- matrix(sim.Y, nrow=R*S, ncol=nlag)
      max.Y  <- apply(sim.Y, 1, max)
      max.Y  <- matrix(max.Y, nrow=R, ncol=S)
      th.tmp <- colSums(max.Y < mesh.L[i])/R
      distr[i,] <- th.tmp
      th[i,]    <- c(mean(th.tmp), median(th.tmp), quantile(th.tmp, levels))
    }## Monte Carlo integration
    else if(grepl("MCi",method[1])){
      mesh.Z <- (mesh.L[i]-sim.L%*%t(a))/exp(log(sim.L)%*%t(b))# [RxS(m-1)]
      mesh.Z <- array(mesh.Z, dim=c(R,S,nlag))
      HZ      <- vapply(1:S, p.res, z=mesh.Z, mu=m, sig=s, w=w, FUN.VALUE=numeric(R))# [RxS]
      th.samp <- colMeans(HZ)# [S]
      distr[i,] <- t(th.samp)
      th[i,]    <- cbind(colMeans(th.samp), apply(th.samp, 2, median),
                             t(apply(th.samp, 2, quantile, levels)))
    }
  }
  # return and summaries
  ret$theta <- th
  ret$distr <- distr
  time <- toc(silent=TRUE)
  if(!silent){
    print(paste("Time elapsed on estimating theta: ",
                time%/%60," min ",round(time-60*(time%/%60), 1)," sec.", sep=""))
    print(paste("Time per x-value: ",round(time/nbr.vert, 1)," sec.", sep=""))
  }
  return(ret)
}



##################################################
## CHI(X) ESTIMATION

## >   ts: vector of reals, time series which to estimate theta(x,m) for
## >   lapl: boolean, is [ts] already with Laplace margins?
## >   m: integer>2, maximum lag for which to compute chi
## >   R: integer>0, nbr of samples used in a single sweep for the estimation
## >   S: integer>0, nbr of posterior sweeps used for the estimation
## >   u.mar: probability, threshold used for marginal threshold
## >   u.dep: probability, threshold used for Heffernan-Tawn model
## >   probs: vector of probabilities, x in theta(x,m)
## >   method.mar: string, "mle" for max likelihood or "mom" for the method of moments or "pwm" for proba weighted moments
## >   method: (vector of) string(s), either "prop" or "MCi"
## >   silent: boolean, verbosity
## >   fit: boolean or list returned from htfit, TRUE means the htfit must be called
## >   par: list, contains all arguments needed for the call to htfit - if empty, default values are assumed
## >   submodel: string, structure imposed on (a,b) - defaults to "fom"
## >   levels: vector of probabilities, specifies which posterior quantiles to compute (+mean+median)
## <   ret: list, ht fit, chi(.MC) [(nbr of mesh nodes) x (nbr of quantiles)]
## .   called by user
chifit <- function(ts, lapl=FALSE, nlag=1, R=1000, S=500,
                   u.mar=0, u.dep, probs=seq(u.dep,0.9999,length.out=30),
                   method.mar=c("mle","mom","pwm"), method=c("prop","MCi"),
                   silent=FALSE,
                   fit=TRUE, prev.fit=bayesfit(), par=bayesparams(),
                   submodel=c("fom","none"), levels=c(.025,.975)){
  # pre-process data
  data.up <- format.ts(ts=ts, u.mar=u.mar, u.dep=u.dep, method=method.mar,
                       nlag=nlag, lapl=lapl)
  n       <- dim(data.up)[1]
  mesh    <- probs
  mesh.O  <- scale.to.original(p=mesh, ts=ts, u=u.mar, gpd.pars=gpd(ts, u=u.mar, method=method.mar)$pars)
  # H+T fit
  if(fit){
    if(!is.bayesparams(par)) stop("par must be of class 'bayesparams'.")
    fit <- htfit(data.up, prop.a=par$prop.a, prop.b=par$prop.b,
                prior.mu=par$prior.mu, prior.nu=par$prior.nu, prior.eta=par$prior.eta,
                trunc=par$trunc, comp.saved=par$comp.saved, maxit=par$maxit,
                burn=par$burn, thin=par$thin,
                adapt=par$adapt, batch.size=par$batch.size, mode=par$mode, submodel=submodel)
    if(!silent)
      print(paste("Fit of DDP done. Mean value for alpha & beta:",
                  signif(apply(fit$a, 2, mean),3),"and",
                  signif(apply(fit$b, 2, mean),3)))
  }else{
    if(!is.bayesfit(prev.fit)) stop("prev.fit must be of class 'bayesfit'.")
    if(prev.fit$len == 0) stop("prev.fit must contain a proper instance of class 'bayesfit'.")
    fit <- prev.fit
  }
  tic()
  ret <- depmeasure("chi")
  ret$bayesfit <- fit
  ret$probs    <- mesh
  ret$levels   <- mesh.O
  ret$nlag     <- nlag
  # call etfit on (X_1,X_j) to get distr and/or distr.MC
  len       <- dim(fit$a)[1]
  n.comp    <- dim(fit$mu)[2]
  nbr.vert  <- length(mesh)
  nbr.quant <- length(levels)+2
  chi <- matrix(0, nrow=nbr.vert, ncol=nbr.quant)
  dimnames(chi) <- list(NULL,c("mean","median",paste(levels*100,"%",sep="")))
  subfit <- bayesfit()
  subfit$a <- fit$a[,nlag, drop=FALSE]
  subfit$b <- fit$b[,nlag, drop=FALSE]
  subfit$sd   <- fit$sd[,,nlag, drop=FALSE]
  subfit$mean <- fit$mean[,,nlag, drop=FALSE]
  subfit$w   <- fit$w
  subfit$len <- fit$len
  fit.sample <- etfit(data=data.up[,c(1,nlag+1)], R=R, S=S, probs=mesh, method=method,
                      silent=silent,
                      fit=FALSE, prev.fit=subfit, submodel=submodel, levels=levels)
  # compute chi_j(x) = Pr(X_j>x | X_1>x)
  chi[,"mean"]   <- 1-fit.sample$theta[,"mean"]
  chi[,"median"] <- 1-fit.sample$theta[,"median"]
  if(nbr.quant>2){
    chi[,3:nbr.quant] <- t(vapply(seq_along(probs), function(i){
      quantile(1-fit.sample$distr[i,], probs=levels)
    }, levels))
  }
  distr <- 1-fit.sample$distr
  ret$chi   <- chi
  ret$distr <- distr
  ## print time
  time <- toc(silent=TRUE)
  if(!silent){
    print(paste("Time elapsed on estimating chi: ",
                time%/%60," min ",round(time-60*(time%/%60), 1)," sec.", sep=""))
    print(paste("Time per x-value: ",round(time/nbr.vert, 1)," sec.", sep=""))
  }
  return(ret)
}