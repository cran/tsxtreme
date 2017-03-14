## Copyright (C) 2017 Thomas Lugrin
## user inteface for htfit class
## (Scope:) Fit H+T model using Bayesian semiparametrics
## List of functions: - depfit
##                    - htfit
##################################################

##################################################
## H+T R interface

depfit <- function(ts, u.mar=0, u.dep=u.mar, lapl=FALSE, method.mar=c("mle","mom","pwm"), nlag=1,
                   par=bayesparams(), submodel=c("fom","none","ugm")){
  if(submodel[1]=="ugm"){
    data.up <- cbind(1,ts)
  }else if(submodel[1] %in% c("fom","none")){
    data.up <- format.ts(ts=ts, u.mar=u.mar, u.dep=u.dep, method=method.mar,
                         nlag=nlag, lapl=lapl)
  }else{
    stop(paste("submodel",submodel,"unknown in depfit()"))
  }
  htfit(data=data.up, prop.a=par$prop.a, prop.b=par$prop.b,
        prior.mu=par$prior.mu, prior.nu=par$prior.nu, prior.eta=par$prior.eta,
        trunc=par$trunc, comp.saved=par$comp.saved,
        maxit=par$maxit, burn=par$burn, thin=par$thin,
        adapt=par$adapt, batch.size=par$batch.size,
        mode=par$mode, submodel=submodel)
}


## >   data: matrix, first column is > threshold -- other columns are lagged versions of the first one
## >   prop.a,prop.b: scalar or vector of length 5 (resp. 3), standard deviation of Gaussian proposal
## >   prior.mu, prior.nu, prior.eta: vectors of length 2, normal - inverse gamma - gamma (shape,scale) parameters
## >   trunc: integer, DP approximation of infinite sum
## >   comp.saved: integer, characteristics only of the first [comp.saved] components are saved in output
## >   maxit: integer, number of sweeps in MCMC
## >   burn: integer, burn-in before saving traces
## >   thin: integer, thinning (1 means keeping all values)
## >   adapt: integer, used for RAMA on (a,b): 0 means no adaption, otherwise specifies when adaption stops
## >   batch.size: integer, size of batches for RAMA, between 100 and 1000 typically
## >   mode: {0,1,2}, console output during computations, 0: debug mode - 1: normal mode - 2: silent mode
## >   submodel: string, structure imposed on (a,b) - defaults to "fom" (first order Markov)
## <   ret: list, MCMC traces
## .   called by user
## ... wrap function for etfit_externC

htfit <- function(data,
                  prop.a, prop.b,
                  prior.mu=c(0,10), prior.nu=c(2,1/2), prior.eta=c(4,1),
                  trunc=100, comp.saved=10,
                  maxit=10000, burn=2000, thin=4,
                  adapt=2000, batch.size=125,
                  mode=1, submodel=c("fom","none","ugm")){# "ugm": univariate Gaussian mixture
  n    <- dim(data)[1]
  nlag <- dim(data)[2]-1
  if(length(prop.a) < 5){ prop.a <- rep(prop.a[1],5) }
  if(length(prop.b) < 3){ prop.b <- rep(prop.b[1],3) }
  ## build trace containers
  tr.len <- floor((maxit-burn)/thin)
  t.a   <- t.b <- matrix(0, nrow=tr.len, ncol=nlag)
  t.sig <- t.mu <- array(0, dim = c(tr.len,comp.saved,nlag))
  t.w   <- matrix(0, nrow=tr.len, ncol=comp.saved)
  t.g   <- numeric(tr.len)
  t.ci  <- matrix(0, nrow=tr.len, ncol=n)
  t.noo <- matrix(0, nrow=tr.len, ncol=comp.saved)
  t.noc <- numeric(tr.len)
  t.sd  <- array(0, dim=c(tr.len,8,nlag))#for RAMA
  if(submodel[1] == "ugm"){ submodel <- 0 }
  else if(submodel[1]=="fom"){ submodel <- 1 }
  else if(submodel[1]=="none"){ submodel <- 2 }
  else{ stop("In htfit(): invalid submodel.") }
  ## call C(++) wrapper
  fit = .C(C_et_interface, as.double(data), as.integer(n), as.integer(nlag),
           as.integer(trunc), as.integer(comp.saved), as.integer(maxit),
           as.integer(burn), as.integer(thin), as.integer(adapt), as.integer(batch.size),
           as.double(prop.a), as.double(prop.b),
           as.double(prior.mu), as.double(prior.nu), as.double(prior.eta),
           as.integer(mode), as.integer(submodel),
           as.double(t.a), as.double(t.b), as.double(t.sig),
           as.double(t.mu), as.double(t.w), as.double(t.g),
           as.integer(t.ci), as.integer(t.noo), as.integer(t.noc), as.double(t.sd))
  ## format C output
  ret <- bayesfit()
  ret$a    <- matrix(fit[[18]], nrow=tr.len, ncol=nlag)
  ret$b    <- matrix(fit[[19]], nrow=tr.len, ncol=nlag)
  ret$sd   <- array(fit[[20]], dim = c(tr.len,comp.saved,nlag))
  ret$mean <- array(fit[[21]], dim = c(tr.len,comp.saved,nlag))
  ret$w       <- matrix(fit[[22]], nrow=tr.len, ncol=comp.saved)
  ret$prec    <- fit[[23]]
  ret$ci      <- matrix(fit[[24]], nrow=tr.len, ncol=n)
  ret$noo     <- matrix(fit[[25]], nrow=tr.len, ncol=comp.saved)
  ret$noc     <- fit[[26]]
  ret$prop.sd <- array(fit[[27]], dim=c(tr.len, 8, nlag))
  ret$len <- tr.len
  ret$nlag <- nlag
  ## add names
  ind <- 1:nlag
  colnames(t.a) <- paste("alpha(",ind,")", sep="")
  colnames(t.b) <- paste("beta(",ind,")", sep="")
  #return in R format
  return(ret)
}