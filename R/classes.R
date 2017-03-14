## Copyright (C) 2017 Thomas Lugrin
## Definition of methods and related functions
## (Scope:) single framework for 1- and 2-stage functions
## (Note:) summary is tricky especially when multi-dimensional!
## List of functions: - stepfit
##                    - is.stepfit
##                    - print/summary.stepfit
##                    - plot.stepfit
##                    - bayesfit
##                    - is.bayesfit
##                    - print/summary.bayesfit
##                    - plot.bayesfit
##                    - bayesparams
##                    - is.bayesparams
##                    - print/summary.bayesparams
##                    - depmeasure
##                    - is.depmeasure
##                    - print/summary.depmeasure
##                    - plot.depmeasure
##                    - is.int
##################################################

##################################################
## CLASS "STEPFIT"
## constructor
stepfit <- function(){
  x <- list(a=numeric(0), b=numeric(0), res=matrix(0,nrow=1,ncol=1), pars.se=matrix(0,nrow=1,ncol=2), nlag=0)
  class(x) <- "stepfit"
  return(x)
}

## validator
is.stepfit <- function(x){
  if(length(names(x))){
    conds <- names(stepfit()) %in% names(x)# minimal requirement
    conds <- sum(conds) == length(names(stepfit()))
    if(conds){
      conds <- is.numeric(x$a) && is.numeric(x$b) && is.matrix(x$res) && is.matrix(x$pars.se)
      if(!conds) return(FALSE)
      conds <- dim(x$pars.se)[1]==2 && length(x$a)==length(x$b) && dim(x$pars.se)[2]==length(x$a) && length(x$nlag) == 1
      return(inherits(x, "stepfit") && conds)
    }
  }
  return(FALSE)
}

## viewer
print.stepfit <- function(x, ...){
  summary.stepfit(x, ...)
}

summary.stepfit <- function(object, ...){
  cat(" Parameters:\n")
  for(ind in 1:object$nlag){
    cat(paste("alpha(",ind,")  ",signif(object$a,3)," (",signif(object$pars.se[1,ind],3),")\n", sep=""), sep="")
    cat(paste("beta(",ind,")   ",signif(object$b,3)," (",signif(object$pars.se[2,ind],3),")\n", sep=""), sep="")
  }
  cat(" Residuals:\n")
  print(summary(object$res))
}

## plotter
plot.stepfit <- function(x, ...){
  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))
  par(ask=TRUE)
  def <- list(xlab="z", ylab=bquote(h(z)), main="")
  if(length(list(...)) > 0){
    pm <- !pmatch(names(def), names(list(...)), 0)# all non-matches
    def <- def[pm]
  }
  for(j in 1:x$nlag){
    if("ylab" %in% names(def)) def$ylab <- bquote(bquote(h(z[.(j)])))
    if("xlab" %in% names(def)) def$xlab <- bquote(bquote(z[.(j)]))
    den <- density(x$res[,j])
    def$x <- x$res[,j]
    do.call(hist, c(def, prob=TRUE, list(...)))
    lines(den, lty="dashed")
  }
}


##################################################
## CLASS "BAYESFIT"
## constructor
bayesfit <- function(){
  z.mat <- matrix(0, nrow=1, ncol=1)
  z.ar  <- array(0, dim=c(1,1,1))
  x <- list(a=z.mat, b=z.mat, sd=z.ar, mean=z.ar, w=z.mat, prec=0,
              ci=z.mat, noo=z.mat, noc=0, prop.sd=z.ar,
              len=0, nlag=0)
  class(x) <- "bayesfit"
  return(x)
}

## validator
is.bayesfit <- function(x){
  if(length(names(x))){
    conds <- names(bayesfit()) %in% names(x)
    conds <- sum(conds) == length(names(bayesfit()))
    if(all(conds)){
      conds <- is.matrix(x$a) && is.matrix(x$b) &&
        is.array(x$sd) && is.array(x$mean) &&
        is.matrix(x$w) && is.numeric(x$prec) && is.array(x$prop.sd) &&
        is.matrix(x$ci) && is.matrix(x$noo) && is.numeric(x$noc) && is.numeric(x$len) && is.numeric(x$nlag) &&
        is.int(x$len) && length(x$len)==1 && is.int(x$nlag) && length(x$nlag)==1 &&
        is.int(x$noc) && is.int(x$noo) && is.int(x$ci)
      return(inherits(x, "bayesfit") && conds)
    }
  }
  return(FALSE)
}

## viewer
print.bayesfit <- function(x, ...){
  summary.bayesfit(x, ...)
}

summary.bayesfit <- function(object, ...){
  cat(" Posterior medians:\n")
  for(j in 1:object$nlag){
    cat("alpha(",j,")  ",median(object$a[,j]),"\n", sep="")
    cat("beta(",j,")   ",median(object$b[,j]),"\n", sep="")
    cat("mean(",j,",",1,") ",median(object$mean[,1,j]),"\n", sep="")
    cat("sd(",j,",",1,")   ",median(object$sd[,1,j]),"\n", sep="")
  }
}

plot.bayesfit <- function(x, which=1:3, ...){
  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))
  stopifnot(is.bayesfit(x))
  stopifnot(is.int(which))
  stopifnot(min(which) >= 1 && max(which) <= 3)
  which <- sort(unique(which))
  ugm <- all(x$a==0) && all(x$b==0)
  if(ugm) which <- intersect(which,1:2)
  if(x$nlag > 1) par(ask=TRUE)
  if(length(which) == 2) par(mfrow=c(1,2))
  else if(length(which) > 2) par(mfrow=c(2,2))
  def <- list(xlab="x", ylab="y", main="")
  if(length(list(...)) > 0){
    pm <- !pmatch(names(def), names(list(...)), 0)# all non-matches
    def <- def[pm]
  }
  for(j in 1:x$nlag){
    if(1 %in% which){
      if("xlab" %in% names(def)) def$xlab <- bquote(bquote(z[.(j)]))
      if("ylab" %in% names(def)) def$ylab <- bquote(bquote(h(z[.(j)])))
      if("main" %in% names(def)) def$main <- "Sample of residual densities"
        grid <- seq(-100,100,2)
        de <- residual.densities(x,j,grid)
        grid <- seq(min(de$x),max(de$x),length.out=151)
        de <- residual.densities(x,j,grid)
        xy <- list(x=de$x, y=de$y[,1])
        if(!("ylim" %in% names(list(...)))) xy$ylim <- c(0,max(de$y))
        do.call(plot, c(xy, type="l",def,list(...)))
        for(it in 1:dim(de$y)[2])
          lines(de$x, de$y[,it])
    }
    if(2 %in% which){
      if("xlab" %in% names(def)) def$xlab <- bquote(bquote(z[.(j)]))
      if("ylab" %in% names(def)) def$ylab <- bquote(bquote(H(z[.(j)])))
      if("main" %in% names(def)) def$main <- "Residual distribution"
      grid <- seq(-100,100,2)
      ds <- residual.distributions(x,j,grid)
      grid <- seq(min(ds$x),max(ds$x),length.out=151)
      ds <- residual.distributions(x,j,grid)
      xy <- list(x=ds$x, y=ds$y[,1])
      do.call(plot, c(xy, type="l",def,list(...)))
      for(q in 2:3)
        lines(ds$x, ds$y[,q], lty="dashed")
    }
    if(3 %in% which){
      if("xlab" %in% names(def)) def$xlab <- bquote(bquote(alpha[.(j)]))
      if("ylab" %in% names(def)) def$ylab <- bquote(bquote(beta[.(j)]))
      if("main" %in% names(def)) def$main <- "Joint posterior of H-T par."
      lims <- c(max(min(x$a[,j]),-1),min(max(x$a[,j]),1),max(min(x$b[,j]),0),min(max(x$b[,j]),1))
      kd <- kde2d(x$a[,j], x$b[,j], lims=lims)
      do.call(contour, c(kd,nlevels=5,def,list(...)))
    }
  }
}

residual.densities <- function(x, lag, grid){
  nsamp<- min(25,x$len)# nbr of samples
  noc  <- dim(x$mean)[2]
  nobs <- sum(x$noo[1,])
  ind <- sample.int(x$len, nsamp, replace=FALSE)
  dens <- vapply(ind, function(it,j){
    de <- vapply(1:noc, function(c,j,it){
        x$noo[it,c]*dnorm(grid, x$mean[it,c,j], x$sd[it,c,j])
      }, grid, j=j, it=it)
      rowSums(de)# sample density
    }, grid, j=lag)
  dens <- dens/nobs
  maxs <- vapply(1:nsamp, function(j) max(dens[,j]), 0)
  wmin <- which.min(maxs)# expected to be the flattest
  minimax <- min(maxs)
  fact <- 1e3# trim tails which are at least [fact] times smaller than minimax
  lo <- 1; up <- length(grid)
  i <- 1
  while(i < up && dens[i,wmin] < minimax/fact)
    i <- i+1
  lo <- max(lo,i-3); i <- up
  while(i > lo && dens[i,wmin] < minimax/fact)
    i <- i-1
  up <- min(up,i+3)
  return(list(x=grid[lo:up], y=dens[lo:up,]))
}

residual.distributions <- function(x, lag, grid){
  nsamp<- min(500,x$len)# nbr of samples
  noc  <- dim(x$mean)[2]
  nobs <- sum(x$noo[1,])
  ind <- sample.int(x$len, nsamp, replace=FALSE)
  dstr <- vapply(ind, function(it, j){
    ds <- vapply(1:noc, function(c,j,it){
      x$noo[it,c]*pnorm(grid, x$mean[it,c,j], x$sd[it,c,j])
    }, grid, j=j, it=it)
    rowSums(ds)# sample distribution
  }, grid, j=lag)
  dstr <- dstr/nobs
  dstr <- vapply(seq_along(grid), function(gr){
    quantile(dstr[gr,], probs=c(.5,.05,.95))
  }, numeric(3))
  dstr <- t(dstr)
  # trim tails of distribution
  lo <- 1; up <- length(grid)
  if(grid[2]-grid[1] > 1){
    eps.prob <- 1e-3
    i <- 1
    while(i < up && dstr[i,1] < eps.prob)
      i <- i+1
    lo <- max(lo,i-1); i <- up
    while(i > lo && dstr[i,1] > 1-eps.prob)
      i <- i-1
    up <- min(up,i+1)
  }
  return(list(x=grid[lo:up], y=dstr[lo:up,]))
}


##################################################
## CLASS "BAYESPARAMS"
## constructor
bayesparams <- function(prop.a=0.02, prop.b=0.02,
                        prior.mu=c(0,10), prior.nu=c(2,1/2), prior.eta=c(2,2),
                        trunc=100, comp.saved=15,
                        maxit=30000, burn=5000, thin=1,
                        adapt=5000, batch.size=125, mode=1){
  x <- list(prop.a=prop.a, prop.b=prop.b,
              prior.mu=prior.mu, prior.nu=prior.nu, prior.eta=prior.eta,
              trunc=trunc, comp.saved=comp.saved,
              maxit=maxit, burn=burn, thin=thin,
              adapt=adapt, batch.size=batch.size, mode=mode)
  class(x) <- "bayesparams"
  if(is.bayesparams(x)) return(x)
  else stop("Wrong type of argument. See help for details.")
}

## validator
is.bayesparams <- function(x){
  if(length(names(x))){
    names.bp <- c("prop.a","prop.b","prior.mu","prior.nu","prior.eta",
                  "trunc","comp.saved","maxit","burn","thin",
                  "adapt","batch.size","mode")
    conds <- names.bp %in% names(x)
    if(all(conds)){
      conds <- is.numeric(x$prop.a) && is.numeric(x$prop.b) &&
        length(x$prop.a)==1 && length(x$prop.b)==1 && x$prop.a>0 && x$prop.b>0 &&
        is.numeric(x$prior.mu) && is.numeric(x$prior.nu) && is.numeric(x$prior.eta) &&
        length(x$prior.mu)==2 && length(x$prior.nu)==2 && length(x$prior.eta)==2 &&
        x$prior.mu[2] > 0 && all(x$prior.nu>0) && all(x$prior.eta>0) &&
        length(x$trunc)==1 && is.int(x$trunc) && length(x$comp.saved)==1 && is.int(x$comp.saved) &&
        length(x$maxit)==1 && is.int(x$maxit) && length(x$burn)==1 && is.int(x$burn) &&
        length(x$thin)==1 && is.int(x$burn) && length(x$adapt)==1 && is.int(x$adapt) &&
        length(x$batch.size)==1 && is.int(x$batch.size) && length(mode)==1 && x$mode %in% c(0,1,2)
      return(inherits(x, "bayesparams") && conds)
    }     
  }
  return(FALSE)
}

## viewer
print.bayesparams <- function(x, ...){
  summary.bayesparams(x, ...)
}

summary.bayesparams <- function(object, ...){
  if(!is.bayesparams(object)){
    summary.default(object)
  }else{
    cat("proposal for alpha            ",object$prop.a,"\n", sep="")
    cat("proposal for beta             ",object$prop.b,"\n", sep="")
    cat("prior par. for means          (",object$prior.mu[1],", ",object$prior.mu[2],")\n", sep="")
    cat("prior par. for sd             (",object$prior.nu[1],", ",object$prior.nu[2],")\n", sep="")
    cat("prior par. for precision par. (",object$prior.eta[1],", ",object$prior.eta[2],")\n", sep="")
    cat("DP truncation                 ",object$trunc,"\n", sep="")
    cat("nbr of components saved       ",object$comp.saved,"\n", sep="")
    cat("max. nbr of MCMC iterations   ",object$maxit,"\n", sep="")
    cat("length of burn-in             ",object$burn,"\n", sep="")
    cat("thinning par.                 ",object$thin,"\n", sep="")
    cat("length of adaption            ",object$adapt,"\n", sep="")
    cat("size of adaption batches      ",object$batch.size,"\n", sep="")
    cat("debug mode                    ",object$mode,"\n", sep="")
  }
  invisible()
}


##################################################
## CLASS "DEPMEASURE"
## constructor
depmeasure <- function(type=c("theta","chi","steptheta","runs")){
  z.mat <- matrix(0, nrow=1, ncol=1)
  x <- list(probs=0, levels=0, nlag=1)
  if(type[1]=="theta" || type[1]=="chi"){
    x$fit   <- bayesfit()
    x$distr <- z.mat
    if(type[1]=="theta")
      x$theta <- z.mat
    else if(type[1]=="chi")
      x$chi <- z.mat
  }else if(type[1]=="steptheta"){
    x$fit <- stepfit()
    x$theta <- z.mat
  }else if(type[1]=="runs"){
    x$theta <- z.mat
    x$nbr.exc <- 0
  }else{
    stop("type not recognised. Must be either 'theta' or 'chi'.")
  }
  class(x) <- "depmeasure"
  return(x)
}

## validator
is.depmeasure <- function(x){
  if(length(names(x))){
    conds <- names(depmeasure("theta")) %in% names(x)
    if(all(conds)){
      conds <- is.bayesfit(x$fit) && is.matrix(x$theta)
      return(inherits(x, "depmeasure") && conds)
    }
    conds <- names(depmeasure("chi")) %in% names(x)
    if(all(conds)){
      conds <- is.bayesfit(x$fit) && is.matrix(x$chi)
      return(inherits(x, "depmeasure") && conds) 
    }
    conds <- names(depmeasure("steptheta")) %in% names(x)
    if(all(conds)){
      conds <- is.stepfit(x$fit) && is.matrix(x$theta)
      return(inherits(x, "depmeasure") && conds) 
    }
    conds <- names(depmeasure("runs")) %in% names(x)
    if(all(conds)){
      conds <- is.matrix(x$theta)
      return(inherits(x, "depmeasure") && conds)
    }
  }
  return(FALSE)
}

## viewer
print.depmeasure <- function(x, ...){
  summary.depmeasure(x, ...)
}

summary.depmeasure <- function(object, ...){
  cat("Bayes fit\n")
  summary(object$fit)
  cat(paste(names(object)[2],"posterior estimates\n"))
  print(summary(object[[2]]))
  invisible()
}

## plotter
plot.depmeasure <- function(x, ...){
  stopifnot(is.depmeasure(x))
  if("nbr.exc" %in% names(x)){
    ord <- order(colnames(x$theta))
    def <- list(xlab="x", ylab=bquote(bquote(theta(x,.(x$nlag)))), ylim=range(x$theta))
    if(length(list(...)) > 0){
      pm <- !pmatch(names(def), names(list(...)), 0)# all non-matches
      def <- def[pm]
    }
    def$x <- x$levels; def$y <- x$theta[,"estimate"]
    do.call(plot, c(type="l", def, list(...)))
    ord <- ord[-length(ord)]
    # add credible/confidence lines
    if(length(ord) > 1){
      for(i in 0:(floor(length(ord)/2)-1)){
        lines(x$levels, x$theta[,ord[i*2+1]], lty=i+2)
        lines(x$levels, x$theta[,ord[i*2+2]], lty=i+2)
      }
    }
  }else if("chi" %in% names(x)){
    ord <- order(colnames(x$chi))
    def <- list(xlab="x", ylab=bquote(bquote(chi[.(x$nlag)](x))), ylim=range(x$chi))
    if(length(list(...)) > 0){
      pm <- !pmatch(names(def), names(list(...)), 0)
      def <- def[pm]
    }
    def$x <- x$levels; def$y <- x$chi[,"median"]
    do.call(plot, c(def,type="l",list(...)))
    ord <- ord[-c(length(ord)-1,length(ord))]
    # add credible/confidence lines
    if(length(ord) > 1){
      for(i in 0:(floor(length(ord)/2)-1)){
        lines(x$levels, x$chi[,ord[i*2+1]], lty=i+2)
        lines(x$levels, x$chi[,ord[i*2+2]], lty=i+2)
      }
    }
  }else{
    ord <- order(colnames(x$theta))
    def <- list(xlab="x", ylab=bquote(bquote(theta(x,.(x$nlag)))), ylim=range(x$theta))
    if(length(list(...)) > 0){
      pm <- !pmatch(names(def), names(list(...)), 0)
      def <- def[pm]
    }
    def$x <- x$levels; ifelse(is.stepfit(x$fit), def$y <- x$theta[,"estimate"], def$y <- x$theta[,"median"])
    do.call(plot, c(def,type="l",list(...)))
    ord <- ord[-length(ord)]
    if(is.bayesfit(x$fit)) ord <- ord[-length(ord)]
    # add credible/confidence lines
    if(length(ord) > 1){
      for(i in 0:(floor(length(ord)/2)-1)){
        lines(x$levels, x$theta[,ord[i*2+1]], lty=i+2)
        lines(x$levels, x$theta[,ord[i*2+2]], lty=i+2)
      }
    }
  }
}


##################################################
## ADDITIONAL VALIDATORS
is.int <- function(x){
  if(all(x%/%1==x)) return(TRUE)
  else return(FALSE)
}