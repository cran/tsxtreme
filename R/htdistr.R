## Copyright (C) 2017 Thomas Lugrin
## Functions related to H+T model
## (Scope:) Heffernan-Tawn model fit in 1 stage
## List of functions: - p.res
##                    - r.res
##################################################


## >   s: integer, ranges across sweeps
## >   z: array of reals, residuals
## >   mu,sig: array of reals, means and sd of components [SxKx(m-1)]
## >   w: matrix of reals, weights of components - rows sum to 1
## <   ret: matrix of reals, evaluation of H+T residual distribution function
## .   called by [etfit()] for MC integration
p.res <- function(s, z, mu, sig, w){
  R   <- dim(z)[1]; nlag = dim(z)[3]
  ret <- numeric(R)
  for(r in 1:R){
    prev <- 1
    for(lag in 1:nlag){
      prev <- prev*sum(w[s,]*pnorm(z[r,s,lag], mu[s,,lag], sig[s,,lag]))
    }
    ret[r] <- prev
  }
  return(ret)
}


## >   R: integer>0, nbr of samples used in a single sweep for the estimation
## >   S: integer>0, nbr of posterior sweeps used for the estimation
## >   nlag: integer>0, nbr of lags
## >   w,m,s: matrix - arrays, output from depfit
## <   sim.Z: matrix of reals, samples of residuals across sweeps
## .   called by [etfit()] for proportion method
r.res <- function(R,S,nlag,w,m,s){
  sim.Z <- array(0, dim=c(R,S,nlag))
  
  for(sw in 1:S){# loop on selected sweeps
    hwmany <- rmultinom(1, R, w[sw,])
    r      <- 1
    c      <- 1
    
    while(r<R){# sample from mixture
      if(hwmany[c]>0){
        if(dim(s)[3]==1){
          sim.Z[r:(r+hwmany[c]-1),sw,] <- rnorm(hwmany[c], m[sw,c,], s[sw,c,])
        }else{
          sim.Z[r:(r+hwmany[c]-1),sw,] <- rmvnorm(hwmany[c], m[sw,c,], diag(s[sw,c,]^2))
        }
      }
      r <- r+hwmany[c]
      c <- c+1
    }
  }
  sim.Z <- matrix(sim.Z, nrow=R, ncol=S*nlag)
  return(sim.Z)
}