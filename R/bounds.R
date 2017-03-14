## Copyright (C) 2017 Thomas Lugrin
## Compute Keef, Papastathopoulos & Tawn bounds on (a,b)
## (Scope:) Heffernan-Tawn model fit in 2 stages
## List of functions: - conditions.verify
##                    - conditions.bounds
##                    - dicho.b.bounds
##                    - dicho.a.bounds
##################################################

##################################################
## Bounds on (alpha,beta)

## >   a: scalar, alpha parameter, in [-1,1]
## >   b: scalar, beta parameter, in [0,1]
## >   p: scalar, probability (to compute quantiles of residual distribution), in [0,1]
## >   v: scalar, threshold, recommended as beyond the maximum observation
## >   data: bivariate vector, (X,Y) with Y | X>u
## <   boolean, does the couple (a,b) satisfy the conditions?
## .   called by user to get the bounds on the posterior samples of (a,b) and by [conditions.bounds()]
conditions.verify <- function(a, b, p, v=-log(2*(1-0.99999999)), data){
  z.q.pos <- q.res2(p, 1, 0, data)
  z.q     <- q.res2(p, a, b, data)
  z.q.neg <- q.res2(p, -1, 0, data)
  if(a < -1 || a > 1 || b >= 1 || b<0){ return(FALSE) }
  
  # Case I
  case.1 <- FALSE
  if(a <= min(1,  1 - b*z.q*v^{b-1},  1 - v^{b-1}*z.q + z.q.pos/v)){ case.1 <- TRUE }
  else{
    cond.1 <- (1 - b*z.q*v^{b-1} < a)
    if(b*z.q>0) cond.2 <- ( (1-1/b) * exp( log(b*z.q)*(1/(1-b)) + log(1-a)*(-b/(1-b)) ) + z.q.pos > 0 )
    else        cond.2 <- ( (1-1/b) * (b*z.q)^{1/(1-b)} * (1-a)^{-b/(1-b)} + z.q.pos > 0 )
    if(cond.1 && cond.2){ case.1 <- TRUE }
    else{ return(FALSE) }
  }
  # Case II
  case.2 <- FALSE
  if(-a <= min(1 + b*z.q*v^{b-1}, 1 + v^{b-1}*z.q - z.q.neg/v)){ case.2 <- TRUE; return(TRUE) }
  else{
    cond.1 <- ( 1 + b*v^{b-1}*z.q < -a )
    if(-b*z.q>0) cond.2 <- ( (1-1/b) * exp( log(-b*z.q)*(1/(1-b)) + log(1+a)*(-b/(1-b)) ) - z.q.neg > 0 )
    else         cond.2 <- ( (1-1/b) * (-b*z.q)^{1/(1-b)} * (1+a)^{-b/(1-b)} - z.q.neg > 0 )
    if(cond.1 && cond.2){ case.2 <- TRUE; return(TRUE) }
    else{ return(FALSE) }
  }
}

## >   fix.par: scalar, either value of alpha or of beta
## >   fix.alpha: boolean, TRUE if fix.par is for alpha, FALSE if fix.par is for beta
## >   data: bivariate vector, (X,Y) with Y | X>u
## >   eps: scalar, tolerance, typically 0.001
## <   vector of length 2, with the bounds found
## .   called by user to get the bounds and thus compute the normalising constant for truncated normal proposal
conditions.bounds <- function(fix.par, fix.alpha=TRUE, data, eps=0.001, v=-log(2*(1-0.99999999))){
  i <- o <- 2
  bounds <- numeric(2)
  if(fix.alpha){
    if(conditions.verify(fix.par, 0, 1, v, data) &&
         conditions.verify(fix.par, 0, 0, v, data)){ i <- 0 }
    else{ o <- 0 }
    if(conditions.verify(fix.par, 1, 1, v, data) &&
         conditions.verify(fix.par, 1, 0, v, data)){ i <- 1 }
    else{ o <- 1 }
    if(i == 2){
      grid <- seq(0, 1, 0.05)
      pt   <- 1
      while(i == 2 && pt <= length(grid)){
        if(conditions.verify(fix.par, grid[pt], 1, v, data) &&
             conditions.verify(fix.par, grid[pt], 0, v, data)){ i <- grid[pt] }
        pt <- pt+1
      }
      if(i == 2){ stop("Refine grid on beta to get proper results... in [conditions.bounds].") }
    }else if(o == 2){
      bounds[1] <- -1
      bounds[2] <- 1
      return(bounds)
    }
    bounds[1] <- i
    bounds[2] <- dicho.b.bounds(i, o, fix.par, data, eps, v)
  }else{
    if(conditions.verify(-1, fix.par, 1, v, data) &&
         conditions.verify(-1, fix.par, 0, v, data)){ i <- -1 }
    else{ o <- -1 }
    if(conditions.verify(1, fix.par, 1, v, data) &&
         conditions.verify(1, fix.par, 0, v, data)){ i <- 1 }
    else{ o <- 1 }
    if(i == 2){
      grid <- seq(-1, 1, 0.05)
      pt   <- 1
      while(i == 2 && pt <= length(grid)){
        if(conditions.verify(grid[pt], fix.par, 1, v, data) &&
             conditions.verify(grid[pt], fix.par, 0, v, data)){ i <- grid[pt] }
        pt <- pt+1
      }
      if(i == 2){ stop("Refine grid on alpha to get proper results... in [conditions.bounds].") }
    }else if(o == 2){
      bounds[1] <- -1
      bounds[2] <- 1
      return(bounds)
    }
    bounds[1] <- i
    bounds[2] <- dicho.a.bounds(i, o, fix.par, data, eps, v)
  }
  if(bounds[1] > bounds[2]){
    v         <- bounds[2]
    bounds[2] <- bounds[1]
    bounds[1] <- v
  }
  return(bounds)
}

## >   i: current inner bound
## >   o: current outer bound
## >   fix.par: scalar, either value of alpha or of beta
## >   data: bivariate vector, (X,Y) with Y | X>u
## >   eps: scalar, tolerance, typically 0.001
## <   scalar, the bound found
## .   pair of functions called by conditions.bounds to get more accurate bounds by dichotomy
dicho.b.bounds <- function(i, o, fix.par, data, eps, v=-log(2*(1-0.99999999))){
  mid <- mean(c(i,o))
  if(conditions.verify(fix.par, mid, 1, v, data) &&
       conditions.verify(fix.par, mid, 0, v, data)){ i<- mid }
  else{ o <- mid }
  if(abs(i-o) < eps){ return(i) }
  else{ return(dicho.b.bounds(i, o, fix.par, data, eps, v)) }
}

dicho.a.bounds <- function(i, o, fix.par, data, eps, v=-log(2*(1-0.99999999))){
  mid <- mean(c(i,o))
  
  if(conditions.verify(mid, fix.par, 1, v, data) &&
       conditions.verify(mid, fix.par, 0, v, data)){ i <- mid }
  else{ o <- mid }
  
  if(abs(i-o) < eps){ return(i) }
  else return(dicho.a.bounds(i, o, fix.par, data, eps, v))
}