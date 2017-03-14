\name{depfit}

\alias{depfit}

\title{
Dependence model fit
}

\description{
Bayesian semiparametrics are used to fit the Heffernan--Tawn model to time series. Options are available to impose a structure in time on the model.
}

\usage{
depfit(ts, u.mar = 0, u.dep=u.mar,
    lapl = FALSE, method.mar=c("mle","mom","pwm"),  nlag = 1,
    par = bayesparams(),
    submodel = c("fom", "none", "ugm"))
}

\arguments{
  \item{ts}{
numeric vector; time series to be fitted.
}
  \item{u.mar}{
  marginal threshold; used when transforming the time series to Laplace scale.
  }
  \item{u.dep}{
  dependence threshold; level above which the dependence is modelled. \code{u.dep} can be lower than \code{u.mar}.
  }
  \item{lapl}{
  logical; is \code{ts} on the Laplace scale already? The default (FALSE) assumes unknown marginal distribution.
}
  \item{method.mar}{
  a character string defining the method used to estimate the marginal GPD; either \code{"mle"} for maximum likelihood or \code{"mom"} for method of moments or \code{"pwm"} for probability weighted moments. Defaults to \code{"mle"}.
}
  \item{nlag}{
  integer; number of lags to be considered when modelling the dependence in time.
  }
  \item{par}{
  an object of class 'bayesparams'.
  }
  \item{submodel}{
a character string; "fom" for \emph{first order Markov}, "none" for \emph{no particular time structure}, or "ugm" for \emph{univariate Gaussian mixture} (see details).
}
}
\details{
\code{submodel} can be \code{"fom"} to impose a first order Markov structure on the model parameters \eqn{\alpha_j} and \eqn{\beta_j} (see \code{\link{thetafit}} for more details); it can take the value \code{"none"} to impose no particular structure in time; it can also be \code{"ugm"} which can be applied to density estimation, as it corresponds to setting \eqn{\alpha=\beta=0} (see examples).
}

\value{
An object of class 'bayesfit' with elements:

  \item{a }{posterior trace of \eqn{\alpha}.}
  \item{b }{posterior trace of \eqn{\beta}.}
  \item{sd }{posterior trace of the components' standard deviations.}
  \item{mean }{posterior trace of the components' means.}
  \item{w }{posterior trace of the components' assigned weights.}
  \item{prec }{posterior trace of the precision parameter.}
  \item{noo }{posterior trace of the number of observations per component.}
  \item{noc }{posterior trace of the number of components containing at least one observation.}
  \item{prop.sd }{trace of proposal standard deviations in the 5+3 regions of the adaption scheme for \eqn{\alpha} and \eqn{\beta}.}
  \item{len }{length of the returned traces.}
}

\seealso{
\code{\link{thetafit}}, \code{\link{chifit}}
}
\examples{
## generate data from an AR(1)
## with Gaussian marginal distribution
n   <- 10000
dep <- 0.5
ar    <- numeric(n)
ar[1] <- rnorm(1)
for(i in 2:n)
  ar[i] <- rnorm(1, mean=dep*ar[i-1], sd=1-dep^2)
  
## rescale the margin
ar <- qlapl(pnorm(ar))

## fit the data
params <- bayesparams()
params$maxit <- 100# bigger numbers would be
params$burn  <- 10 # more sensible...
params$thin  <- 4
fit <- depfit(ts=ar, u.mar=0.95, u.dep=0.98, par=params)

########
## density estimation with submodel=="ugm"
data <- MASS::galaxies/1e3
dens <- depfit(ts=data, par=params, submodel="ugm")
}

\keyword{models}
\keyword{iteration}
\keyword{multivariate}
\keyword{ts}
\keyword{nonparametric}
