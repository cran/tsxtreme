\name{theta2fit}

\alias{theta2fit}

\title{
Fit time series extremes
}
\description{
Appropriate marginal transforms are done before the fit using standard procedures, before the dependence model is fitted to the data. Then the measure of dependence \eqn{\theta(x,m)} is derived using a method described in Eastoe and Tawn (2012).
}
\usage{
theta2fit(ts, lapl = FALSE, nlag = 1, R = 1000,
          u.mar = 0, u.dep, probs = seq(u.dep, 0.9999, length.out = 30),
          method.mar = c("mle","mom","pwm"), method = c("prop", "MCi"),
          silent = FALSE,
          R.boot = 0, block.length = m * 5, levels = c(.025,.975))
}

\arguments{
  \item{ts}{
numeric vector; time series to be fitted.
}
  \item{lapl}{
logical; is \code{ts} on the Laplace scale already? The default (FALSE) assumes unknown marginal distribution.
}
  \item{nlag}{
integer; number of lags to be considered when modelling the dependence in time.
}
  \item{R}{
integer; the number of samples used for estimating \eqn{\theta(x,m)}.
}
  \item{u.mar}{
marginal threshold; used when transforming the time series to Laplace scale if \code{lapl} is FALSE; not used otherwise.
}
  \item{u.dep}{
dependence threshold; level above which the dependence is modelled. \code{u.dep} can be lower than \code{u.mar}.
}
  \item{probs}{
vector of probabilities; the values of \eqn{x} for which to evaluate \eqn{\theta(x,m)}.
}
  \item{method.mar}{
  a character string defining the method used to estimate the marginal GPD; either \code{"mle"} for maximum likelihood of \code{"mom"} for method of moments or \code{"pwm"} for probability weighted moments methods. Defaults to \code{"mle"}.
}
  \item{method}{
a character string defining the method used to estimate the dependence measure; either \code{"prop"} for proportions or \code{"MCi"} for Monte Carlo integration (see Details).
}
  \item{silent}{
logical (\code{FALSE}); verbosity.
}
  \item{R.boot}{
integer; the number of samples used for the block bootstrap for the confidence intervals.
}
  \item{block.length}{
integer; the block length used for the block-bootstrapped confidence intervals.
}
  \item{levels}{
vector of probabilities; the quantiles of the bootstrap distribution of the extremal measure to be computed.
}
}
\details{
The standard procedure (\code{method="prop"}) to estimating probabilities from a Heffernan-Tawn fit best illustrated in the bivariate context (\eqn{Y\mid X>u}{Y | X > u}):

1. sample \eqn{X} from an exponential distribution above \eqn{v \ge u},

2. sample \eqn{Z} (residuals) from their empirical distribution,

3. compute \eqn{Y} using the relation \eqn{Y = \alpha\times X + X^\beta\times Z}{Y = \alpha * X + X^\beta * Z},

4. estimate \eqn{Pr(X > v_x, Y > v_y)} by calculating the proportion \eqn{p} of \eqn{Y} samples above \eqn{v_y} and multiply \eqn{p} with the marginal survival distribution evaluated at \eqn{v_x}.

With \code{method="MCi"} a Monte Carlo integration approach is used, where the survivor distribution of \eqn{Z} is evaluated at pseudo-residuals of the form
\deqn{\frac{v_y - \alpha\times X}{X^\beta},}{(v_y - \alpha * X)/X^\beta,}
where \eqn{X} is sampled from an exponential distribution above \eqn{v_x}. Taking the mean of these survival probabilities, we get the Monte Carlo equivalent of \eqn{p} in the proportion approach.
}
\value{
 List containing:
  \item{depfit }{an object of class '\link{stepfit}'}
  \item{probs }{\code{probs}}
  \item{levels }{\code{probs} transformed to original scale of \code{ts}}
  \item{theta }{a matrix with \emph{proportion} or \emph{Monte Carlo} estimates of \eqn{\theta(x,m)}. Rows correspond to values in \code{probs}, columns are point estimates and bootstrap quantiles}
}

\seealso{
\code{\link{dep2fit}}, \code{\link{thetafit}}, \code{\link{thetaruns}}
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
plot(ar, type="l")
plot(density(ar))
grid <- seq(-3,3,0.01)
lines(grid, dnorm(grid), col="blue")

## rescale the margin (focus on dependence)
ar <- qlapl(pnorm(ar))

## fit the data
fit <- theta2fit(ts=ar, u.mar=0.95, u.dep=0.98)

## plot theta(x,1)
plot(fit)
abline(h=1, lty="dotted")
}

\keyword{models}
\keyword{ts}
\keyword{nonparametric}
\keyword{multivariate}