\name{depmeasures}
\alias{depmeasures}
\alias{thetafit}
\alias{chifit}

\title{
Estimate dependence measures
}

\description{
 Appropriate marginal transforms are done before the fit using standard procedures, before the dependence model is fitted to the data. Then the posterior distribution of a measure of dependence is derived. \code{thetafit} gives posterior samples for the extremal index \eqn{\theta(x,m)} and \code{chifit} does the same for the coefficient of extremal dependence \eqn{\chi_m(x)}.
}

\usage{
thetafit(ts, lapl = FALSE, nlag = 1,
      R = 1000, S = 500,
      u.mar = 0, u.dep,
      probs = seq(u.dep, 0.9999, length.out = 30),
      method.mar = c("mle", "mom","pwm"), method = c("prop", "MCi"),
      silent = FALSE,
      fit = TRUE, prev.fit=bayesfit(), par = bayesparams(),
      submodel = c("fom", "none"), levels=c(.025,.975))

chifit(ts, lapl = FALSE, nlag = 1,
      R = 1000, S = 500,
      u.mar = 0, u.dep,
      probs = seq(u.dep, 0.9999, length.out = 30),
      method.mar = c("mle", "mom","pwm"), method = c("prop", "MCi"),
      silent = FALSE,
      fit = TRUE, prev.fit=bayesfit(), par = bayesparams(),
      submodel = c("fom", "none"), levels=c(.025,.975))
}

\arguments{
  \item{ts}{
a vector, the time series for which to estimate the extremal index \eqn{\theta(x,m)} or the coefficient of extremal dependence \eqn{\chi_m(x)}, with \eqn{x} a probability level and \eqn{m} a run-length (see details).
}
  \item{lapl}{
logical; \code{TRUE} indicates that \code{ts} has a marginal Laplace distribution. If \code{FALSE} (default), \code{method.mar} is used to transform the marginal distribution of \code{ts} to Laplace.
}
  \item{nlag}{
the run-length; an integer larger or equal to 1.
}
  \item{R}{
the number of samples per MCMC iteration drawn from the sampled posterior distributions; used for the estimation of the dependence measure.
}
  \item{S}{
the number of posterior distributions sampled to be used for the estimation of the dependence measure.
}
  \item{u.mar}{
probability; threshold used for marginal transformation if \code{lapl} is \code{FALSE}. Not used otherwise.
}
  \item{u.dep}{
probability; threshold used for the extremal dependence model.
}
  \item{probs}{
vector of probabilities; the values of \eqn{x} for which to evaluate \eqn{\theta(x,m)} or \eqn{\chi_m(x)}.
}
  \item{method.mar}{
a character string defining the method used to estimate the marginal GPD; either \code{"mle"} for maximum likelihood of \code{"mom"} for method of moments or \code{"pwm"} for probability weighted moments methods. Defaults to \code{"mle"}.
}
  \item{method}{
a character string defining the method used to estimate the dependence measure; either \code{"prop"} for proportions or \code{"MCi"} for Monte Carlo integration (see details).
}
  \item{silent}{
logical (\code{FALSE}); verbosity.
}
  \item{fit}{
logical; \code{TRUE} means that the dependence model must be fitted and the values in \code{par} are used. Otherwise the result from a previous call to \code{\link{depfit}}.
}
  \item{prev.fit}{
an object of class 'bayesfit'. Needed if \code{fit} is FALSE. Typically returned by a previous call to \code{\link{depfit}}.
}
  \item{par}{
an object of class '\code{\link{bayesparams}}' to be used for the fit of dependence model.
}
  \item{submodel}{
a character string, either \code{"fom"} for \emph{first order Markov} or \code{"none"} for no specification.
}
\item{levels}{
vector of probabilites; the quantiles of the posterior distribution of the extremal measure to be computed.
}
}
\details{
The sub-asymptotic extremal index is defined as
\deqn{\theta(x,m) = Pr(X_1 < x,\ldots,X_m < x | X_0 > x),}
whose limit as \eqn{x} and \eqn{m} go to \eqn{\infty}appropriately is the extremal index \eqn{\theta}. The extremal index can be interpreted as the inverse of the asymptotic mean cluster size (see \code{\link{thetaruns})}.

The sub-asymptotic coefficient of extremal dependence is
\deqn{\chi_m(x) = Pr(X_m > x | X_0 > x),}
whose limit \eqn{\chi} defines asymptotic dependence (\eqn{\chi > 0}) or asymptotic independence (\eqn{\chi = 0}).

Both types of extremal dependence measures can be estimated either using a

* proportion method (\code{method == "prop"}), sampling from the conditional probability given \eqn{X_0 > x} and counting the proportion of sampled points falling in the region of interest, or

* Monte Carlo integration (\code{method == "MCi"}), sampling replicates from the marginal exponential tail distribution and evaluating the conditional tail distribution in these replicates, then taking their mean as an approximation of the integral.

\code{submodel == "fom"} imposes a first order Markov structure to the model, namely a geometrical decrease in \eqn{\alpha} and a constant \eqn{\beta} across lags, i.e. \eqn{\alpha_j = \alpha^j} and \eqn{\beta_j = \beta}, \eqn{j=1,\ldots,m}.
}
\value{
An object of class 'depmeasure', containing a subset of:
  \item{bayesfit }{An object of class 'bayesfit'}
  \item{theta }{An array with dimensions \code{m} \eqn{\times}{*} \code{length(probs)} \eqn{\times}{*} (2+\code{length(levels)}), with the last dimension listing the posterior mean and median, and the \code{level} posterior quantiles}
  \item{distr }{An array with dimensions \code{m} \eqn{\times}{*} \code{length(probs)} \eqn{\times}{*} \code{S}; posterior samples of \code{theta}}
  \item{chi }{An array with dimensions \code{m} \eqn{\times}{*} \code{length(probs)} \eqn{\times}{*} (2+\code{length(levels)}), with the last dimension listing the posterior mean and median, and the \code{level} posterior quantiles}
  \item{probs }{\code{probs}}
  \item{levels }{\code{probs} transformed to original scale of \code{ts}}
}

\seealso{
\code{\link{depfit}}, \code{\link{theta2fit}}, \code{\link{thetaruns}}
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
params <- bayesparams()
params$maxit <- 100 # bigger numbers would be
params$burn  <- 10  # more sensible...
params$thin  <- 4
theta <- thetafit(ts=ar, R=500, S=100, u.mar=0.95, u.dep=0.98,
                  probs = c(0.98, 0.999), par=params)
## or, same thing in two steps to control fit output before computing theta:
fit <- depfit(ts=ar, u.mar=0.95, u.dep=0.98, par=params)
plot(fit)
theta <- thetafit(ts=ar, R=500, S=100, u.mar=0.95, u.dep=0.98,
                  probs = c(0.98, 0.999), fit=FALSE, prev.fit=fit)
}

\keyword{models}
\keyword{iteration}
\keyword{multivariate}
\keyword{ts}
\keyword{nonparametric}
