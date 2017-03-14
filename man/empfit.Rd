\name{thetaruns}

\alias{thetaruns}

\title{
Runs estimator
}
\description{
Compute the empirical estimator of the extremal index using the runs method (Smith & Weissman, 1994, JRSSB).
}
\usage{
thetaruns(ts, lapl = FALSE, nlag = 1,
          u.mar = 0, probs = seq(u.mar, 0.995, length.out = 30),
          method.mar = c("mle", "mom", "pwm"),
          R.boot = 0, block.length = (nlag+1) * 5, levels = c(0.025, 0.975))
}

\arguments{
  \item{ts}{
a vector, the time series for which to estimate the threshold-based extremal index \eqn{\theta(x,m)}, with \eqn{x} a probability level and \eqn{m} a run-length (see details).
}
  \item{lapl}{
logical; is \code{ts} on the Laplace scale already? The default (FALSE) assumes unknown marginal distribution. 
}
  \item{nlag}{
the run-length; an integer larger or equal to 1.
}
  \item{u.mar}{
marginal threshold (probability); used when transforming the time series to Laplace scale if \code{lapl} is FALSE; if \code{lapl} is TRUE, it is nevertheless used when bootstrapping, since the bootstrapped series generally do not have Laplace marginal distributions.
}
  \item{probs}{
vector of probabilities; the values of \eqn{x} for which to evaluate \eqn{\theta(x,m)}.
}
  \item{method.mar}{
a character string defining the method used to estimate the marginal GPD; either \code{"mle"} for maximum likelihood or \code{"mom"} for method of moments or \code{"pwm"} for probability weighted moments methods. Defaults to \code{"mle"}.
}
  \item{block.length}{
integer; the block length used for the block-bootstrapped confidence intervals.
}
  \item{R.boot}{
integer; the number of samples used for the block bootstrap.
}
  \item{levels}{
vector of probabilites; the quantiles of the posterior distribution of the extremal index \eqn{\theta(x,m)} to output.
}
}
\details{
Consider a stationary time series \eqn{(X_t)}. A characterisation of the extremal index is
\deqn{\theta(x,m) = Pr(X_1\le x,\ldots,X_m\le x \mid X_0\ge x).}{\theta(x,m) = Pr(X_1\le x,\dots,X_m\le x | X_0\ge x).}
In the limit when \eqn{x} and \eqn{m} tend to \eqn{\infty} appropriately, \eqn{\theta} corresponds to the asymptotic inverse mean cluster size. It also links the generalised extreme value distribution of the independent series \eqn{(Y_t)}, with the same marginal distribution as \eqn{(X_t)},
\deqn{G_Y(z)=G_X^\theta(z),}
with \eqn{G_X} and \eqn{G_Y} the extreme value distributions of \eqn{(X_t)} and \eqn{(Y_t)} respectively.

\code{nlag} corresponds to the \emph{run-length} \eqn{m} and \code{probs} is a set of values for \eqn{x}.
The \emph{runs} estimator is computed, which consists of counting the proportion of clusters to the number of exceedances of a threshold \eqn{x}; two exceedances of the threshold belong to different clusters if there are at least \eqn{m+1} non-exceedances inbetween.
}
\value{
An object of class '\code{\link{depmeasure}}' containing:
  \item{theta }{matrix; estimates of the extremal index \eqn{\theta(x,m)} with rows corresponding to the \code{probs} values of \eqn{x} and the columns to the runs estimate and the chosen \code{levels}-quantiles of the bootstrap distribution.}
  \item{nbr.exc }{numeric vector; number of exceedances for each threshold corresponding to the elements in \code{probs}.}
  \item{probs }{\code{probs}.}
  \item{levels }{numeric vector; \code{probs} converted to the original scale of \code{ts}.}
  \item{nlag }{\code{nlag}.}
}

\seealso{
\code{\link{theta2fit}}, \code{\link{thetafit}}
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
## transform to Laplace scale
ar <- qlapl(pnorm(ar))
## compute empirical estimate
theta <- thetaruns(ts=ar, u.mar=.95, probs=c(.95,.98,.99))
## output
plot(theta, ylim=c(.2,1))
abline(h=1, lty="dotted")
}

\keyword{models}
\keyword{ts}
\keyword{nonparametric}
