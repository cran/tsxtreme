\name{dep2fit}

\alias{dep2fit}

\title{
Dependence model fit (stepwise)
}
\description{
The conditional Heffernan--Tawn model is used to fit the dependence in time of a stationary series. A standard 2-stage procedure is used.
}
\usage{
dep2fit(ts, u.mar = 0, u.dep,
        lapl = FALSE, method.mar = c("mle","mom","pwm"),
        nlag = 1, conditions = TRUE)
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
 a character string defining the method used to estimate the marginal GPD; either \code{"mle"} for maximum likelihood of \code{"mom"} for method of moments. Defaults to \code{"mle"}.
}
  \item{nlag}{
integer; number of lags to be considered when modelling the dependence in time.
}
  \item{conditions}{
logical; should conditions on \eqn{\alpha} and \eqn{\beta} be set? (see Details) Defaults to \code{TRUE}.
}
}
\details{
Consider a stationary time series \eqn{(X_t)} with Laplace marginal distribution; the fitting procedure consists of fitting
\deqn{X_t = \alpha_t\times x_0 + x_0^{\beta_t}\times Z_t,\quad t=1,\ldots,m,}{X_t = \alpha_t * x_0 + x_0^{\beta_t} * Z_t, t=1,\ldots,m,}
with \eqn{m} the number of lags considered. A likelihood is maximised assuming \eqn{Z_t\sim N(\mu_t, \sigma^2_t)}{Z_t ~ N(\mu_t, \sigma^2_t)}, then an empirical distribution for the \eqn{Z_t} is derived using the estimates of \eqn{\alpha_t} and \eqn{\beta_t} and the relation
\deqn{\hat Z_t = \frac{X_t - \hat\alpha_t\times x_0}{x_0^{\hat\beta_t}}.}{Z_t = (X_t - \alpha_t * x_0) / x_0^{\beta_t}.}

\code{conditions} implements additional conditions suggested by Keef, Papastathopoulos and Tawn (2013) on the ordering of conditional quantiles. These conditions help with getting a consistent fit by shrinking the domain in which \eqn{(\alpha,\beta)} live.
}
\value{
  \item{alpha }{parameter controlling the conditional extremal expectation.}
  \item{beta }{parameter controlling the conditional extremal expectation and variance.}
  \item{res }{empirical residual of the model.}
  \item{pars.se }{vector of length 2 giving the estimated standard errors for \code{alpha} and \code{beta} given by the hessian matrix of the likelihood function used in the first step of the inference procedure.}
}

\seealso{
\code{\link{depfit}}, \code{\link{theta2fit}}
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

## rescale margin
ar <- qlapl(pnorm(ar))

## fit model without constraints...
fit1 <- dep2fit(ts=ar, u.mar = 0.95, u.dep=0.98, conditions=FALSE)
fit1$a; fit1$b

## ...and compare with a fit with constraints
fit2 <- dep2fit(ts=ar, u.mar = 0.95, u.dep=0.98, conditions=TRUE)
fit2$a; fit2$b# should be similar, as true parameters lie well within the constraints
}

\keyword{models}
\keyword{ts}
\keyword{nonparametric}
\keyword{multivariate}