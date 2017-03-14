\name{tsxtreme-package}
\alias{tsxtreme-package}
\alias{tsxtreme}
\docType{package}
\title{
\packageTitle{tsxtreme}
}
\description{
\packageDescription{tsxtreme}
}
\details{
\packageIndices{tsxtreme}
The Heffernan--Tawn conditional formulation for a stationary time series \eqn{(X_t)} with Laplace marginal distribution states that for a large enough threshold \eqn{u} there exist scale parameters \eqn{-1 \le\alpha_j\le 1} and \eqn{0 \le \beta_j \le 1} such that
\deqn{Pr\left(\frac{X_j-\alpha_j X_0}{(X_j)^{\beta_j}} < z_j, j=1,\ldots,m \mid X_0 > u\right) = H(z_1,\ldots,z_m),}{%
Pr{(X_j-\alpha_j * X_0)/(X_j)^{\beta_j} < z_j, j=1,\ldots,m | X_0 > u} = H(z_1,\ldots,z_m)}

with \eqn{H} non-degenerate; the equality holds exactly only when \eqn{u} tends to infinity.

There are mainly 3 functions provided by this package, which allow estimation of extremal dependence measures and fitting the Heffernan--Tawn model using Dirichlet processes.

\code{\link{depfit}} fits the Heffernan--Tawn model using a Bayesian semi-parametric approach.

\code{\link{thetafit}} computes posterior samples of the threshold-based index of Ledford and Tawn (2003) based on inference in \code{\link{depfit}}.

\code{\link{chifit}} computes posterior samples of the extremal measure of dependence of Coles, Heffernan and Tawn (1999) at any extremal level.

Some corresponding functions using the stepwise approach of Heffernan and Tawn (2004) are also part of the package, namely \code{\link{dep2fit}} and \code{\link{theta2fit}}.

The empirical estimation of the extremal index can be done using \code{\link{thetaruns}} and some basic functions handling the Laplace distribution are also available in \code{\link{dlapl}}.
}
\author{
\packageAuthor{tsxtreme}

Maintainer: \packageMaintainer{tsxtreme}
}
\references{
Coles, S., Heffernan, J. E. and Tawn, J. A. (1999) Dependence measures for extreme value analyses. \emph{Extremes}, \bold{2}, 339--365.

Davison, A. C. and Smith, R. L. (1990) Models for exceedances over high thresholds. \emph{Journal of the Royal Statistical Society Series B}, \bold{52}, 393--442.

Heffernan, J. E. and Tawn, J. A. (2004) A conditional approach for multivariate extreme values. \emph{Journal of the Royal Statistical Society Series B}, \bold{66}, 497--546.

Ledford, W. A. and Tawn, J. A. (2003) Diagnostics for dependence within time series extremes. \emph{Journal of the Royal Statistical Society Series B}, \bold{65}, 521--543.

Lugrin, T., Davison, A. C. and Tawn, J. A. (2016) Bayesian uncertainty management in temporal dependence of extremes. \emph{Extremes}, \bold{19}, 491--515.
}

\keyword{ package }
\seealso{
\code{\link{thetafit}}, \code{\link{chifit}}, \code{\link{depfit}}
}
