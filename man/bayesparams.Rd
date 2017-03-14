\name{bayesparams}
\alias{bayesparams}
\alias{is.bayesparams}
\alias{summary.bayesparams}
\alias{print.bayesparams}

\title{
Parameters for the semi-parametric approach
}
\description{
Create, test or show objects of class "bayesparams".
}
\usage{
bayesparams(prop.a = 0.02, prop.b = 0.02,
  prior.mu = c(0, 10), prior.nu = c(2, 1/2), prior.eta = c(2, 2),
  trunc = 100, comp.saved = 15, maxit = 30000,
  burn = 5000, thin = 1,
  adapt = 5000, batch.size = 125,
  mode = 1)

is.bayesparams(x)
}

\arguments{
  \item{prop.a,
    prop.b}{standard deviation for the Gaussian proposal of the Heffernan--Tawn parameters.}
  \item{prior.mu}{mean and standard deviation of the Gaussian prior for the components' means.}
  \item{prior.nu}{shape and rate of the inverse gamma prior for the components' variances.}
  \item{prior.eta}{shape and scale of the gamma prior for the precision parameter of the Dirichlet process.}
  \item{trunc}{integer; value of the truncation for the approximation of the infinite sum in the stick-breaking representation.}
  \item{comp.saved}{number of first components to be saved and returned.}
  \item{maxit}{maximum number of iterations.}
  \item{burn}{number of first iterations to discard.}
  \item{thin}{positive integer; spacing between iterations to be saved. Default is 1, i.e., all iterations are saved.}
  \item{adapt}{integer; number of iterations during which an adaption algorithm is applied to the proposal variances of \eqn{\alpha} and \eqn{\beta}.}
  \item{batch.size}{size of batches used in the adaption algorithm. It has no effect if \code{adapt==0}.}
  \item{mode}{verbosity; 0 for debug mode, 1 (default) for standard output, and 2 for silent.}
  \item{x}{an arbitrary \R object.}
}

\details{
\code{prop.a} is a vector of length 5 with the standard deviations for each region of the RAMA for the (Gaussian) proposal for \eqn{\alpha}. If a scalar is given, 5 identical values are assumed.

\code{prop.b} is a vector of length 3 with the standard deviations for each region of the RAMA for the (Gaussian) proposal for \eqn{\beta}. If a scalar is provided, 3 identical values are assumed.

\code{comp.saved} has no impact on the calculations: its only purpose is to prevent from storing huge amounts of empty components.

The regional adaption scheme targets a \eqn{0.44} acceptance probability. It splits \eqn{[-1;1]} in \eqn{5} regions for \eqn{\alpha} and \eqn{[0;1]} in \eqn{3} regions for \eqn{\beta}. The decision to increase/decrease the proposal standard deviation is based on the past \code{batch.size} MCMC iterations, so too low values yield inefficient adaption, while too large values yield slow adaption.

Default values for the hyperparameters are chosen to get reasonably uninformative priors.}

\seealso{
\code{\link{bayesfit}}, \code{\link{depmeasure}}
}
\examples{
is.bayesparams(bayesparams()) # TRUE
## use defaults, change max number of iteration of MCMC
par <- bayesparams(maxit=1e5)
}

\keyword{classes}
