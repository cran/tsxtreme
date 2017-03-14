\name{bayesfit}
\alias{bayesfit}
\alias{is.bayesfit}
\alias{summary.bayesfit}
\alias{plot.bayesfit}
\alias{print.bayesfit}

\title{Traces from MCMC output}

\description{Test or show objects of class "bayesfit".}

\usage{
is.bayesfit(x)
}

\arguments{
\item{x}{an arbitrary \R object.}
}

\details{
Default plot shows samples of residual densities (\code{which==1}), residual distribution with credible interval (\eqn{5\%} and \eqn{95\%} posterior quantiles; \code{which==2}), and joint posterior distribution of \eqn{\alpha} and \eqn{\beta} (\code{which==3}) for each lag successively. \code{which} can be any composition of 1,2 and 3.
}

\value{
An object of class "bayesfit" is a list containing MCMC traces for:
\item{a,b }{Heffernan-Tawn parameters.}
\item{sd,mean,w }{standard deviations, means and weights of the mixture components.}
\item{prec }{precision parameter of the Dirichlet process.}
\item{ci }{auxiliary variable; components' indices for each observation.}
\item{noo }{number of observations in each mixture component.}
\item{noc }{number of non-empty components in the mixture.}
\item{prop.sd }{standard deviations of the proposal distributions for \code{a} and \code{b}.}
And \code{len}, the length of the traces, i.e., the number of iterations saved.
}

\seealso{
\code{\link{bayesparams}}, \code{\link{stepfit}}
}

\keyword{classes}
