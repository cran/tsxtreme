\name{depmeasure}
\alias{depmeasure}
\alias{is.depmeasure}
\alias{summary.depmeasure}
\alias{plot.depmeasure}
\alias{print.depmeasure}

\title{Dependence measures estimates}
\description{
Test or show objects of class "depmeasure".
}
\usage{
is.depmeasure(x)
}

\arguments{
  \item{x}{an arbitrary \R object.}
}

\value{
An object of class 'depmeasure' is a list which contains:
\item{fit }{an object of class '\link{bayesfit}'}
\item{distr }{an array with the samples used for the estimation.}
\item{probs,levels }{points  ---probability and original scale respectively--- at which the dependence measure is estimated}
Depending on the dependence measure, \code{theta} or \code{chi}, a matrix with levels on row-entries and mean, median and specified quantiles of the posterior distribution of theta or chi respectively.
}

\seealso{
\code{\link{depmeasures}}
}

\keyword{classes}
