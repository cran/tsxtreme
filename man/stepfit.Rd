\name{stepfit}
\alias{stepfit}
\alias{is.stepfit}
\alias{summary.stepfit}
\alias{plot.stepfit}
\alias{print.stepfit}

\title{
Estimates from stepwise fit
}
\description{
Create, test or show objects of class "stepfit".
}
\usage{
stepfit()

is.stepfit(x)
}

\arguments{
\item{x}{an arbitrary \R object.}
}

\value{
An object of class "stepfit" is a list containing:
\item{a,b }{Heffernan--Tawn parameters.}
\item{res }{fitted residuals.}
\item{pars.se }{estimated standard error of \code{a} and \code{b}.}
}

\seealso{
\code{\link{bayesfit}}, \code{\link{depmeasure}}
}

\keyword{classes}
