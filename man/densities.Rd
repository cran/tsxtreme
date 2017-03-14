\name{Laplace}

\alias{dlapl}
\alias{plapl}
\alias{qlapl}
\alias{rlapl}

\title{
The Laplace Distribution
}
\description{
Density, distribution function, quantile function and random generation for the Laplace distribution with location parameter \code{loc} and scale parameter \code{scale}.
}
\usage{
dlapl(x, loc = 0, scale = 1, log = FALSE)
plapl(q, loc = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
qlapl(p, loc = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
rlapl(n, loc = 0, scale = 1)
}

\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of samples to generate.}
  \item{loc}{vector of location parameters.}
  \item{scale}{vector of scale parameters. These must be positive.}
  \item{lower.tail}{logical; if TRUE (default), probabilities are \eqn{Pr(X\le x)}, otherwise \eqn{Pr(X>x)}.}
  \item{log,log.p}{logical; if TRUE, probabilities \eqn{p} are given as \eqn{\log(p)}{log(p)}.}
}
\details{
If \code{loc} or \code{scale} are not specified, they assume the default values of 0 and 1 respectively.

The Laplace distribution has density
\deqn{f(x) = \exp(- |x-\mu|/\sigma)/(2\sigma)}
where \eqn{\mu} is the location parameter and \eqn{\sigma} is the scale parameter.
}
\value{
\code{dlapl} gives the density, \code{plapl} gives the distribution function, \code{qlapl} gives the quantile function, and \code{rlapl} generates random deviates.

The length of the result is determined by \code{n} in \code{rlapl}, and is the maximum of the lengths of the numerical arguments for the other functions. Standard \code{R} vector operations are to be assumed.

If \code{sd}=0, the limit as \code{sd} decreases to 0 is returned, i.e., a point mass at \code{loc}. The case \code{sd}<0 is an error and nothing is returned.
}

\section{Warning}{
Some checks are done previous to standard evaluation, but vector computations have not yet been tested thoroughly! Typically vectors not having lengths multiple of each other return an error.
}

\seealso{
\link[stats]{dexp} for the exponential distribution which is the positive part of the Laplace distribution.
}

\examples{
## evaluate the density function on a grid of values
x  <- seq(from=-5, to=5, by=0.1)
fx <- dlapl(x, loc=1, scale=.5)

## generate random samples of a mixture of Laplace distributions
rnd <- rlapl(1000, loc=c(-5,-3,2), scale=0.5)

## an alternative:
rnd <- runif(1000)
rnd <- qlapl(rnd, loc=c(-5,-3,2), scale=0.5)

## integrate the Laplace density on [a,b]
a <- -1
b <- 7
integral <- plapl(b)-plapl(a)
}

\keyword{datagen}
\keyword{distribution}
