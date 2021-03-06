% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/order_pig.R
\name{order_pig}
\alias{order_pig}
\title{Random Sampling of Order Statistics from a Poisson-inverse Gaussian Distribution}
\usage{
order_pig(size, k, mu, sigma, n, ...)
}
\arguments{
\item{size}{numeric, represents the size of the sample.}

\item{k}{numeric, represents the Kth smallest value from a sample.}

\item{mu}{numeric, represents the location parameter values.}

\item{sigma}{numeric, represents scale parameter values.}

\item{n}{numeric, represents the size of the sample to compute the order statistic from.}

\item{...}{represents others parameters of a Poisson-inverse Gaussian distribution.}
}
\value{
A list with a random sample of order statistics from a Poisson-inverse Gaussian Distribution and the value of its join probability density function evaluated in the random sample.
}
\description{
\code{order_pig} is used to obtain a random sample of order statistics from a Poisson-inverse Gaussian Distribution.
}
\examples{
library(orders)
# A sample of size 10 of the 3-th order statistics from a Poisson-inverse Gaussian Distribution
order_pig(size=10,k=3,mu=6,sigma=1,n=30)
}
\references{
Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.

Ribgy, R. and Stasinopoulos, M. (2005) Generalized Additive Models for Location Scale and Shape, Journal of the Royal Statistical Society. Applied Statistics, Series C.
}
\author{
Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
}
