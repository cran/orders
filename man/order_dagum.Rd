% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/order_dagum.R
\name{order_dagum}
\alias{order_dagum}
\title{Random Sampling of k-th Order Statistics from a Dagum Distribution}
\usage{
order_dagum(size, k, shape1, shape2, scale, n, p = 0.5, alpha = 0.05, ...)
}
\arguments{
\item{size}{numeric, represents the size of the sample.}

\item{k}{numeric, represents the k-th smallest value from a sample.}

\item{shape1}{numeric, represents a first shape parameter value. Must be strictly positive.}

\item{shape2}{numeric, represents a second shape parameter value. Must be strictly positive.}

\item{scale}{numeric, represents scale parameter values. Must be strictly positive.}

\item{n}{numeric, represents the size of the sample to compute the order statistic from.}

\item{p}{numeric, represents the 100p percentile for the distribution of the k-th order statistic. Default value is population median, p = 0.5.}

\item{alpha}{numeric, (1 - alpha) represents the confidence of an interval for the population percentile p of the distribution of the k-th order statistic. Default value is 0.05.}

\item{...}{represents others parameters of a Dagum distribution.}
}
\value{
A list with a random sample of order statistics from a Dagum Distribution, the value of its join probability density function evaluated in the random sample
and an approximate (1 - alpha) confidence interval for the population percentile p of the distribution of the k-th order statistic.
}
\description{
\code{order_dagum} is used to obtain a random sample of the k-th order statistic from a Dagum distribution and some associated quantities of interest.
}
\examples{
library(orders)
# A sample of size 10 of the 3-th order statistics from a Dagum Distribution
order_dagum(size=10,shape1=1,shape2=2,scale=1,k=3,n=50,p=0.5,alpha=0.02)
}
\references{
Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.

Kleiber, C. and Kotz, S. (2003). Statistical Size Distributions in Economics and Actuarial Sciences, Hoboken, NJ, USA: Wiley-Interscience.
}
\author{
Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
}
