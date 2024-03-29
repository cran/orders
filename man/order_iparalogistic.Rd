% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/order_iparalogistic.R
\name{order_iparalogistic}
\alias{order_iparalogistic}
\title{Random Sampling of k-th Order Statistics from a Inverse Paralogistic Distribution}
\usage{
order_iparalogistic(size, k, shape, scale, n, p = 0.5, alpha = 0.05, ...)
}
\arguments{
\item{size}{numeric, represents the size of the sample.}

\item{k}{numeric, represents the k-th smallest value from a sample.}

\item{shape}{numeric, represents a first shape parameter value. Must be strictly positive.}

\item{scale}{numeric, represents scale parameter values. Must be strictly positive.}

\item{n}{numeric, represents the size of the sample to compute the order statistic from.}

\item{p}{numeric, represents the 100p percentile for the distribution of the k-th order statistic. Default value is population median, p = 0.5.}

\item{alpha}{numeric, (1 - alpha) represents the confidence of an interval for the population percentile p of the distribution of the k-th order statistic. Default value is 0.05.}

\item{...}{represents others parameters of a Inverse Paralogistic distribution.}
}
\value{
A list with a random sample of order statistics from a Inverse Paralogistic Distribution, the value of its join probability density function evaluated in the random sample
and an approximate (1 - alpha) confidence interval for the population percentile p of the distribution of the k-th order statistic.
}
\description{
\code{order_iparalogistic} is used to obtain a random sample of the k-th order statistic from a Inverse Paralogistic distribution and some associated quantities of interest.
}
\examples{
library(orders)
# A sample of size 10 of the 3-th order statistics from a Inverse Paralogistic Distribution
order_iparalogistic(size=10,shape=1.5,scale=0.5,k=3,n=50,p=0.5,alpha=0.02)
}
\references{
Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.

Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012), Loss Models, From Data to Decisions, Fourth Edition, Wiley.
}
\author{
Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
}
