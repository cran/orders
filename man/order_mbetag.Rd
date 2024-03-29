% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/order_mbetag.R
\name{order_mbetag}
\alias{order_mbetag}
\title{Random Sampling of k-th Order Statistics from a Modified Beta G Distribution}
\usage{
order_mbetag(size, spec, beta, a, b, k, n, p = 0.5, alpha = 0.05, ...)
}
\arguments{
\item{size}{numeric, represents the size of the sample.}

\item{spec}{character, represents an specific G distribution. Possible values "norm", "exp","lnorm","chisq".}

\item{beta}{numeric, represents the scale parameter. Default value is 1.}

\item{a}{numeric, represents a shape parameter must be positive. Default value is 1.}

\item{b}{numeric, represents a shape parameter must be positive. Default value is 1.}

\item{k}{numeric, represents the k-th smallest value from a sample.}

\item{n}{numeric, represents the size of the sample to compute the order statistic from.}

\item{p}{numeric, represents the 100p percentile for the distribution of the k-th order statistic. Default value is population median, p = 0.5.}

\item{alpha}{numeric, (1 - alpha) represents the confidence of an interval for the population percentile p of the distribution of the k-th order statistic. Default value is 0.05.}

\item{...}{represents others parameters of the G distribution.}
}
\value{
A list with a random sample of order statistics from a Modified Beta G Distribution, the value of its join probability density function evaluated in the random sample
and an approximate (1 - alpha) confidence interval for the population percentile p of the distribution of the k-th order statistic.
}
\description{
\code{order_mbetag} is used to obtain a random sample of k-th order statistic from a Modified Beta G distribution.
}
\examples{
library(orders)
# A sample of size 10 of the 3-th order statistics from
# a Modified Beta Exponential Distribution
order_mbetag(10,"exp",1,1,1,k=3,n=50,p=0.5,alpha=0.02)
# A sample of size 10 of the 3-th order statistics from
# a Modified Beta Normal Distribution
order_mbetag(10,"norm",1,1,1,k=3,n=50,p=0.5,)
# A sample of size 10 of the 3-th order statistics from
# a Modified Beta Log-normal Distribution
order_mbetag(10,"lnorm",1,1,1,k=3,n=50,p=0.5)
}
\references{
Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.

Naradajah, S. and Rocha, R. (2016) Newdistns: An R Package for New Families of Distributions, Journal of Statistical Software.
}
\author{
Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
}
