% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/order_tessg.R
\name{order_tessg}
\alias{order_tessg}
\title{Random Sampling of k-th Order Statistics from a Truncated-Exponential Skew-Symmetric G Distribution}
\usage{
order_tessg(size, spec, lambda, k, n, p = 0.5, alpha = 0.05, ...)
}
\arguments{
\item{size}{numeric, represents the size of the sample.}

\item{spec}{character, represents an specific G distribution. Possible values "norm", "exp","lnorm","chisq".}

\item{lambda}{numeric, represents the skewness parameter. Default value is 1.}

\item{k}{numeric, represents the k-th smallest value from a sample.}

\item{n}{numeric, represents the size of the sample to compute the order statistic from.}

\item{p}{numeric, represents the 100p percentile for the distribution of the k-th order statistic. Default value is population median, p = 0.5.}

\item{alpha}{numeric, (1 - alpha) represents the confidence of an interval for the population percentile p of the distribution of the k-th order statistic. Default value is 0.05.}

\item{...}{represents others parameters of the G distribution.}
}
\value{
A list with a random sample of order statistics from a Truncated-Exponential Skew-Symmetric G Distribution, the value of its join probability density function evaluated in the random sample
and an approximate (1 - alpha) confidence interval for the population percentile p of the distribution of the k-th order statistic.
}
\description{
\code{order_tessg} is used to obtain a random sample of the k-th order statistic from a Truncated-Exponential Skew-Symmetric G distribution.
}
\examples{
library(orders)
# A sample of size 10 of the 3-th order statistics from
# a Truncated-Exponential Skew-Symmetric Exponential Distribution
order_tessg(10,"exp",1,k=3,n=50,p=0.5,alpha=0.02)
# A sample of size 10 of the 3-th order statistics from
# a Truncated-Exponential Skew-Symmetric Normal Distribution
order_tessg(10,"norm",1,k=3,n=50,p=0.5,)
# A sample of size 10 of the 3-th order statistics from
# a Truncated-Exponential Skew-Symmetric Log-normal Distribution
order_tessg(10,"lnorm",1,k=3,n=50,p=0.5,)
# A sample of size 10 of the 3-th order statistics from
# a Truncated-Exponential Skew-Symmetric Chi-square Distribution
order_tessg(10,"chisq",1,k=3,n=50,p=0.5,df=3)
}
\references{
Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.

Naradajah, S. and Rocha, R. (2016) Newdistns: An R Package for New Families of Distributions, Journal of Statistical Software.
}
\author{
Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
}
