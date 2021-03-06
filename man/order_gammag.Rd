% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/order_gammag.R
\name{order_gammag}
\alias{order_gammag}
\title{Random Sampling of Order Statistics from a Gamma Uniform G Distribution}
\usage{
order_gammag(size, spec, a, k, n, ...)
}
\arguments{
\item{size}{numeric, represents the size of the sample.}

\item{spec}{character, represents an specific G distribution. Possible values "norm", "exp","lnorm","chisq".}

\item{a}{numeric, represents the first shape parameter. Default value is 1.}

\item{k}{numeric, represents the Kth smallest value from a sample.}

\item{n}{numeric, represents the size of the sample to compute the order statistic from.}

\item{...}{represents others parameters of the G distribution.}
}
\value{
A list with a random sample of order statistics from a Gamma Uniform G Distribution and the value of its join probability density function evaluated in the random sample.
}
\description{
\code{order_gammag} is used to obtain a random sample of order statistics from a Gamma Uniform G Distribution.
}
\examples{
library(orders)
# A sample of size 10 of the 3-th order statistics from a Gamma Uniform Exponential Distribution
order_gammag(10,"exp",1,k=3,50)
# A sample of size 10 of the 3-th order statistics from a Gamma Uniform Normal Distribution
order_gammag(10,"norm",1,k=3,50)
# A sample of size 10 of the 3-th order statistics from a Gamma Uniform Log-normal Distribution
order_gammag(10,"lnorm",1,k=3,50)
# A sample of size 10 of the 3-th order statistics from a Gamma Uniform Chi-square Distribution
order_gammag(10,"chisq",1,k=3,50,df=3)
}
\references{
Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.

Naradajah, S. and Rocha, R. (2016) Newdistns: An R Package for New Families of Distributions, Journal of Statistical Software.
}
\author{
Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
}
