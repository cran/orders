#' Random Sampling of Order Statistics from a Poisson-inverse Gaussian Distribution
#'
#'\code{order_pig} is used to obtain a random sample of order statistics from a Poisson-inverse Gaussian Distribution.
#' @param size numeric, represents the size of the sample.
#' @param mu numeric, represents the location parameter values.
#' @param sigma numeric, represents scale parameter values.
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param ... represents others parameters of a Poisson-inverse Gaussian distribution.
#' @return A list with a random sample of order statistics from a Poisson-inverse Gaussian Distribution and the value of its join probability density function evaluated in the random sample.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Ribgy, R. and Stasinopoulos, M. (2005) Generalized Additive Models for Location Scale and Shape, Journal of the Royal Statistical Society. Applied Statistics, Series C.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of size 10 of the 3-th order statistics from a Poisson-inverse Gaussian Distribution
#' order_pig(size=10,k=3,mu=6,sigma=1,n=30)
#' @importFrom gamlss.dist qPIG dPIG
#' @importFrom stats rbeta
#' @export order_pig

order_pig <- function(size,k,mu,sigma,n,...){
  initial <- rbeta(size, k, n + 1 - k)
  sample  <- qPIG(initial,mu,sigma,...)
  pdf     <- factorial(size)*cumprod(dPIG(sample,mu,sigma,...))[size]
  return(list(sample=sample,pdf=pdf))
}
