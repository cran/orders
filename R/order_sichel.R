#' Random Sampling of Order Statistics from a Sichel Distribution
#'
#'\code{order_sichel} is used to obtain a random sample of order statistics from a Sichel Distribution.
#' @param size numeric, represents the size of the sample.
#' @param mu numeric, represents the location parameter values.
#' @param sigma numeric, represents scale parameter values.
#' @param nu numeric, represents skewness parameter values
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param ... represents others parameters of a Skew normal type 1 distribution.
#' @return A list with a random sample of order statistics from a Sichel Distribution and the value of its join probability density function evaluated in the random sample.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Ribgy, R. and Stasinopoulos, M. (2005) Generalized Additive Models for Location Scale and Shape, Journal of the Royal Statistical Society. Applied Statistics, Series C.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of size 10 of the 3-th order statistics from a Sichel Distribution
#' order_sichel(size=10,k=3,mu=5,sigma=1,nu=1,n=30)
#' @importFrom gamlss.dist qSICHEL dSICHEL
#' @importFrom stats rbeta
#' @export order_sichel

order_sichel <- function(size,k,mu,sigma,nu,n,...){
  initial <- rbeta(size, k, n + 1 - k)
  sample  <- qSICHEL(initial,mu,sigma,nu,...)
  pdf     <- factorial(size)*cumprod(dSICHEL(sample,mu,sigma,nu,...))[size]
  return(list(sample=sample,pdf=pdf))
}
