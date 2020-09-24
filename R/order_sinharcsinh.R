#' Random Sampling of Order Statistics from a Sinh-Arcsinh Distribution
#'
#'\code{order_sinharcsinh} is used to obtain a random sample of order statistics from a Sinh-Arcsinh Distribution.
#' @param size numeric, represents the size of the sample.
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param mu numeric, represents the location parameter values.
#' @param sigma numeric, represents scale parameter values.
#' @param nu numeric, represents skewness parameter values
#' @param tau numeric, represents kurtosis tau parameter values.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param ... represents others parameters of a Sinh-Arcsinh distribution.
#'
#' @return A list with a random sample of order statistics from a Sinh-Arcsinh Distribution and the value of its join probability density function evaluated in the random sample.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Ribgy, R. and Stasinopoulos, M. (2005) Generalized Additive Models for Location Scale and Shape, Journal of the Royal Statistical Society. Applied Statistics, Series C.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of size 10 of the 3-th order statistics from a Sinh-Arcsinh Distribution
#' order_sinharcsinh(size=10,k=3,mu=0,sigma=1,nu=1,tau=2,n=30)
#' @importFrom gamlss.dist qSHASH dSHASH
#' @importFrom stats rbeta
#' @export order_sinharcsinh

order_sinharcsinh <- function(size,k,mu,sigma,nu,tau,n,...){
  initial <- rbeta(size, k, n + 1 - k)
  sample  <- qSHASH(initial,mu,sigma,nu,tau,...)
  pdf     <- factorial(size)*cumprod(dSHASH(sample,mu,sigma,nu,tau,...))[size]
  return(list(sample=sample,pdf=pdf))
}
