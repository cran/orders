#' Random Sampling of Order Statistics from a Zero Inflated Poisson Distribution
#'
#'\code{order_zip} is used to obtain a random sample of order statistics from a Zero Inflated Poisson Distribution.
#' @param size numeric, represents the size of the sample.
#' @param mu numeric, represents the location parameter values.
#' @param sigma numeric, represents scale parameter values.
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param ... represents others parameters of a Zero Inflated Poisson distribution.
#' @return A list with a random sample of order statistics from a Zero Inflated Poisson Distribution and the value of its join probability density function evaluated in the random sample.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Ribgy, R. and Stasinopoulos, M. (2005) Generalized Additive Models for Location Scale and Shape, Journal of the Royal Statistical Society. Applied Statistics, Series C.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of size 10 of the 3-th order statistics from a Zero Inflated Poisson Distribution
#' order_zip(size=10,k=3,mu=5,sigma=0.1,n=30)
#' @importFrom gamlss.dist qZIP dZIP
#' @importFrom stats rbeta
#' @export order_zip

order_zip <- function(size,k,mu,sigma,n,...){
  initial <- rbeta(size, k, n + 1 - k)
  sample  <- qZIP(initial,mu,sigma,...)
  pdf     <- factorial(size)*cumprod(dZIP(sample,mu,sigma,...))[size]
  return(list(sample=sample,pdf=pdf))
}
