#' Random Sampling of Order Statistics from a Skew normal type 1 Distribution
#'
#'\code{order_snormal1} is used to obtain a random sample of order statistics from a Skew normal type 1 Distribution.
#' @param size numeric, represents the size of the sample.
#' @param mu numeric, represents the location parameter values.
#' @param sigma numeric, represents scale parameter values.
#' @param nu numeric, represents skewness parameter values
#' @param tau numeric, represents kurtosis tau parameter values.
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param ... represents others parameters of a Skew normal type 1 distribution.
#' @return A list with a random sample of order statistics from a Skew normal type 1 Distribution and the value of its join probability density function evaluated in the random sample.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Ribgy, R. and Stasinopoulos, M. (2005) Generalized Additive Models for Location Scale and Shape, Journal of the Royal Statistical Society. Applied Statistics, Series C.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of size 10 of the 3-th order statistics from a Skew normal type 1 Distribution
#' order_snormal1(size=10,mu=0,sigma=1,nu=0,tau=2,k=3,n=50)
#' @importFrom gamlss.dist qST1 dST1
#' @importFrom stats rbeta
#' @export order_snormal1

order_snormal1 <- function(size,k,mu,sigma,nu,tau,n,...){
  initial <- rbeta(size, k, n + 1 - k)
  sample  <- qST1(initial,mu,sigma,nu,tau,...)
  pdf     <- factorial(size)*cumprod(dST1(sample,mu,sigma,nu,tau,...))[size]
  return(list(sample=sample,pdf=pdf))
}
