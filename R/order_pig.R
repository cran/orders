#' Random Sampling of k-th Order Statistics from a Poisson-inverse Gaussian Distribution
#'
#'\code{order_pig} is used to obtain a random sample of the k-th order statistic from a Poisson-inverse Gaussian distribution and some associated quantities of interest.
#' @param size numeric, represents the size of the sample.
#' @param mu numeric, represents the location parameter values.
#' @param sigma numeric, represents scale parameter values.
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param alpha numeric, (1 - alpha) represents the confidence of an interval for the population median of the distribution of the k-th order statistic. Default value is 0.05.
#' @param ... represents others parameters of a Poisson-inverse Gaussian distribution.
#' @return A list with a random sample of order statistics from a Poisson-inverse Gaussian Distribution, the value of its join probability density function evaluated in the random sample and
#' a (1 - alpha) confidence interval for the population median of the distribution of the k-th order statistic.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Ribgy, R. and Stasinopoulos, M. (2005) Generalized Additive Models for Location Scale and Shape, Journal of the Royal Statistical Society. Applied Statistics, Series C.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of size 10 of the 3-th order statistics from a Poisson-inverse Gaussian Distribution
#' # order_pig(size=20,k=5,mu=6,sigma=1,n=30,alpha=0.02)
#' @importFrom gamlss.dist qPIG dPIG
#' @importFrom stats rbeta

order_pig <- function(size,k,mu,sigma,n,alpha=0.05,...){
  sample  <- qPIG(initial_order(size,k,n),mu,sigma,...)
  pdf     <- factorial(size)*cumprod(dPIG(sample,mu,sigma,...))[size]
  if(size>5){
    return(list(sample=sample,pdf=pdf,ci_median=interval_median(size,sample,alpha)))
  }
  cat("---------------------------------------------------------------------------------------------\n")
  cat("We cannot report the confidence interval. The size of the sample is less or equal than five.\n")
  return(list(sample=sample,pdf=pdf))
}
