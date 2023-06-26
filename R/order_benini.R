#' Random Sampling of k-th Order Statistics from a Benini Distribution
#'
#'\code{order_benini} is used to obtain a random sample of the k-th order statistic from a Benini distribution and some associated quantities of interest.
#' @param size numeric, represents the size of the sample.
#' @param shape1 numeric, represents a first shape parameter value. Must be strictly positive.
#' @param scale numeric, represents scale parameter values. Must be strictly positive.
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param alpha numeric, (1 - alpha) represents the confidence of an interval for the population median of the distribution of the k-th order statistic. Default value is 0.05.
#' @param ... represents others parameters of a Benini distribution.
#' @return A list with a random sample of order statistics from a Benini Distribution, the value of its join probability density function evaluated in the random sample and
#' a (1 - alpha) confidence interval for the population median of the distribution of the k-th order statistic.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Kleiber, C. and Kotz, S. (2003). Statistical Size Distributions in Economics and Actuarial Sciences, Hoboken, NJ, USA: Wiley-Interscience.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of size 10 of the 3-th order statistics from a Benini Distribution
#' order_benini(size=10,shape1=0.75,scale=1,k=3,n=50,alpha=0.02)
#' @importFrom VGAM qbenini dbenini
#' @importFrom stats rbeta
#' @export order_benini

order_benini <- function(size,k,shape1,scale,n,alpha=0.05,...){
  sample  <- qbenini(initial_order(size,k,n),shape1,scale,...)
  pdf     <- factorial(size)*cumprod(dbenini(sample,shape1,scale,...))[size]
  if(size>5){
    return(list(sample=sample,pdf=pdf,ci_median=interval_median(size,sample,alpha)))
  }
  cat("---------------------------------------------------------------------------------------------\n")
  cat("We cannot report the confidence interval. The size of the sample is less or equal than five.\n")
  return(list(sample=sample,pdf=pdf))
}
