#' Random Sampling of k-th Order Statistics from a Modified Beta G Distribution
#'
#'\code{order_mbetag} is used to obtain a random sample of k-th order statistic from a Modified Beta G distribution.
#' @param size numeric, represents the size of the sample.
#' @param spec character, represents an specific G distribution. Possible values "norm", "exp","lnorm","chisq".
#' @param beta numeric, represents the scale parameter. Default value is 1.
#' @param a numeric, represents a shape parameter must be positive. Default value is 1.
#' @param b numeric, represents a shape parameter must be positive. Default value is 1.
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param alpha numeric, (1 - alpha) represents the confidence of an interval for the population median of the distribution of the k-th order statistic. Default value is 0.05.
#' @param ... represents others parameters of the G distribution.
#' @return A list with a random sample of order statistics from a Modified Beta G Distribution, the value of its join probability density function evaluated in the random sample and
#' a (1 - alpha) confidence interval for the population median of the distribution of the k-th order statistic.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Naradajah, S. and Rocha, R. (2016) Newdistns: An R Package for New Families of Distributions, Journal of Statistical Software.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Modified Beta Exponential Distribution
#' order_mbetag(10,"exp",1,1,1,k=3,50,alpha=0.02)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Modified Beta Normal Distribution
#' order_mbetag(10,"norm",1,1,1,k=3,50)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Modified Beta Log-normal Distribution
#' order_mbetag(10,"lnorm",1,1,1,k=3,50)
#' @importFrom Newdistns qmbetag dmbetag
#' @importFrom stats rbeta
#' @export order_mbetag

order_mbetag <- function(size,spec,beta,a,b,k,n,alpha=0.05,...){
  sample  <- qmbetag(initial_order(size,k,n),spec,beta,a,b,...)
  pdf     <- factorial(size)*cumprod(dmbetag(sample,spec,beta,a,b,...))[size]
  if(size>5){
    return(list(sample=sample,pdf=pdf,ci_median=interval_median(size,sample,alpha)))
  }
  cat("---------------------------------------------------------------------------------------------\n")
  cat("We cannot report the confidence interval. The size of the sample is less or equal than five.\n")
  return(list(sample=sample,pdf=pdf))
}
