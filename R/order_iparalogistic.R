#' Random Sampling of k-th Order Statistics from a Inverse Paralogistic Distribution
#'
#'\code{order_iparalogistic} is used to obtain a random sample of the k-th order statistic from a Inverse Paralogistic distribution and some associated quantities of interest.
#' @param size numeric, represents the size of the sample.
#' @param shape numeric, represents a first shape parameter value. Must be strictly positive.
#' @param scale numeric, represents scale parameter values. Must be strictly positive.
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param alpha numeric, (1 - alpha) represents the confidence of an interval for the population median of the distribution of the k-th order statistic. Default value is 0.05.
#' @param ... represents others parameters of a Inverse Paralogistic distribution.
#' @return A list with a random sample of order statistics from a Inverse Paralogistic Distribution, the value of its join probability density function evaluated in the random sample and
#' a (1 - alpha) confidence interval for the population median of the distribution of the k-th order statistic.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012), Loss Models, From Data to Decisions, Fourth Edition, Wiley.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of size 10 of the 3-th order statistics from a Inverse Paralogistic Distribution
#' order_iparalogistic(size=10,shape=1.5,scale=0.5,k=3,n=50,alpha=0.02)
#' @importFrom actuar qinvparalogis dinvparalogis
#' @importFrom stats rbeta
#' @export order_iparalogistic

order_iparalogistic <- function(size,k,shape,scale,n,alpha=0.05,...){
  sample  <- qinvparalogis(initial_order(size,k,n),shape,scale,...)
  pdf     <- factorial(size)*cumprod(dinvparalogis(sample,shape,scale,...))[size]
  if(size>5){
    return(list(sample=sample,pdf=pdf,ci_median=interval_median(size,sample,alpha)))
  }
  cat("---------------------------------------------------------------------------------------------\n")
  cat("We cannot report the confidence interval. The size of the sample is less or equal than five.\n")
  return(list(sample=sample,pdf=pdf))
}
