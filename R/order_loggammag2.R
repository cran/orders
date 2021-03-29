#' Random Sampling of k-th Order Statistics from a Log Gamma G II Distribution
#'
#'\code{order_loggammag2} is used to obtain a random sample of the k-th order statistic from a Log Gamma G II distribution.
#' @param size numeric, represents the size of the sample.
#' @param spec character, represents an specific G distribution. Possible values "norm", "exp","lnorm","chisq".
#' @param a numeric, represents the first shape parameter. Default value is 1.
#' @param b numeric, represents the second shape parameter. Default value is 1.
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param alpha numeric, (1 - alpha) represents the confidence of an interval for the population median of the distribution of the k-th order statistic. Default value is 0.05.
#' @param ... represents others parameters of the G distribution.
#' @return A list with a random sample of order statistics from a Log Gamma G II Distribution, the value of its join probability density function evaluated in the random sample and
#' a (1 - alpha) confidence interval for the population median of the distribution of the k-th order statistic.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Naradajah, S. and Rocha, R. (2016) Newdistns: An R Package for New Families of Distributions, Journal of Statistical Software.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Log Gamma Exponential II Distribution
#' order_loggammag2(10,"exp",1,1,k=3,50,alpha=0.02)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Log Gamma Normal II Distribution
#' order_loggammag2(10,"norm",1,1,k=3,50)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Log Gamma Log-normal II Distribution
#' order_loggammag2(10,"lnorm",1,1,k=3,50)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Log Gamma Chi-square II Distribution
#' order_loggammag2(10,"chisq",1,1,k=3,50,df=3)
#' @importFrom Newdistns qloggammag2 dloggammag2
#' @importFrom stats rbeta
#' @export order_loggammag2

order_loggammag2 <- function(size,spec,a,b,k,n,alpha=0.05,...){
  sample  <- qloggammag2(initial_order(size,k,n),spec,a,b,...)
  pdf     <- factorial(size)*cumprod(dloggammag2(sample,spec,a,b,...))[size]
  if(size>5){
    return(list(sample=sample,pdf=pdf,ci_median=interval_median(size,sample,alpha)))
  }
  cat("---------------------------------------------------------------------------------------------\n")
  cat("We cannot report the confidence interval. The size of the sample is less or equal than five.\n")
  return(list(sample=sample,pdf=pdf))
}
