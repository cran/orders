#' Random Sampling of Order Statistics from a Log Gamma G I Distribution
#'
#'\code{order_loggammag1} is used to obtain a random sample of order statistics from a Log Gamma G I Distribution.
#' @param size numeric, represents the size of the sample.
#' @param spec character, represents an specific G distribution. Possible values "norm", "exp","lnorm","chisq".
#' @param a numeric, represents the first shape parameter. Default value is 1.
#' @param b numeric, represents the second shape parameter. Default value is 1.
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param ... represents others parameters of the G distribution.
#' @return A list with a random sample of order statistics from a Log Gamma G I Distribution and the value of its join probability density function evaluated in the random sample.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Naradajah, S. and Rocha, R. (2016) Newdistns: An R Package for New Families of Distributions, Journal of Statistical Software.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Log Gamma Exponential I Distribution
#' order_loggammag1(10,"exp",1,1,k=3,50)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Log Gamma Normal I Distribution
#' order_loggammag1(10,"norm",1,1,k=3,50)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Log Gamma Log-normal I Distribution
#' order_loggammag1(10,"lnorm",1,1,k=3,50)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Log Gamma Chi-square I Distribution
#' order_loggammag1(10,"chisq",1,1,k=3,50,df=3)
#' @importFrom Newdistns qloggammag1 dloggammag1
#' @importFrom stats rbeta
#' @export order_loggammag1

order_loggammag1 <- function(size,spec,a,b,k,n,...){
  initial <- rbeta(size, k, n + 1 - k)
  sample  <- qloggammag1(initial,spec,a,b,...)
  pdf     <- factorial(size)*cumprod(dloggammag1(sample,spec,a,b,...))[size]
  return(list(sample=sample,pdf=pdf))
}
