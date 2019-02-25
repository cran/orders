#' Random Sampling of Order Statistics from a Exponentiated Exponential Poisson G Distribution
#'
#'\code{order_eepg} is used to obtain a random sample of order statistics from a Exponentiated Exponential Poisson G Distribution.
#' @param size numeric, represents the size of the sample.
#' @param spec character, represents an specific distribution. Possible values "norm", "exp","lnorm","chisq".
#' @param lambda numeric, represents a scale parameter. Default value is 1.
#' @param a numeric, represents the shape parameter. Default value is 1.
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param ... represents others parameters of the G distribution.
#' @return A list with a random sample of order statistics from a Exponentiated Exponential Poisson G Distribution and the value of its join probability density function evaluated in the random sample.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Naradajah, S. and Rocha, R. (2016) Newdistns: An R Package for New Families of Distributions, Journal of Statistical Software.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of order statistics from a Exponentiated Exponential Poisson Exponential Distribution
#' order_eepg(10,"exp",1,1,1,50)
#' # A sample of order statistics from a Exponentiated Exponential Poisson Normal Distribution
#' order_eepg(10,"norm",1,1,1,50)
#' # A sample of order statistics from a Exponentiated Exponential Poisson Log-normal Distribution
#' order_eepg(10,"lnorm",1,1,1,50)
#' # A sample of order statistics from a Exponentiated Exponential Poisson Chi-square Distribution
#' order_eepg(10,"chisq",1,1,1,50,df=3)
#' @importFrom Newdistns qeepg deepg
#' @importFrom stats rbeta
#' @export order_eepg

order_eepg <- function(size,spec,lambda,a,k,n,...){
  initial <- rbeta(size, k, n + 1 - k)
  sample  <- qeepg(initial,spec,lambda,a,...)
  pdf     <- factorial(size)*cumprod(deepg(sample,spec,lambda,a,...))[size]
  return(list(sample=sample,pdf=pdf))
}
