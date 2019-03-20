#' Random Sampling of Order Statistics from a Truncated-Exponential Skew-Symmetric G Distribution
#'
#'\code{order_tessg} is used to obtain a random sample of order statistics from a Truncated-Exponential Skew-Symmetric G Distribution.
#' @param size numeric, represents the size of the sample.
#' @param spec character, represents an specific G distribution. Possible values "norm", "exp","lnorm","chisq".
#' @param lambda numeric, represents the skewness parameter. Default value is 1.
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param ... represents others parameters of the G distribution.
#' @return A list with a random sample of order statistics from a Truncated-Exponential Skew-Symmetric G Distribution and the value of its join probability density function evaluated in the random sample.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Naradajah, S. and Rocha, R. (2016) Newdistns: An R Package for New Families of Distributions, Journal of Statistical Software.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Truncated-Exponential Skew-Symmetric Exponential Distribution
#' order_tessg(10,"exp",1,k=3,50)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Truncated-Exponential Skew-Symmetric Normal Distribution
#' order_tessg(10,"norm",1,k=3,50)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Truncated-Exponential Skew-Symmetric Log-normal Distribution
#' order_tessg(10,"lnorm",1,k=3,50)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Truncated-Exponential Skew-Symmetric Chi-square Distribution
#' order_tessg(10,"chisq",1,k=3,50,df=3)
#' @importFrom Newdistns qtessg dtessg
#' @importFrom stats rbeta
#' @export order_tessg

order_tessg <- function(size,spec,lambda,k,n,...){
  initial <- rbeta(size, k, n + 1 - k)
  sample  <- qtessg(initial,spec,lambda,...)
  pdf     <- factorial(size)*cumprod(dtessg(sample,spec,lambda,...))[size]
  return(list(sample=sample,pdf=pdf))
}
