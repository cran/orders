#' Random Sampling of Order Statistics from a Generalized Beta G Distribution
#'
#'\code{order_gbg} is used to obtain a random sample of order statistics from a Generalized Beta G Distribution.
#' @param size numeric, represents the size of the sample.
#' @param spec character, represents an specific G distribution. Possible values "norm", "exp","lnorm","chisq".
#' @param a numeric, represents the first shape parameter. Default value is 1.
#' @param b numeric, represents the second shape parameter. Default value is 1.
#' @param c numeric, represents the third shape parameter. Default value is 1.
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param ... represents others parameters of the G distribution.
#' @return A list with a random sample of order statistics from a Generalized Beta G Distribution and the value of its join probability density function evaluated in the random sample.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Naradajah, S. and Rocha, R. (2016) Newdistns: An R Package for New Families of Distributions, Journal of Statistical Software.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Generalized Beta Exponential Distribution
#' order_gbg(10,"exp",1,1,1,k=3,50)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Generalized Beta Normal Distribution
#' order_gbg(10,"norm",1,1,1,k=3,50)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Generalized Beta Log-normal Distribution
#' order_gbg(10,"lnorm",1,1,1,k=3,50)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Generalized  Beta Chi-square Distribution
#' order_gbg(10,"chisq",1,1,1,k=3,50,df=3)
#' @importFrom Newdistns qgbg dgbg
#' @importFrom stats rbeta
#' @export order_gbg

order_gbg <- function(size,spec,a,b,c,k,n,...){
  initial      <- rbeta(size, k, n + 1 - k)
  sample       <- qgbg(initial,spec,a,b,c,...)
  pdf          <- factorial(size)*cumprod(dgbg(sample,spec,...))[size]
  return(list(sample=sample,pdf=pdf))
}
