#' Random Sampling of Order Statistics from a Exponentiated Kumaraswamy G Distribution
#'
#'\code{order_expkumg} is used to obtain a random sample of order statistics from a Exponentiated Kumaraswamy G Distribution.
#' @param size numeric, represents the size of the sample.
#' @param spec character, represents an specific G distribution. Possible values "norm", "exp","lnorm","chisq".
#' @param a numeric, represents the first shape parameter. Default value is 1.
#' @param b numeric, represents the second shape parameter. Default value is 1.
#' @param c numeric, represents the third shape parameter. Default value is 1.
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param ... represents others parameters of the G distribution.
#' @return A list with a random sample of order statistics from a Exponentiated Kumaraswamy G Distribution and the value of its join probability density function evaluated in the random sample.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Naradajah, S. and Rocha, R. (2016) Newdistns: An R Package for New Families of Distributions, Journal of Statistical Software.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Exponentiated Kumaraswamy Exponential Distribution
#' order_expkumg(10,"exp",1,1,1,k=3,50)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Exponentiated Kumaraswamy Normal Distribution
#' order_expkumg(10,"norm",1,1,1,k=3,50)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Exponentiated Kumaraswamy Log-normal Distribution
#' order_expkumg(10,"lnorm",1,1,1,k=3,50)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Exponentiated Kumaraswamy Chi-square Distribution
#' order_expkumg(10,"chisq",1,1,1,k=3,50,df=3)
#' @importFrom Newdistns qexpkumg dexpkumg
#' @importFrom stats rbeta
#' @export order_expkumg

order_expkumg <- function(size,spec,a,b,c,k,n,...){
  initial <- rbeta(size, k, n + 1 - k)
  sample  <- qexpkumg(initial,spec,a,b,c,...)
  pdf     <- factorial(size)*cumprod(dexpkumg(sample,spec,a,b,c,...))[size]
  return(list(sample=sample,pdf=pdf))
}
