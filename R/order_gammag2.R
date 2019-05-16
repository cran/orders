#' Random Sampling of Order Statistics from a Gamma G II Distribution
#'
#'\code{order_gammag2} is used to obtain a random sample of order statistics from a Gamma G II Distribution.
#' @param size numeric, represents the size of the sample.
#' @param spec character, represents an specific G distribution. Possible values "norm", "exp","lnorm","chisq".
#' @param a numeric, represents the first shape parameter. Default value is 1.
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param ... represents others parameters of the G distribution.
#' @return A list with a random sample of order statistics from a Gamma G II Distribution and the value of its join probability density function evaluated in the random sample.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Naradajah, S. and Rocha, R. (2016) Newdistns: An R Package for New Families of Distributions, Journal of Statistical Software.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of size 10 of the 3-th order statistics from a Gamma Exponential II Distribution
#' order_gammag2(10,"exp",1,k=3,50)
#' # A sample of size 10 of the 3-th order statistics from a Gamma Normal II Distribution
#' order_gammag2(10,"norm",1,k=3,50)
#' # A sample of size 10 of the 3-th order statistics from a Gamma Log-normal II Distribution
#' order_gammag2(10,"lnorm",1,k=3,50)
#' # A sample of size 10 of the 3-th order statistics from a Gamma Chi-square II Distribution
#' order_gammag2(10,"chisq",1,k=3,50,df=3)
#' @importFrom Newdistns qgammag2 dgammag2
#' @importFrom stats rbeta
#' @export order_gammag2

order_gammag2 <- function(size,spec,a,k,n,...){
  initial <- rbeta(size, k, n + 1 - k)
  sample  <- qgammag2(initial,spec,a,...)
  pdf     <- factorial(size)*cumprod(dgammag2(sample,spec,a,...))[size]
  return(list(sample=sample,pdf=pdf))
}
