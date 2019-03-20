#' Random Sampling of Order Statistics from a Marshall Olkin G Distribution
#'
#'\code{order_mog} is used to obtain a random sample of order statistics from a Marshall Olkin G Distribution.
#' @param size numeric, represents the size of the sample.
#' @param spec character, represents an specific G distribution. Possible values "norm", "exp","lnorm","chisq".
#' @param beta numeric, represents the scale parameter. Default value is 1.
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param ... represents others parameters of the G distribution.
#' @return A list with a random sample of order statistics from a Marshall Olkin G Distribution and the value of its join probability density function evaluated in the random sample.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Naradajah, S. and Rocha, R. (2016) Newdistns: An R Package for New Families of Distributions, Journal of Statistical Software.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Marshall Olkin Exponential Distribution
#' order_mog(10,"exp",1,k=3,50)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Marshall Olkin Normal Distribution
#' order_mog(10,"norm",1,k=3,50)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Marshall Olkin Log-normal Distribution
#' order_mog(10,"lnorm",1,k=3,50)
#' @importFrom Newdistns qmog dmog
#' @importFrom stats rbeta
#' @export order_mog

order_mog <- function(size,spec,beta,k,n,...){
  initial <- rbeta(size, k, n + 1 - k)
  sample  <- qmog(initial,spec,beta,...)
  pdf     <- factorial(size)*cumprod(dmog(sample,spec,beta,...))[size]
  return(list(sample=sample,pdf=pdf))
}
