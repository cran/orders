#' Random Sampling of Order Statistics from a Beta G Distribution
#'
#'\code{order_betag} is used to obtain a random sample of order statistics from a Beta G Distribution.
#' @param size numeric, represents the size of the sample.
#' @param spec character, represents an specific distribution. Possible values "norm", "exp","lnorm","chisq".
#' @param a numeric, represents the first shape parameter. Default value is 1.
#' @param b numeric, represents the first shape parameter. Default value is 1.
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param ... represents others parameters of the G distribution.
#' @return A list with a random sample of order statistics from a Beta G Distribution and the value of its join probability density function evaluated in the random sample.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Naradajah, S. and Rocha, R. (2016) Newdistns: An R Package for New Families of Distributions, Journal of Statistical Software.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of size 10 of order statistics from a Beta Exponential Distribution
#' order_betag(10,"exp",1,1,1,50)
#' # A sample of size 10 of order statistics from a Beta Normal Distribution
#' order_betag(10,"norm",1,1,1,50)
#' # A sample of size 10 of order statistics from a Beta Log-normal Distribution
#' order_betag(10,"lnorm",1,1,1,50)
#' # A sample of size 10 of order statistics from a Beta Chis-square Distribution
#' order_betag(10,"chisq",1,1,1,50,df=3)
#' @importFrom Newdistns qbetag dbetag
#' @importFrom stats rbeta
#' @export order_betag

order_betag <- function(size,spec,a,b,k,n,...){
  initial   <- rbeta(size, k, n + 1 - k)
  sample    <- qbetag(initial,spec,a,b,...)
  pdf       <- factorial(size)*cumprod(dbetag(sample,spec,a,b,...))[size]
  return(list(sample=sample,pdf=pdf))
}
