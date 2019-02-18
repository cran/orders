#' Random Sampling of Order Statistics from a Beta Extended G Distribution
#'
#'\code{order_beg} is used to obtain a random sample of order statistics from a Beta Extended G Distribution.
#' @param size numeric, represents the size of the sample.
#' @param spec character, represents an specific distribution. Possible values "norm", "exp","lnorm".
#' @param alpha numeric, represents the scale parameter. Default value is 1.
#' @param a numeric, represents the first shape parameter. Default value is 1.
#' @param b numeric, represents the second shape parameter. Default value is 1.
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Naradajah, S. and Rocha, R. (2016) Newdistns: An R Package for New Families of Distributions, Journal of Statiscal Software.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' order_beg(10,"exp",1,1,1,1,50)
#' order_beg(10,"norm",1,1,1,1,50)
#' order_beg(10,"lnorm",1,1,1,1,50)
#' @importFrom Newdistns qbeg
#' @importFrom stats rbeta
#' @export order_beg

order_beg <- function(size,spec,alpha,a,b,k,n){
  initial <- rbeta(size, k, n + 1 - k)
  return(qbeg(initial,spec,alpha,a,b))
}
