#' Random Sampling of Order Statistics from a Exponentiated Exponential Poisson G Distribution
#'
#'\code{order_eepg} is used to obtain a random sample of order statistics from a Exponentiated Exponential Poisson G Distribution.
#' @param size numeric, represents the size of the sample.
#' @param spec character, represents an specific distribution. Possible values "norm", "exp","lnorm".
#' @param lambda numeric, represents a scale parameter. Default value is 1.
#' @param a numeric, represents the shape parameter. Default value is 1.
#' @param k numeric, represents the error bound of the approximation of P(t). Default values is 0.01.
#' @param n numeric, represents the error bound of the approximation of P(t). Default values is 0.01.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Naradajah, S. and Rocha, R. (2016) Newdistns: An R Package for New Families of Distributions, Journal of Statiscal Software.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' order_eepg(10,"exp",1,1,1,50)
#' order_eepg(10,"norm",1,1,1,50)
#' order_eepg(10,"norm",1,1,1,50)
#' @importFrom Newdistns qeepg
#' @importFrom stats rbeta
#' @export order_eepg

order_eepg <- function(size,spec,lambda,a,k,n){
  initial <- rbeta(size, k, n + 1 - k)
  return(qeepg(initial,spec,lambda,a))
}
