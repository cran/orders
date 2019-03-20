#' Random Sampling of Order Statistics from a Weibull G Distribution
#'
#'\code{order_weibullg} is used to obtain a random sample of order statistics from a Weibull G Distribution.
#' @param size numeric, represents the size of the sample.
#' @param spec character, represents an specific G distribution. Possible values "norm", "exp","lnorm","chisq".
#' @param beta numeric, represents the scale parameter. Default value is 1.
#' @param c numeric, represents the shape parameter. Default value is 1.
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param ... represents others parameters of the G distribution.
#' @return A list with a random sample of order statistics from a Weibull G Distribution and the value of its join probability density function evaluated in the random sample.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Naradajah, S. and Rocha, R. (2016) Newdistns: An R Package for New Families of Distributions, Journal of Statistical Software.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of size 10 of the 3-th order statistics from a Weibull Exponential Distribution
#' order_weibullg(10,"exp",beta=1,c=1,k=3,n=50)
#' # A sample of size 10 of the 3-th order statistics from a Weibull Normal Distribution
#' order_weibullg(10,"norm",beta=1,c=1,k=3,n=50)
#' # A sample of size 10 of the 3-th order statistics from a Weibull Log-normal Distribution
#' order_weibullg(10,"lnorm",beta=1,c=1,k=3,n=50)
#' # A sample of size 10 of the 3-th order statistics from a Weibull Chi-square Distribution
#' order_weibullg(10,"chisq",beta=1,c=1,k=3,n=50,df=3)
#' @importFrom Newdistns qweibullg dweibullg
#' @importFrom stats rbeta
#' @export order_weibullg

order_weibullg <- function(size,spec,beta,c,k,n,...){
  initial <- rbeta(size, k, n + 1 - k)
  sample  <- qweibullg(initial,spec,beta,c,...)
  pdf     <- factorial(size)*cumprod(dweibullg(sample,spec,beta,c,...))[size]
  return(list(sample=sample,pdf=pdf))
}
