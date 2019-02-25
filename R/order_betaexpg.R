#' Random Sampling of Order Statistics from a Beta Exponential G Distribution
#'
#'\code{order_betaexpg} is used to obtain a random sample of order statistics from a Beta Exponential G Distribution.
#' @param size numeric, represents the size of the sample.
#' @param spec character, represents an specific distribution. Possible values "norm", "exp","lnorm","chisq".
#' @param lambda numeric, represents the first shape parameter. Default value is 1.
#' @param a numeric, represents the second shape parameter. Default value is 1.
#' @param b numeric, represents the third shape parameter. Default value is 1.
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param ... represents others parameters of the G distribution.
#' @return A list with a random sample of order statistics from a Beta Exponential G Distribution and the value of its join probability density function evaluated in the random sample.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Naradajah, S. and Rocha, R. (2016) Newdistns: An R Package for New Families of Distributions, Journal of Statistical Software.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of order statistics from a Beta Exp Exponential Distribution
#' order_betaexpg(10,"exp",1,1,1,1,50)
#' # A sample of order statistics from a Beta Exp Normal Distribution
#' order_betaexpg(10,"norm",1,1,1,1,50)
#' # A sample of order statistics from a Beta Exp Log-normal Distribution
#' order_betaexpg(10,"lnorm",1,1,1,1,50)
#' # A sample of order statistics from a Beta Exp Chi-square  Distribution
#' order_betaexpg(10,"chisq",1,1,1,1,50,df=3)
#' @importFrom Newdistns qbetaexpg dbetaexpg
#' @importFrom stats rbeta
#' @export order_betaexpg

order_betaexpg <- function(size,spec,lambda,a,b,k,n,...){
  initial      <- rbeta(size, k, n + 1 - k)
  sample       <- qbetaexpg(initial,spec,lambda,a,b,...)
  pdf          <- factorial(size)*cumprod(dbetaexpg(sample,spec,...))[size]
  return(list(sample=sample,pdf=pdf))
}
