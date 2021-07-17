#' Random Sampling of k-th Order Statistics from a Feller-Pareto Distribution
#'
#'\code{order_fellerpareto} is used to obtain a random sample of the k-th order statistic from a Feller-Pareto distribution and some associated quantities of interest.
#' @param size numeric, represents the size of the sample.
#' @param min numeric, represents the lower bound of the support of the distribution.
#' @param shape1 numeric, represents a first shape parameter value. Must be strictly positive.
#' @param shape2 numeric, represents a second shape parameter value. Must be strictly positive.
#' @param shape3 numeric, represents a third shape parameter value. Must be strictly positive.
#' @param scale numeric, represents scale parameter values. Must be strictly positive.
#' @param k numeric, represents the Kth smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param alpha numeric, (1 - alpha) represents the confidence of an interval for the population median of the distribution of the k-th order statistic. Default value is 0.05.
#' @param ... represents others parameters of a Feller-Pareto distribution.
#' @return A list with a random sample of order statistics from a Feller-Pareto Distribution, the value of its join probability density function evaluated in the random sample and
#' a (1 - alpha) confidence interval for the population median of the distribution of the k-th order statistic.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012), Loss Models, From Data to Decisions, Fourth Edition, Wiley.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of size 10 of the 3-th order statistics from a Feller-Pareto Distribution
#' order_fellerpareto(size=10,min=0.5,shape1=0.75,shape2=1,shape3=1.25,scale=0.5,k=3,n=50,alpha=0.02)
#' @importFrom actuar qfpareto dfpareto
#' @importFrom stats rbeta
#' @export order_fellerpareto

order_fellerpareto <- function(size,k,min,shape1,shape2,shape3,scale,n,alpha=0.05,...){
  sample  <- qfpareto(initial_order(size,k,n),min,shape1,shape2,shape3,scale,...)
  pdf     <- factorial(size)*cumprod(dfpareto(sample,min,shape1,shape2,shape3,scale,...))[size]
  if(size>5){
    return(list(sample=sample,pdf=pdf,ci_median=interval_median(size,sample,alpha)))
  }
  cat("---------------------------------------------------------------------------------------------\n")
  cat("We cannot report the confidence interval. The size of the sample is less or equal than five.\n")
  return(list(sample=sample,pdf=pdf))
}
