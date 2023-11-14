#' Random Sampling of k-th Order Statistics from a Inverse Pareto Distribution
#'
#'\code{order_invpareto} is used to obtain a random sample of the k-th order statistic from a Inverse Pareto distribution and some associated quantities of interest.
#' @param size numeric, represents the size of the sample.
#' @param shape1 numeric, represents a first shape parameter value. Must be strictly positive.
#' @param scale numeric, represents scale parameter values. Must be strictly positive.
#' @param k numeric, represents the k-th smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param p numeric, represents the 100p percentile for the distribution of the k-th order statistic. Default value is population median, p = 0.5.
#' @param alpha numeric, (1 - alpha) represents the confidence of an interval for the population percentile p of the distribution of the k-th order statistic. Default value is 0.05.
#' @param ... represents others parameters of a Inverse Pareto distribution.
#' @return A list with a random sample of order statistics from a Inverse Pareto Distribution, the value of its join probability density function evaluated in the random sample
#' and an approximate (1 - alpha) confidence interval for the population percentile p of the distribution of the k-th order statistic.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012), Loss Models, From Data to Decisions, Fourth Edition, Wiley.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of size 10 of the 3-th order statistics from a Inverse Pareto Distribution
#' order_invpareto(size=10,shape1=0.75,scale=0.5,k=3,n=50,p=0.5,alpha=0.02)
#' @importFrom actuar qinvpareto dinvpareto
#' @importFrom stats rbeta
#' @export order_invpareto

order_invpareto <- function(size,k,shape1,scale,n,p=0.5,alpha=0.05,...){
  sample  <- qinvpareto(initial_order(size,k,n),shape1,scale,...)
  pdf     <- factorial(size)*cumprod(dinvpareto(sample,shape1,scale,...))[size]
  log_pdf     <- sum(log(2:size)) + sum(log(dinvpareto(sample,shape1,scale,...)))
  if(size>5){
    int_perc_est <- interval_percentile_est(p,size,sample,alpha)
    return(list(sample = sample,
                pdf = pdf,
                log_pdf = log_pdf,
                point_percentile_est = point_percentile_est(p,size,sample),
                confidence_percentile_est = int_perc_est[1:2],
                aprox_coverage_prob = int_perc_est[3]))
  }
  cat("---------------------------------------------------------------------------------------------\n")
  cat("We cannot report the confidence interval. The size of the sample is less or equal than five.\n")
  return(list(sample=sample,pdf=pdf))
}
