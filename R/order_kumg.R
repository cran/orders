#' Random Sampling of k-th Order Statistics from a Kumaraswamy G Distribution
#'
#'\code{order_expkumg} is used to obtain a random sample of he k-th order statistic from a Kumaraswamy G distribution.
#' @param size numeric, represents the size of the sample.
#' @param spec character, represents an specific G distribution. Possible values "norm", "exp","lnorm","chisq".
#' @param a numeric, represents the first shape parameter. Default value is 1.
#' @param b numeric, represents the second shape parameter. Default value is 1.
#' @param k numeric, represents the k-th smallest value from a sample.
#' @param n numeric, represents the size of the sample to compute the order statistic from.
#' @param p numeric, represents the 100p percentile for the distribution of the k-th order statistic. Default value is population median, p = 0.5.
#' @param alpha numeric, (1 - alpha) represents the confidence of an interval for the population percentile p of the distribution of the k-th order statistic. Default value is 0.05.
#' @param ... represents others parameters of the G distribution.
#' @return A list with a random sample of order statistics from a Kumaraswamy G Distribution, the value of its join probability density function evaluated in the random sample
#' and an approximate (1 - alpha) confidence interval for the population percentile p of the distribution of the k-th order statistic.
#' @references Gentle, J, Computational Statistics, First Edition. Springer - Verlag, 2009.
#' @references Naradajah, S. and Rocha, R. (2016) Newdistns: An R Package for New Families of Distributions, Journal of Statistical Software.
#' @author Carlos Alberto Cardozo Delgado <cardozorpackages@gmail.com>.
#' @examples
#' library(orders)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Kumaraswamy Exponential Distribution
#' order_kumg(10,"exp",1,1,k=3,n=50,p=0.5,alpha=0.02)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Kumaraswamy Normal Distribution
#' order_kumg(10,"norm",1,1,k=3,n=50,p=0.5)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Kumaraswamy Log-normal Distribution
#' order_kumg(10,"lnorm",1,1,k=3,n=50,p=0.5)
#' # A sample of size 10 of the 3-th order statistics from
#' # a Kumaraswamy Chi-square Distribution
#' order_kumg(10,"chisq",1,1,k=3,n=50,p=0.5,df=3)
#' @importFrom Newdistns qkumg dkumg
#' @importFrom stats rbeta
#' @export order_kumg

order_kumg <- function(size,spec,a,b,k,n,p=0.5,alpha=0.05,...){
  sample  <- qkumg(initial_order(size,k,n),spec,a,b,...)
  pdf     <- factorial(size)*cumprod(dkumg(sample,spec,a,b,...))[size]
  log_pdf     <- sum(log(2:size)) + sum(log(dkumg(sample,spec,a,b,...)))
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
