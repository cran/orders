#' @importFrom stats qbinom rbeta pbinom

point_percentile_est <- function(p,size,sample){
  sort_sample <- sort(sample)
  point_est <- sort_sample[floor((size+1)*p)]
  return(point_est)
}

interval_percentile_est <- function(p,size,sample,alpha){
  interval <- sort(sample)[qbinom(c(alpha/2,1-(alpha/2)),size,prob=p)]
  approx_p <- pbinom(qbinom(1-(alpha/2),size,prob=p),size,p) - pbinom(qbinom(alpha/2,size,prob=p),size,p)
  return(c(interval,approx_p))
}

initial_order <- function(size,k,n){return(rbeta(size, k, n + 1 - k))}
