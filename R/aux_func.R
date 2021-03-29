#' @importFrom stats qbinom rbeta

interval_median <- function(size,sample,alpha){
  output <- sort(sample)[qbinom(c(alpha/2,1-alpha/2),size,0.5)]
  return(output)
}

initial_order <- function(size,k,n){return(rbeta(size, k, n + 1 - k))}
