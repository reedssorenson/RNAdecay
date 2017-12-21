
#' sigma^2 estimation
#'
#' calculates the variance (sigma^2) estimate from the sum of the squared errors from the fit model.
#'
#' @param x is the minimum SSE (sum of the squared errors) from the slsqp() fitting
#' @param n is the total number of observations 8 time points, 4 replicates, 4 treatments (or genotypes)
#'
#' @export
#'
#' @return returns the sigma squared estimate
#'
#' @examples
#' sig.fun(1,128)

sig.fun <- function(x,n=128){
  x/n
}
