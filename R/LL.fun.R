
#' log likelihood
#'
#' Calculates the log likelihood value from the sum of the squared errors, sigma2, and the total number of data points.
#'
#' @param x is the minimum SSE from the slsqp() fitting
#' @param y is sigma2 cooresponding to the minimum SSE
#' @param n is the total number of observations of a single gene (e.g., 8 time points X 4 replicates X 4 treatments/genotypes = 128)
#'
#' @export
#'
#' @return returns Log Likelihood
#'
#' @examples
#' LL.fun(1,1/128,128)

LL.fun <- function(x,y,n){
  -(n/2)*log(2*pi*y) - (1/(2*y))*x
}
