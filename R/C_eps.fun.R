
#' combined adjusted coefficient of variation
#'
#' Calculates the sum of the column standard devation divided by the sum of the column mean and a small value to avoid dividing by 0 (\code{eps})
#'
#' @param X data.frame or matrix of numeric data
#' @param eps small value to add to the mean to avoid dividing by 0; defaults to 1e-4
#'
#' @export
#'
#' @return returns the sum of the coefficients of variation for all columns of \code{X}
#'
#' @examples
#' C_eps.fun(data.frame(test1=rep(0,5),test2=c(0.2,0.3,0.35,0.27,0.21)))


C_eps.fun <- function(X, eps = 1e-04) {
  Y <- as.matrix(X)
  sd.Y <- apply(Y, 2, stats::sd)
  mean.Y <- apply(Y, 2, mean)
  return(sum(sd.Y / (mean.Y + eps)))
}

