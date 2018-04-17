
#' combined adjusted coefficient of variation
#'
#' Calculates the sum of the column standard devation divided by the sum of
#'     the column mean and a small value to avoid dividing by 0 (\code{eps})
#'
#' @param X data.frame or matrix of numeric data
#' @param eps small value to add to the mean to avoid dividing by 0; defaults
#'     to 1e-4
#'
#' @return returns the sum of the coefficients of variation for all columns of \code{X}
#'
#' @export
#'
<<<<<<< HEAD:R/comb_cv.R
#' @keywords internal
=======
#' @return returns the sum of the coefficients of variation for all columns of
#'      \code{X}
>>>>>>> 6348625a4176d22804ebafedf857ed5938765d20:R/C_eps.fun.R
#'
#' @examples
#' comb_cv( data.frame( test1=rep(0,5), test2=c(0.2,0.3,0.35,0.27,0.21) ) )

comb_cv <- function(X, eps = 1e-04) {
  Y <- as.matrix(X)
  sd.Y <- apply(Y, 2, stats::sd)
  mean.Y <- apply(Y, 2, mean)
  return(sum(sd.Y / (mean.Y + eps)))
}

