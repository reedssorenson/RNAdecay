
#' sigma^2 estimation
#'
#' Calculates the variance (sigma^2) estimate from the sum of the squared errors
#'     from the fit model.
#'
<<<<<<< HEAD:R/fit_var.R
#' @param sse is the minimum SSE (sum of the squared errors) from the slsqp() fitting
#' @param n is the total number of observations of a single gene (e.g., 8 time points X 4 replicates X 4 treatments/genotypes = 128)
=======
#' @param x is the minimum SSE (sum of the squared errors) from the slsqp()
#'     fitting
#' @param n is the total number of observations of a single gene (e.g., 8 time
#'     points X 4 replicates X 4 treatments/genotypes = 128)
>>>>>>> 6348625a4176d22804ebafedf857ed5938765d20:R/sig.fun.R
#'
#' @return returns the sigma squared estimate
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#' sig.fun(1,128)


fit_var <- function(sse, n) {
  sse / n
}
