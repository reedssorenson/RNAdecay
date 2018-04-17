
#' sigma^2 estimation
#'
#' Calculates the variance (sigma^2) estimate from the sum of the squared errors from the fit model.
#'
#' @param sse is the minimum SSE (sum of the squared errors) from the slsqp() fitting
#' @param n is the total number of observations of a single gene (e.g., 8 time points X 4 replicates X 4 treatments/genotypes = 128)
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
