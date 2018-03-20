
#' Akaike information criterion (with correction)
#'
#' Calculates AIC or AICc.
#'
#' @param maxlLik maximum log likelihood value identified upon model convergence
#' @param p number of parameters in the model
#' @param n is the total number of observations of a single gene (e.g., 8 time
#'     points X 4 replicates X 4 treatments/genotypes = 128)
#'
#' @export
#'
#' @return returns the AIC or AICc values
#'
#' @examples
#' fun_aicc(100,5,15)


fun_aic = function(maxlLik, p) {
  2 * p - 2 * maxlLik
}

#' @rdname fun_aic
#' @export
fun_aicc = function(maxlLik, p, n) {
  2 * p - 2 * maxlLik + ((2 * p * (1 + p)) / (n - p - 1))
}
