
#' Akaike information criterion (with correction)
#'
#' calculates AIC or AICc
#'
#' @param maxlLik maximum log likelihood value identified upon model convergence
#' @param p number of parameters in the model
#' @param n number of measured values for a single gene
#'
#' @export
#'
#' @return returns the AIC or AICc values
#'
#' @examples
#' fun_aicc(100,5,15)

fun_aic = function(maxlLik,p) {2 * p - 2 * maxlLik}

#' @rdname fun_aic
#' @export
fun_aicc = function(maxlLik, p, n) {
  2 * p - 2 * maxlLik   +   ((2 * p * (1 + p)) / (n - p - 1))
}
