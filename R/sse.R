
#' sum of the squared errors for null models
#'
#' For a model that uses a constant decay rate or a decaying decay rate, calculates the sum of the squared errors (differences between the supplied data points and the modeled values based on alpha and/or beta values). For these models all treatments are assumed to have the same \code{a} (alpha) and/or \code{b} (beta).
#'
#' @param a alpha value
#' @param b beta value
#' @param m mRNA abundance values
#' @param t time points of \code{m}
#'
#' @return  Returns the sum of the squared errors
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
<<<<<<< HEAD
#' sse_null_decaying_decay(a=0.05, b = 0.001,
#'            m = c(1,1,1,0.99,0.5,0.5,0.5,0.49,0.25,0.25,0.25,0.24,0.12,0.125,0.125,0.126),
#'            t = rep(c(0,10,20,30),each = 4))
#' sse_null_const_decay(a=0.05,
#'            m = c(1,1,1,0.99,0.5,0.5,0.5,0.49,0.25,0.25,0.25,0.24,0.12,0.125,0.125,0.126),
=======
#' sse.nulldExp.fun(a=0.05, b = 0.001,
#'            m = c(1,1,1,0.99,0.5,0.5,0.5,0.49,0.25,0.25,0.25,0.24,0.12,0.125,
#'            0.125,0.126),
#'            t = rep(c(0,10,20,30),each = 4))
#' sse.nullExp.fun(a=0.05,
#'            m = c(1,1,1,0.99,0.5,0.5,0.5,0.49,0.25,0.25,0.25,0.24,0.12,0.125,
#'            0.125,0.126),
>>>>>>> dc1765d3bcd0efd62ae2886244fbbd464064a245
#'            t = rep(c(0,10,20,30),each = 4))


sse_null_decaying_decay <- function(a, b, m, t) {
  sum((m - exp(-(a / b) * (1 - exp(
    -b * t
  )))) ^ 2)
}

#' @rdname sse_null_decaying_decay
#' @export
sse_null_const_decay <- function(a, m, t) {
  sum((m - exp(-a * t)) ^ 2)
}
