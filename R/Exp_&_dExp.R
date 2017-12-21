
#' exponential decay functions
#'
#' single exponential decay function (case when betas=0) e^-a*t; double exponential decay function  e^(-(a/b)*(1-e^(-b*t))); (normalized so at t=0 the function is 1)
#'
#' @param t time (in minutes)
#' @param a alpha (in per time, thus in per minute when time is in minutes)
#' @param	par vector of length 2 containing alpha (par[1]) and beta (par[2]) values; alpha=initial decay rate, beta=decay of decay rate (both in per time, thus in per minute when time is in minutes)
#'
#' @return returns abundance after time \code{t}  at alpha initial decay rate and beta decay of decay rate relative to an initial abundance of 1
#'
#' @export
#'
#' @examples
#' Exp(10,log(2)/10) ## returns 0.5
#' dExp(10,c(log(2)/10,0.01)) ##returns 0.5170495

Exp <- function(t, a) {
  exp(-a * t)
}


#' @rdname Exp
#' @export
dExp <- function(t, par) {
  a = par[1]
  b = par[2]
  exp(-(a / b) * (1 - exp(-b * t)))
}
