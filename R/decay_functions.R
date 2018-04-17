
#' exponential decay functions
#'
<<<<<<< HEAD:R/decay_functions.R
#' Constant decay rate function (const_decay(), case when betas=0) e^-a*t; decaying decay rate function (decaying_decay())  e^(-(a/b)*(1-e^(-b*t))). Functions are normalized so at t=0 the function is 1.
=======
#' Constant decay rate function (Exp(), case when betas=0) e^-a\*t; decaying
#' decay rate function (dExp())  e^(-(a/b)\*(1-e^(-b\*t))). Functions are
#' normalized so at t=0 the function is 1.
>>>>>>> 6348625a4176d22804ebafedf857ed5938765d20:R/Exp_&_dExp.R
#'
#' @param t time (in minutes)
#' @param a alpha (in per time, thus in per minute when time is in minutes)
#' @param par vector of length 2 containing alpha (par[1]) and beta (par[2])
#'     values; alpha=initial decay rate, beta=decay of decay rate (both in per
#'     time, thus in per minute when time is in minutes)
#'
#' @return returns abundance after time \code{t}  at alpha initial decay rate
#'     and beta decay of decay rate relative to an initial abundance of 1
#'
#' @export
#'
#' @examples
#' const_decay(10,log(2)/10) ## returns 0.5
#' decaying_decay(10,c(log(2)/10,0.01)) ##returns 0.5170495

const_decay <- function(t, a) {
  exp(-a * t)
}


#' @rdname const_decay
#' @export
decaying_decay <- function(t, par) {
  a <- par[1]
  b <- par[2]
  exp(-(a / b) * (1 - exp(-b * t)))
}
