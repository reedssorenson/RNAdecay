
#' calculates bounds for modeled parameters
#'
#' calculates maximum and minimum bounds for parameters alpha and beta based on experimental time points. These are based on the time points in each experiment. For example if an RNA degrades to 0.1% of the initial (T0) amount by 10 min but the 1st time point measured is at 2 h, it will not be possible to estimate that decay rate.
#'
#' @param t_min time of first experiemtal time point after inhibition of transcription (not T0)
#' @param t_max time of last experimental time point
#'
#' @export
#'
#' @return returns the lowest/highest parameter values to be used as bounds on modeled parameters
#'
#' @examples
#' a.high(7.5)
#' a.low(480)
#' b.low(480)

a.high <- function(t_min){
  signif(-(1/t_min)*log(0.005)	,2)
}

#' @rdname a.high
#' @export

a.low <- function(t_max){
  10^(floor(log(-(1/t_max)*log(0.95),base=10)))
}

#' @rdname a.high
#' @export

b.high <- function() {
  0.075
}

#' @rdname a.high
#' @export

b.low <- function(t_max){
  10^(floor(log((1/t_max),base=10)))
}




