
#' calculates bounds for modeled parameters
#'
#' Calculates maximum and minimum bounds for parameter alpha based on experimental time points (t_0, t_1, t_2, t_3, ..., t_max). If RNA level is too low at t_1, then the decay has happened before our observations began - there is an upper bound to the decay rate we can detect (a.high).
#' If RNA level is too high at t_max, then relatively little decay has happened and we can not distinguish the decay rate and the decay of the decay rate  - there is a lower bound to the base decay rate of the decaying decay model (a.low).
#'
#' Similarly, limits on beta are required to prevent precude ranges in which the decay rate and decaing decay are indistinguishable. See vignette 2: _RNA Decay Modeling_ for more information.
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

b.low <- function(t_max){
  10^(floor(log((1/t_max),base=10)))
}




