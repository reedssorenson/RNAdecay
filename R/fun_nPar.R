
#' number of Parameters function
#'
#' Calculates number of parameters for a specified model given the model
#'     parameter constraints.
#'
#' @param model model name (e.g. 'mod1')
#' @param mod two column data.frame with combination of alpha and beta grouping
#'     numbers in each model with rownames 'modX'
#' @param group matrix of all treatment alpha or beta equivalence groups
#'
#' @export
#'
#' @return returns the integer value of number of parameters in \code{model}
#'
#' @examples
#' fun_nPar('mod1',data.frame('a'=1:5,'b'=rep(2,5),row.names=paste0('mod',1:5)),
#'          t(matrix(c(1,2,3,4,1,2,2,2,1,1,2,2,1,1,1,2,1,2,1,2,1,2,2,1),nrow=4)))


fun_nPar = function(model, mod, group) {
  length(unique(group[mod[model, "a"],])) + 1 + if (
    any(is.na(group[mod[model, "b"],]))
    ) {
    0
  } else {
    length(unique(unlist(group[mod[model, "b"],])))
  }
}
