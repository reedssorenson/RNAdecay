<<<<<<< HEAD

#' Combinatorial groups matrix generator
#'
#' Generates a combinatorial grouping matrix based on the \code{decaydata} data.frame.
#'
#' The resulting matrix of indices is used to constrain treatment alphas or treatment betas in combination. For example, in one model, treatment alphas might be allowed to vary independently (gp1), but the beta models might be constrained to be equal for some treatments indicated by haveing the same index number (other gp).
#'
#' @param decaydata a data.frame with column names:
#'     'geneID','treatment','t.decay','rep','value' with classes
#'     \code{factor}, \code{factor}, \code{numeric}, \code{factor}, \code{numeric}
#'
#' @return returns a matrix of equivalence group indicies based on the number of levels in the 'treatment' column (max of 4).
#'
#' @export
#'
#' @examples
#' groupings(data.frame(geneID=paste0('gene',1:4),treatment=as.factor(paste0('treat',1:4)),
#'                      t.decay=0:3,rep=rep('rep1'),value=c(1,0.5,0.25,0.12)))
#'
groupings <- function(decaydata) {
  nTreat <- length(as.character(unique(decaydata$treatment)))
  nEquivGrp <- if (nTreat == 2)
    2
  else if (nTreat == 3)
    5
  else if (nTreat == 4)
    15
  groupings <- matrix(nrow = nEquivGrp + 1, ncol = nTreat)
  colnames(groupings) <- as.character(unique(decaydata$treatment))
  rownames(groupings) <- paste0("gp", 1:(nEquivGrp + 1))
  if (nTreat == 2) {
    groupings[1,] <- c(1, 2)
    groupings[2,] <- c(1, 1)
    groupings[3,] <- c(NA, NA)
  }
  if (nTreat == 3) {
    groupings[1,] <- c(1, 2, 3)
    groupings[2,] <- c(1, 2, 2)
    groupings[3,] <- c(1, 2, 1)
    groupings[4,] <- c(1, 1, 2)
    groupings[5,] <- c(1, 1, 1)
    groupings[6,] <- c(NA, NA, NA)
  }
  if (nTreat == 4) {
    groupings[1,] <- c(1, 2, 3, 4)
    groupings[2,] <- c(1, 2, 2, 2)
    groupings[3,] <- c(2, 1, 2, 2)
    groupings[4,] <- c(2, 2, 1, 2)
    groupings[5,] <- c(2, 2, 2, 1)
    groupings[6,] <- c(1, 2, 3, 3)
    groupings[7,] <- c(1, 3, 2, 3)
    groupings[8,] <- c(1, 3, 3, 2)
    groupings[9,] <- c(3, 1, 2, 3)
    groupings[10,] <- c(3, 1, 3, 2)
    groupings[11,] <- c(3, 3, 1, 2)
    groupings[12,] <- c(1, 1, 2, 2)
    groupings[13,] <- c(1, 2, 1, 2)
    groupings[14,] <- c(1, 2, 2, 1)
    groupings[15,] <- c(1, 1, 1, 1)
    groupings[16,] <- c(NA, NA, NA, NA)
  }
  return(groupings)
}
=======

#' Combinatorial groups matrix generator
#'
#' Generates a combinatorial grouping matrix based on the \code{decaydata}
#'     data.frame.
#'
#' The resulting matrix of indices is used to constrain treatment alphas or
#'     treatment betas in combination. For example, in one model, treatment
#'    alphas might be allowed to vary independently (gp1), but the beta models
#'    might be constrained to be equal for some treatments indicated by haveing the same index number (other gp).
#'
#' @param decaydata a data.frame with column names:
#'     'geneID','treatment','t.decay','rep','value' with classes
#'     \code{factor}, \code{factor}, \code{numeric}, \code{factor},
#'     \code{numeric}
#'
#' @return returns a matrix of equivalence group indicies based on the number of
#'     levels in the 'treatment' column (max of 4).
#'
#' @export
#'
#' @examples
#' groupings(data.frame(geneID=paste0('gene',seq_len(4)),
#'                      treatment=as.factor(paste0('treat',seq_len(4))),
#'                      t.decay=0:3,rep=rep('rep1'),value=c(1,0.5,0.25,0.12)))
#'
groupings <- function(decaydata) {
  nTreat <- length(as.character(unique(decaydata$treatment)))
  nEquivGrp <- if (nTreat == 2)
    2
  else if (nTreat == 3)
    5
  else if (nTreat == 4)
    15
<<<<<<< HEAD
  groupings <- matrix(nrow = nEquivGrp + 1, ncol = nTreat)
  colnames(groupings) <- as.character(unique(decaydata$treatment))
  rownames(groupings) <- paste0("gp", 1:(nEquivGrp + 1))
=======
  groupings = matrix(nrow = nEquivGrp + 1, ncol = nTreat)
  colnames(groupings) = as.character(unique(decaydata$treatment))
  rownames(groupings) = paste0("gp", seq_len((nEquivGrp + 1)))
>>>>>>> dc1765d3bcd0efd62ae2886244fbbd464064a245
  if (nTreat == 2) {
    groupings[1,] <- c(1, 2)
    groupings[2,] <- c(1, 1)
    groupings[3,] <- c(NA, NA)
  }
  if (nTreat == 3) {
    groupings[1,] <- c(1, 2, 3)
    groupings[2,] <- c(1, 2, 2)
    groupings[3,] <- c(1, 2, 1)
    groupings[4,] <- c(1, 1, 2)
    groupings[5,] <- c(1, 1, 1)
    groupings[6,] <- c(NA, NA, NA)
  }
  if (nTreat == 4) {
    groupings[1,] <- c(1, 2, 3, 4)
    groupings[2,] <- c(1, 2, 2, 2)
    groupings[3,] <- c(2, 1, 2, 2)
    groupings[4,] <- c(2, 2, 1, 2)
    groupings[5,] <- c(2, 2, 2, 1)
    groupings[6,] <- c(1, 2, 3, 3)
    groupings[7,] <- c(1, 3, 2, 3)
    groupings[8,] <- c(1, 3, 3, 2)
    groupings[9,] <- c(3, 1, 2, 3)
    groupings[10,] <- c(3, 1, 3, 2)
    groupings[11,] <- c(3, 3, 1, 2)
    groupings[12,] <- c(1, 1, 2, 2)
    groupings[13,] <- c(1, 2, 1, 2)
    groupings[14,] <- c(1, 2, 2, 1)
    groupings[15,] <- c(1, 1, 1, 1)
    groupings[16,] <- c(NA, NA, NA, NA)
  }
  return(groupings)
}
>>>>>>> 6348625a4176d22804ebafedf857ed5938765d20
