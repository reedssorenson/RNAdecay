
#' Indexes column names with multiple patterns (i.e., multigrep)
#'
#' Identifies dataframe column names that have all of the pattern arguments
#'     (up to 4).
#'
#' Be aware that column data labels that are part of another data label are
#'     not advisable (e.g. mut1, mut2, mut1.mut2; cols(df,'mut1') will return
#'      indices for both 'mut1' and 'mut1.mut2' labeled columns
#'
#' @param df a dataframe with column names to parse
#' @param w,x,y,z character string or regular expression passed to grep pattern
#'     argument
#'
#' @return  returns a vector of integer indices of the column names of \code{df}
#'     that match to all of \code{w}, \code{x}, \code{y}, and \code{z}
#'
#' @export
#'
#' @examples
#' cols(data.frame(xyz=1:5,zay=6:10,ybz=11:15,tuv=16:20),'y','z') ##
#' # returns 1 2 3

cols = function(df,
                w,
                x = NA,
                y = NA,
                z = NA) {
  if (is.na(x))
    grep(w, names(df))
  else if (is.na(y))
    intersect(grep(w, names(df)), grep(x, names(df)))
  else if (is.na(z))
    intersect(
      intersect(
        grep(w, names(df)), grep(x, names(df))),
      grep(y, names(df)
           )
      )
  else
    intersect(intersect(intersect(grep(w, names(
      df
    )), grep(x, names(
      df
    ))), grep(y,
              names(df))), grep(z, names(df)))
}

