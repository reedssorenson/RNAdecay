
#' contraint function list maker
#'
#' Individual double exponential models are all nested within model number 1 in
#' which alpha and beta parameters vary independently for each treatment. Models
#' that assume no difference in parameters between specific treatments manifest
#' as constraints in the modeling. These constraints are coded as functions
#' that are passed to the optimization process. Each model has a distinct
#' constraint function.
#'
#' @param mods data.frame specifying alpha and beta group pairs for each model
#' @param groups grouping matrix for alphas or betas
#'
#' @return  Returns a list of constraint functions to be passed to the
#' optimization function.
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#' constraint_fun_list_maker(mods = data.frame(a = c(1,1,1,2,2,2), b = c(1,2,3,1,2,3),
#'                           row.names = paste0('mod',1:6)),
#'                           groups = data.frame(treat1 = c(1,1,NA), treat2 = c(2,1,NA)))
#'

constraint_fun_list_maker <- function(mods, groups) {
  # function to make a function with an extra argument 'body.subset' to select lines
  # of the function body by numerical index
  make_function <-
    function(body.subset,
             args,
             body,
             last.line,
             env = parent.frame()) {
      f <- function() {

      }
      formals(f) <- args
      body(f) <- body[c(1:2, body.subset + 2, last.line - 1)]
      environment(f) <- env
      f
    }

  if (ncol(groups) == 4) {
    const_Ind = list(
      c(NULL),
      c(4, 5),
      c(2, 3),
      c(1, 3),
      c(1, 2),
      6,
      5,
      4,
      3,
      2,
      1,
      c(1, 6),
      c(2, 5),
      c(3, 4),
      c(1, 2, 3),
      NULL
    )
    names(const_Ind) <- paste0("gp", 1:16)
    const_Ind <- unlist(lapply(const_Ind[1:15], function(x)
      lapply(lapply(const_Ind,
                    function(y)
                      y + 6), function(z)
                        c(x, z))), recursive = FALSE)
    names(const_Ind) <- paste0("mod", 1:240)
    constraint_fun_list <- lapply(const_Ind, function(x) {
      make_function(
        body.subset = x,
        args = alist(pars =),
        body = as.call(c(
          as.name("{"),
          expression(
            g <-
              numeric(),
            g[length(g) + 1] <- pars[1] - pars[2],
            g[length(g) + 1] <- pars[1] - pars[3],
            g[length(g) + 1] <- pars[1] - pars[4],
            g[length(g) + 1] <- pars[2] - pars[3],
            g[length(g) + 1] <- pars[2] - pars[4],
            g[length(g) + 1] <- pars[3] - pars[4],
            g[length(g) + 1] <- pars[5] - pars[6],
            g[length(g) + 1] <- pars[5] - pars[7],
            g[length(g) + 1] <- pars[5] - pars[8],
            g[length(g) + 1] <- pars[6] - pars[7],
            g[length(g) + 1] <- pars[6] - pars[8],
            g[length(g) + 1] <- pars[7] - pars[8],
            return(g)
          )
        )),
        last.line = 16
      )
    })
  }

  if (ncol(groups) == 3) {
    const_Ind <- list(c(NULL), 3, 2, 1, 1:2, NULL)
    names(const_Ind) <- paste0("gp", 1:6)
    const_Ind <- unlist(lapply(const_Ind[1:5], function(x)
      lapply(lapply(const_Ind,
                    function(y)
                      y + 3), function(z)
                        c(x, z))), recursive = FALSE)
    names(const_Ind) <- paste0("mod", 1:30)
    constraint_fun_list <- lapply(const_Ind, function(x) {
      make_function(
        body.subset = x,
        args = alist(pars =),
        body = as.call(c(
          as.name("{"),
          expression(g <- numeric(),
            g[length(g) + 1] <- pars[1] - pars[2],
            g[length(g) + 1] <- pars[1] - pars[3],
            g[length(g) + 1] <- pars[2] - pars[3],
            g[length(g) + 1] <- pars[4] - pars[5],
            g[length(g) + 1] <- pars[4] - pars[6],
            g[length(g) + 1] <- pars[5] - pars[6],
            return(g)
          )
        )),
        last.line = 10
      )
    })
  }

  if (ncol(groups) == 2) {
    const_Ind <- list(c(NULL), 1, NULL)
    names(const_Ind) <- paste0("gp", 1:3)
    const_Ind <- unlist(lapply(const_Ind[1:2], function(x)
      lapply(lapply(const_Ind,
                    function(y)
                      y + 1), function(z)
                        c(x, z))), recursive = FALSE)
    names(const_Ind) <- paste0("mod", 1:6)
    constraint_fun_list <- lapply(const_Ind, function(x) {
      make_function(
        body.subset = x,
        args = alist(pars =),
        body = as.call(c(
          as.name("{"),
          expression(g <- numeric(),
                     g[length(g) + 1] <- pars[1] - pars[2],
                     g[length(g) + 1] <- pars[3] - pars[4],
                     return(g)
                     )
        )),
        last.line = 6
      )
    })


  }

  if (ncol(groups) == 1) {
    const_Ind <- list(NULL)
    names(const_Ind) <- paste0("gp", 1)
    const_Ind <- list(c(NULL), NULL)
    names(const_Ind) <- paste0("gp", 1)
    const_Ind <- unlist(lapply(const_Ind[1], function(x)
      lapply(lapply(const_Ind,
                    function(y)
                      y), function(z)
                        c(x, z))), recursive = FALSE)
    names(const_Ind) <- paste0("mod", 1:2)
    constraint_fun_list <- lapply(const_Ind, function(x) {
      make_function(
        body.subset = x,
        args = alist(pars =),
        body = as.call(c(
          as.name("{"),
          expression(g <- numeric(),
                     return(g)
                     )
        )),
        last.line = 4
      )
    })
  }
  return(constraint_fun_list)
}
