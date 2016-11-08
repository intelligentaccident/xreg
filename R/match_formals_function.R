#' Match arguments sent to current function to arguments found in separtae function, and return as list.
#'
#' @param fun Function whose arguments should be matched
#' @param ... named arguments that should override or supplement those provided to current function.
#' @description Match arguments sent to current function to arguments found in separtae function, and return as list.
#' @examples
#' fun1 <- function(a, b) print(paste(a, b, sep = ", "))
#' fun2 <- function(a, c = 2, d) do.call(fun1, match_formals(fun1, b = c))
#' fun2(a = "a", d = "d")

match_formals <- function(fun = as.character(do.call("match.call", args = list(), envir = parent.frame())[1]), ...) {
  f <- formals(fun)
  mcl <- as.list(do.call("match.call", args = list(), envir = parent.frame()))[-1]
  mcln <- names(mcl)
  fn <- names(f)
  mclnx <- mcln[mcln %in% fn]

  f[mclnx] <- mcl[mclnx]

  argl <- list(...)
  argln <- names(argl)
  arglnx <- argln[argln %in% fn]
  f[arglnx] <- argl[arglnx]

  return(f)
}
