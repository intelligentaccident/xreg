% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/match_formals_function.R
\name{match_formals}
\alias{match_formals}
\title{Match arguments sent to current function to arguments found in separtae function, and return as list.}
\usage{
match_formals(
  fun = as.character(do.call("match.call", args = list(), envir = parent.frame())[1]),
  ...
)
}
\arguments{
\item{fun}{Function whose arguments should be matched}

\item{...}{named arguments that should override or supplement those provided to current function.}
}
\description{
Match arguments sent to current function to arguments found in separtae function, and return as list.
}
\examples{
fun1 <- function(a, b) print(paste(a, b, sep = ", "))
fun2 <- function(a, c = 2, d) do.call(fun1, match_formals(fun1, b = c))
fun2(a = "a", d = "d")
}
