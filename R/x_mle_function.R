#' Modification of the mle function from the stats4 package
#'
#' @param minuslogl Function to be minimized
#' @param start list or named vector of start values for parametres to be fitted
#' @param method Method to be used by \code{\link[stats]{optim}}.
#' @param fixed list of fixed named parameters (will not be passed to \code{\link[stats]{optim}} for fitting).
#' @param nobs number of observations. Not used presently
#' @param solve_hessian If false, will not try to calculate standard errors. Default is TRUE
#' @param return_first If set to TRUE, will return result of minusloglik using values supplied to start
#' @param lower named list of lower bounds for fitting of values. Requires method L-BFGS-B
#' @param upper named list of upper bounds for fitting of values. Requires method L_BFGS-B
#' @param ... Arguments to be passed on to minuslogl or optim
#' @description Modified version of the mle function from stats4. Allows arbitrary parameters in function provided to parmaeter minuslogl. Intended for internal use by xreg.
#' @author Original funciton by R Core Team and contributors. Edits to work properly with xreg by Kim Rand
#' @export
x_mle <- function (minuslogl, start = formals(minuslogl), method = "L-BFGS-B",
                   fixed = list(), nobs, solve_hessian = T, return_first = FALSE, lower = -Inf, upper = Inf, ...) {
  
  call <- match.call()
  optim_string <- "optim"
  if(exists('optimParallel')) if(!is.null(getDefaultCluster())) optim_string <- "optimParallel"
  n <- names(fixed)
  fullcoef <- formals(minuslogl)
  #if (any(!n %in% names(fullcoef)))
  #  stop("some named arguments in 'fixed' are not arguments to the supplied log-likelihood")
  fullcoef[n] <- fixed
  if (!missing(start) && (!is.list(start) || is.null(names(start))))
    stop("'start' must be a named list")
  start[n] <- NULL
  start <- sapply(start, eval.parent)
  nm <- names(start)
  oo <- match(nm, names(fullcoef))
  #if (anyNA(oo))
  #  stop("some named arguments in 'start' are not arguments to the supplied log-likelihood")
  start <- start[order(oo)]
  nm <- names(start)
  f <- function(p, ...) {
    l <- as.list(p)
    names(l) <- nm
    l[n] <- fixed
    l <- c(l, list(...))
    # print(names(l))
    do.call("minuslogl", l)
  }
  oout <- if (length(start)) {
    if(method == "L-BFGS-B") {
      lowerd <- upperd <- start
      lowerd[] <- -Inf
      upperd[] <- Inf
      
      if(is.null(names(lower))) lower <- rep(lower, length(start))
      else lowerd[names(lower)[names(lower) %in% names(lowerd)]] <- lower[names(lower)[names(lower) %in% names(lowerd)]]
      if(is.null(names(upper))) upper <- rep(upper, length(start))
      else upperd[names(upper)[names(upper) %in% names(upperd)]] <- upper[names(upper)[names(upper) %in% names(upperd)]]
      
      
      
      do.call(optim_string, c(list(par = start, fn = f, method = method, hessian = solve_hessian, lower = lowerd, upper = upperd), list(...)))
    } else {
      do.call(optim_string, c(list(par = start, fn = f, method = method, hessian = solve_hessian), list(...)))
    }
    
    
  }
  
  else {
    
    list(par = numeric(), value = f(start, ...))
    
  }
  
  
  coef <- oout$par
  vcov <- if (length(coef) & solve_hessian)
    solve(oout$hessian)
  else matrix(numeric(), 0L, 0L)
  min <- oout$value
  fullcoef[nm] <- coef
  
  fullcoef <- data.frame(Estimate = unlist(fullcoef[nm]))
  
  
  
  if(NROW(fullcoef)){
    fullcoef[, 'Std. Error'] <- NA
    if(solve_hessian) fullcoef[, 'Std. Error'] <-   sqrt(diag(vcov))
  }
  
  
  retval <- list( call = call, coef = coef, fullcoef = fullcoef,
                  vcov = vcov, min = min, details = oout, minuslogl = minuslogl, fixed_values = unlist(fixed),
                  nobs = if (missing(nobs))
                    NA_integer_
                  else nobs, method = method)
  class(retval) <- c("x_mle", "list")
  return(retval)
}