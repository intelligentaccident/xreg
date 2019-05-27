#' Formula, start value, and custom criterion function fitting for xreg#'
#' @param formulas List of formulas to be fitted for data.frames with same name. If more thann one formula is provided, a column named like the the left-hand side will be generated  in the data.frame during calculation, which can then be used in subsequent formulas. See examples
#' @param start_values List or named vector of start values for parmaeters to be fitted
#' @param fixed_values List or named vector of fixed values to be used in formulas.
#' @param p_fun function to generate likelihood for each observation or interval of observations. Some standard functions are provided. Defaults to least-squares using normally distributed maximum likelihood.
#' @param p_aggregatio_fun function to aggregate likelihoods (or similar), to be returned to optim
#' @param weights_var optional parameter. If provided, should match a column in provided data.frame, and will be used as a frequency weight when aggregating likelihoods.
#' @param name Name of control-set. If none is provided, names will be set by the order of controlsets. data.frames will be matched to control-sets using these names (or order).
#' @param censor_bounds Optional observation bounds for interval regression, e.g. c(1, 2). Defaults to c(-Inf, Inf)
#' @param lower Named vector of optional box constraint lower bounds for fitting of parameters. E.g. c(VAR1 = 0) will force VAR1 to be >=0.
#' @param upper Named vector of optional box constraint upper bounds for fitting of parameters. E.g. c(VAR1 = 1) will force VAR1 to be <=1.
#' @description Allows a list of formulas and fitting functions to be used for each type or set of data,
#' @aliases c.xregControl, c.xregControlList
#' @author Kim Rand
#' @examples
#' control_continuous <- xregControl(formulas = list(x ~ y * YVAR + z + ZVAR, value ~ INTERCEPT + x * XVAR), start_values = c(INTERCEPT = 0, XVAR = 1, YVAR = 1, ZVAR = 1), p_fun = cont_normal, name = "CONTINUOUS")
#' control_dichotomos <- xregControl(formulas = list(value ~ INTERCEPT2 + z * ZVAR), start_values = c(INTERCEPT2 = 0, ZVAR = 1), p_fun = dich_logistic, name = "DICHOTOMOUS")
#' joint_control <- c(control_continuous, control_dichotomous)
#' xreg(controlList = joint_control, dataList = list(CONTINUOUS = df1, DICHOTOMUS = df2))
#' @export


xregControl <- function(formulas,
                        start_values = numeric(),
                        fixed_values = numeric(),
                        p_fun = cont_normal,
                        p_aggregation_fun = function(d_df) return(-d_df$p * d_df$internal_count),
                        weights_var = NA,
                        name = NA,
                        censor_bounds = c(-Inf, Inf),
                        lower = NA,
                        upper = NA) {
  if(class(formulas)[1] %in% c("xregControl", "xregControlList")) return(c(formulas))
  if(class(formulas)[1] == "formula") formulas <- list(formulas)
  required <- names(formals(p_fun))


  defined_vars <- character()
  for(formula in formulas) {
    targetName <- as.character(formula[[2]])
    defined_vars <- c(defined_vars, targetName)
  }

  if(length(required)) {
    required <- names(formals(p_fun))
    if(length(required) > 1) {
      req <- as.list(rep(0, length(required)))
      names(req) <- required
      for(i in 2:length(required)) {
        if(!names(req)[i] %in% defined_vars) {
          required <- req
          defined_vars <- c(defined_vars, names(req)[i])
          required[[i]] <- 1
          req_return <- do.call(p_fun, args = required)
          formulas <- c(formulas, list(req_return$formula))
          if(exists("start_values", req_return)) start_values  <- c(start_values, req_return$start_values[!req_return$start_values %in% start_values])
          if(exists("fixed_values", req_return)) fixed_values  <- c(fixed_values, req_return$fixed_values[!req_return$fixed_values %in% fixed_values])
        }
      }
    }
  }
  start_values <- start_values[unique(names(start_values))]


  retObj <- mget(c(names(formals()), "defined_vars"))

  retObj$columnNames <- character()
  for(formula in formulas) {
    retObj$columnNames <- c(retObj$columnNames, all.vars(formula)[!all.vars(formula) %in% retObj$columnNames])
  }

  retObj$columnNames <- retObj$columnNames[!retObj$columnNames %in% names(start_values)]

  class(retObj) <- c("xregControl", "list")
  #print("GOTHERE")
  retObj <- c(retObj)

  return(retObj)
}

c.xregControl <- function(...) {
  args <- list(...)

  if(is.null(names(args) )) names(args) <- NA
  for(i in 1:length(args)) {
    if(is.na(names(args)[i])) {
      if(is.na(args[[i]]$name)) {
        names(args)[i] <- i
      } else {
        names(args)[i] <- args[[i]]$name
      }
    }
    if(!"xregControl" %in% class(args[[i]])) stop("At least one non-xregControl object among arguments.")
  }
  class(args) <- c("xregControlList", "list")


  columnNames <- character()
  startValues <- numeric()
  fixedValues <- numeric()
  lower <- numeric()
  upper <- numeric()
  for(xregControl in args) {
    columnNames <- c(columnNames, xregControl$columnNames[!xregControl$columnNames %in% columnNames])
    startValues <- c(startValues, xregControl$start_values[!names(xregControl$start_values) %in% names(startValues)])
    fixedValues <- c(fixedValues, xregControl$fixed_values[!names(xregControl$fixed_values) %in% names(fixedValues)])
    lower <- c(lower, xregControl$lower[!xregControl$lower %in% lower])
    upper <- c(upper, xregControl$upper[!xregControl$upper %in% upper])
  }


  attr(args, which = "columnNames") <- columnNames
  attr(args, which = "startValues") <- startValues
  attr(args, which = "fixedValues") <- fixedValues
  attr(args, which = "lower") <- lower
  attr(args, which = "upper") <- upper


  return(args)
}


c.xregControlList <- function(...) {  #return(do.call(c.xregControl, do.call(c.list, list(...))))
  args <- list(...)
  for(narg in 1:length(args)) {
    class(args[[narg]]) <- "list"
  }
  #return(args)
  return(do.call(c.xregControl, do.call(c, args)))
}
