#' Flexible maximum-likelihood-based regression function
#' @param formula a formula or list of formulas to be passed to xregControl and xreg. If a single formula is provided, the left side argument must match the dependent variable. If interval regression is intended, the upper and lower bounds should be assigned to two columns starting with a common name, e.g. value.lb and value.ub. For point estimates, these should be identical.
#' @param df data.frame containing both the continuous and dichotomous data. Will be split into two data.frames using a column matching the name provided in the datatype argument.
#' @param datatype character defining the column name separating dichotomous from continuous data. Continuous data should be interpretable as TRUE in as.logical, and dichotomous as FALSE
#' @param init Optional named list or vector of numeric start values, to override those provided in the xregControl obects. Previously estimated xreg-objects can be used, in which case fitted values will be extracted.
#' @param contdist character option to determine if the continuous model should use normal or logistic distribution. Options are "normal", and "logistic". "normal" instructs xregControl to provide cont_normal(), logistic provides cont_logistic(). Default is "normal".
#' @param dichdist character option to determine if the dichotomous model should use normal or logistic distribution. Options are "normal", and "logistic". "normal" instructs xregControl to provide dich_normal(), logistic provides dich_logistic(). Default is "logistic".
#' @param hetcont optional formula for fitting hetreoscedastic standard deviations for the continuous model. Parameters sharing names with the ones provided in the formula argument will be the same. If this is not intended, different parmaeter names are required.
#' @param hetdich optional formula for fitting hetreoscedastic standard deviations for the dichotomous model. Requires dichotomous distributionn to be "normal". Parameters sharing names with the ones provided in the formula argument will be the same. If this is not intended, different parmaeter names are required.
#' @param ul optional upper bound value for the dependent variable. Values over ul will be set to an interval between the observed and Inf, making the fitting similar to tobit regression. Defaults to Inf.
#' @param ll optional lower bound value for the dependent variable. Values below ll will be set to an interval between the observed and -Inf, making the fitting similar to tobit regression. Defaults to -Inf.
#' @param lntheta If TRUE, parameters or formulas for theta will be exponentiated using exp() prior to fitting, to match the default in the STATA hyreg command. Defaults to TRUE.
#' @param lnsigma If TRUE, parameters or formulas for sigma will be exponentiated using exp() prior to fitting, to match the default in the STATA hyreg command. Defaults to TRUE.
#' @param dichformula optional formula argument for the dichotomous dataset. If not provided, the formula for dichotomous data will be the same as for continuous.
#' @param ... optional arguments to be forwarded to xregControl, xreg, x_mle, or optim
#' @description Wrapper for xregControl and xreg to match the functionality of the STATA hyreg function.
#' @author Kim Rand
#' @examples
#' control_continuous <- xregControl(formulas = list(x ~ y * YVAR + z + ZVAR, value ~ INTERCEPT + x * XVAR), start_values = c(INTERCEPT = 0, XVAR = 1, YVAR = 1, ZVAR = 1), p_fun = cont_normal, name = "CONTINUOUS")
#' control_dichotomos <- xregControl(formulas = list(value ~ INTERCEPT2 + z * ZVAR), start_values = c(INTERCEPT2 = 0, ZVAR = 1), p_fun = dich_logistic, name = "DICHOTOMOUS")
#' joint_control <- c(control_continuous, control_dichotomous)
#' xreg(controlList = joint_control, dataList = list(CONTINUOUS = df1, DICHOTOMUS = df2))
#' @export
hyreg <- function(formula, df, datatype = "method", init = numeric(), contdist = "normal", dichdist = "logistic", hetcont = NULL, hetdich = NULL, ul = Inf, ll = -Inf, lntheta = T, lnsigma = T, dichformula = NULL, ...) {
  if(!"data.frame" %in% class(df)) stop(paste0(class(df)[1], " provided to the df argument. Should be data.frame."))
  if(!datatype %in% colnames(df)) stop(paste0("No column named \"", datatype, "\" in provided dataframe."))
  if(!all(as.logical(unique(df[,datatype])) %in% c(TRUE, FALSE))) stop("Datatype column in provided data.frame must be interpretable as boolean.")
  if(!contdist %in% c("normal", "logistic")) stop(paste0("contdist-argument has to be either \"normal\" or \"logistic\". "))
  if(!dichdist %in% c("normal", "logistic")) stop(paste0("dichdist-argument has to be either \"normal\" or \"logistic\". "))
  if(is.null(dichformula)) dichformula <- formula
  if(class(formula) == "formula") formula <- list(formula)
  if(class(dichformula) == "formula") dichformula <- list(dichformula)
  dataList <- list(Continuous = (cont_df <- df[as.logical(df[, datatype] == TRUE),]), Dichotomous = (dich_df <- df[as.logical(df[, datatype] == FALSE),]))

  dot_args <- list(...)
  dot_args[c("p_fun", "p_aggrgation_fun")] <- NULL


  if(!is.null(hetcont)) {
    if(class(hetcont) == "character") {
      hetcont_tmp <- hetcont
      tmp <- try0(as.formula(hetcont_tmp))
      if(is.null(tmp)){
        hetcont_tmp <- paste0("~ ", hetcont_tmp)
        tmp <- try0(as.formula(hetcont_tmp))
      }
      if(is.null(tmp)) stop("Coult not generate useful formula from provided hetcont-argument.")
      hetcont <- tmp
    }
    hetcont <- change_target(hetcont, "sigma_est")
  } else {
    hetcont <- as.formula("sigma_est ~ SIGMA")
  }

  if(!is.null(hetdich)) {
    if(class(hetdich) == "character") {
      hetdich_tmp <- hetdich
      tmp <- try0(as.formula(hetdich_tmp))
      if(is.null(tmp)){
        hetdich_tmp <- paste0("~ ", hetdich_tmp)
        tmp <- try0(as.formula(hetdich_tmp))
      }
      if(is.null(tmp)) stop("Coult not generate useful formula from provided hetdich-argument.")
      hetdich <- tmp
    }
    hetdich <- change_target(hetdich, "theta_est")
  } else {
    hetdich <- as.formula("theta_est ~ THETA")
  }

  if(lnsigma) hetcont <- appendmodel(hetcont, " exp(", ")")
  if(lntheta) hetdich <- appendmodel(hetdich, " exp(", ")")


  cont_formula <- c(formula, hetcont)
  dich_formula <- c(dichformula, hetdich)


  cont_p_fun <- cont_normal
  dich_p_fun <- dich_logistic
  if(contdist == "logistic") cont_p_fun <- cont_logistic
  if(dichdist == "normal") dich_p_fun <- dich_normal
  cont_control <- xregControl(formulas = cont_formula, p_fun = cont_p_fun, name = "Continuous", censor_bounds = c(ll, ul))
  dich_control <- xregControl(formulas = dich_formula, p_fun = dich_p_fun, name = "Dichotomous")
  control <- c(cont_control, dich_control)

  #if(!is.na(retn)) return(get(retn))

  ret <- do.call("xreg", c(list(
    controlList = control,
    dataList = dataList,
    start_values = init,
    run_from = "hyreg"),
    dot_args))


  return(ret)
}
