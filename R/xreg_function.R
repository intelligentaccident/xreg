#' Flexible maximum-likelihood-based regression function
#' @param controlList one or more xregControl-objects (concatenated using c() ) to control maximum-likelihood estimation for different types or sets of data provided in dataList. Can alternatively provide a single formula or list of formulas, in which case an xregControl object will be generated on the fly. An xreg-object can be provided, in which case the xregControlList will be extracted. See examples
#' @param dataList Named list of data.frames for which maximum sum-likelihood should be estimated. Names must correspond to names of xregControl objects. If no names are provided, order will be used to decide which xregControl-object to use.
#' @param start_values Optional named list or vector of numeric start values, to override those provided in the xregControl obects. Previously estimated xreg-objects can be used, in which case fitted values will be extracted.
#' @param fixed_values Optional named list or vector of numeric values that will be held as fixed, and will be available in the same way as the fitted parameters. An xreg-object can be used.
#' @param latent_classes Optional integer parameter defining the number of latent classes to be fitted
#' @param latent_class_parameter Optional character vector with names corresponding to the parameters that are to be varied in the latent classes (if latent_classes > 1)
#' @param latent_id_colname Optional parameter defining the name of the id column defining the individuals (or groups) over which latent classes are to be fitted
#' @param return_type Optional debug function. Default value, "fit", returns and xreg-object with fitted values. "df" returns the final data.frame of estimated values. "first" returns likelihood for providedd start values. "first_df" returns the initial data.frame with start_values. "predict" returns final data.frame. "precalc_df" returns data.frame prior to calculation of likelihood. Useful for debugging of likelihood-function.
#' @param print_sum Optional argument to force continuous printing of the progress of fitting. Can be useful for debugging.
#' @param ... optional arguments to be forwarded to x_mle, optim, xregControl, etc.
#' @description Allows a list of formulas and fitting functions to be used for each type or set of data.
#' @author Kim Rand-Hendriksen
#' @examples
#' control_continuous <- xregControl(formulas = list(x ~ y * YVAR + z + ZVAR, value ~ INTERCEPT + x * XVAR), start_values = c(INTERCEPT = 0, XVAR = 1, YVAR = 1, ZVAR = 1), p_fun = cont_normal, name = "CONTINUOUS")
#' control_dichotomos <- xregControl(formulas = list(value ~ INTERCEPT2 + z * ZVAR), start_values = c(INTERCEPT2 = 0, ZVAR = 1), p_fun = dich_logistic, name = "DICHOTOMOUS")
#' joint_control <- c(control_continuous, control_dichotomous)
#' xreg(controlList = joint_control, dataList = list(CONTINUOUS = df1, DICHOTOMUS = df2))
xreg <- function(controlList,
                 dataList,
                 start_values = numeric(),
                 fixed_values = numeric(),
                 latent_classes = 0,
                 latent_class_parameters = character(),
                 latent_id_colname = character(),
                 return_type = "fit",
                 print_sum = F,
                 ...) {


  dot_args <- list(...)
  dot_args[c("p_fun", "p_aggrgation_fun")] <- NULL

  if(latent_classes < 2) latent_class_parameters <- character()
  if(!is.character(latent_class_parameters)) {
    warning("Latent_class_parameters not character. Latent classes will be disregarded.")
    latent_class_parameters <- character()
    latent_classes <- 0
    latent_id_colname <- character()
  }

  if(!is.character(latent_id_colname)) {
    warning("latent_id_colname not character. Latent classes will be disregarded.")
    latent_class_parameters <- character()
    latent_classes <- 0
    latent_id_colname <- character()
  }
  if(length(latent_id_colname)) latent_id_colname <- latent_id_colname[1]



  xlik <- function(...,
                   controlList,
                   dataList,
                   print_sum = FALSE,
                   return_lik_df = FALSE,
                   return_first_df = FALSE,
                   latent_classes = 0,
                   latent_class_parameters = character()
  ) {
    args <- list(...)



    inner_loop <- function(...,
                           formulas,
                           valueVar,
                           data_df,
                           p_fun,
                           p_aggregation_fun) {

      formula_loop <- function(cur_suffix = "") {

        for(formulanum in 1:length(formulas)) {
          formula <- formulas[[formulanum]]
          new_var <- as.character(formula[[2]])
          if(new_var == valueVar) new_var = "Xb"
          d_df[, new_var] <- with(args, with(d_df, eval(parse(text = formula))))
        }
        d_df$p <- p_fun(d_df)
        d_df$wp <- p_aggregation_fun(d_df)
        colnames(d_df)[match(c("Xb", "p", "wp"), colnames(d_df))] <- paste0(c("Xb", "p", "wp"), cur_suffix)
        return(d_df)
      }


      args <- list(...)
      #print(names(args))
      d_df <- data_df


      if(latent_classes > 0 && length(latent_class_parameters)>0) {
        for(i in 1:latent_classes) {

          args[[latent_class_parameters]] <- args[[paste(latent_class_parameters, i, sep = ".")]]
          d_df <- formula_loop(paste0(".", i))
        }
      } else {
        d_df <- formula_loop()
      }

      return(d_df)

    }


    resList <- list()


    for(dataName in names(dataList)) {
      data_df <- dataList[[dataName]]
      #print(colnames(data_df))
      control <- controlList[[dataName]]
      resList[[dataName]] <- with(control, {
        return(inner_loop(...,
                          valueVar = valueVar,
                          formulas = formulas,
                          data_df = data_df,
                          p_fun = p_fun,
                          p_aggregation_fun = p_aggregation_fun))
      })
    }

    sumList <- list()
    ret_sum <- 0

    if(latent_classes > 0 && length(latent_class_parameters)>0)  {

      for(i in 2:latent_classes) {
        if(any(as.logical(unlist(args[paste0(latent_class_parameters, ".", i)])) < unlist(args[paste0(latent_class_parameters, ".", i-1)]))) return(999999999)
      }


      for(i in 1:length(resList)) {
        if(i == 1) lc_df <- resList[[i]][, c("internal_classgroup", paste0("wp.", 1:latent_classes))]
        else lc_df <- rbind(lc_df, resList[[i]][, c("internal_classgroup", paste0("wp.", 1:latent_classes))])
      }
      latent_sums <- aggregate(lc_df[,2:(latent_classes +1)], by = list(internal_classgroup = lc_df$internal_classgroup), FUN = sum)

      latent_sums$latent_class <- apply(latent_sums[, 2:NCOL(latent_sums)], 1, function(x) match(min(x), x))
      for(i in 1:length(resList)) {
        resList[[i]]$latent_class <- latent_sums$latent_class[match(resList[[i]]$internal_classgroup, latent_sums$internal_classgroup)]
        resList[[i]]$Xb <- resList[[i]][, paste0("Xb.", 1:latent_classes)][tmp <- matrix(c(1:NROW(resList[[i]]), resList[[i]]$latent_class), ncol = 2)]
        resList[[i]]$p <- resList[[i]][, paste0("p.", 1:latent_classes)][tmp]
        resList[[i]]$wp <- resList[[i]][, paste0("wp.", 1:latent_classes)][tmp]
      }
    }


    for(i in 1:length(resList)) ret_sum <- ret_sum + (sumList[[names(resList)[i]]] <- sum(resList[[i]]$wp))

    if(return_first_df) return(resList)

    if(print_sum) {

      this_call <- match.call()
      #print(names(this_call))

      #print(names(args))
      this_call[!names(this_call) %in% names(args)] <- NULL

      print(unlist(as.list(this_call)))

      # print(resList)

      print(paste("Sum: ",ret_sum))

    }

    if(is.na(ret_sum)) ret_sum <- 9999999999999999999

    return(ret_sum)
  }



  lik_fun = xlik
  startValues <- attr(controlList, "startValues")
  if("xreg" %in% class(controlList)) {
    if(is.na(start_values[1])) start_values <- controlList$coef
    controlList <- controlList$controlList
  }
  formulas <- controlList
  #return(match_formals(fun = xregControl, formulas = controlList))
  controlList <- do.call(xregControl, match_formals(fun = xregControl, formulas = controlList))
  #print(controlList)
  #if("p_fun" %in% names(list(...))) list(...)[["p_fun"]] <- NULL
  if("xreg" %in% class(start_values)) startValues <- start_values$coef
  if("numeric" %in% class(start_values) & length(start_values) > 0) startValues <- start_values

  if(class(dataList) == "data.frame") {
    dataList <- list(dataList)
    names(dataList)[1] <- names(controlList)[1]
    warning(paste0("Dataframe provided where list expected in dataList argument. The dataframe has been enclosed in a list, and named \"", names(controlList)[1], "\" to match the first xregControl-object."))
  }

  if(!all(names(dataList) %in% names(controlList))) stop(paste("No control list provided for ", paste(names(dataList)[!names(dataList) %in% names(controlList)], collapse = ", ")))


  if(is.null(names(dataList))) {
    if(length(dataList) <= length(controlList)) {
      warning("Provided datalist was unnamed. Names coerced from controlList")
      names(dataList) <- names(controlList)[1:length(datalist)]
    } else {
      stop("Provided datalist was unnamed, and longer than providedd controllist. Exiting.")
    }
  }
  if(length(unique(names(startValues))) < length(startValues)) {
    warning(paste0("More than one start value assigned for parameters (", paste(names(startValues)[-match(unique(names(startValues)), names(startValues))], collapse = ", "), "). First value used."))
    startValues <- startValues[unique(names(startValues))]
  }

  if(!firstNA(fixed_values)) {
    if("xreg" %in% class(fixed_values)) {
      fixed_values <- fixed_values$pars
      #print("Gothere")
    }
    if("data.frame" %in% class(fixed_values)) {
      fixed_df <- fixed_values[, 1:2]
      fixed_values <- fixed_values[,1]
      names(fixed_values) <- rownames(fixed_df)
    } else {
      #print(fixed_values)
      fixed_df <- data.frame(Estimate = fixed_values, `Std. Error` = rep(NA, NROW(fixed_values)))
      rownames(fixed_df) <- names(fixed_values)
    }
  }

  fixedValues <- numeric()


  this_call <- match.call()
  initControlList <- controlList

  aggr <- TRUE
  return_df <- (return_type == "df")
  return_first <- (return_type == "first")
  return_lik_df <- FALSE
  if(return_type == "predict") {
    aggr <- FALSE
    return_lik_df = TRUE
    return_type = "first_df"
    #print("GOTHERE")
    # print(fixed_values)
  }
  if(return_type == "first_df") {
    return_first <- TRUE
    return_first_df <- TRUE
  } else {
    return_first_df <- FALSE
  }
  if(return_type == "precalc_df")  {
    return_first <- TRUE
    return_lik_df <- FALSE
    return_first_df <- TRUE
  }

  #print(return_lik_df)

  start_vals = vector()
  prep_out <- list()

  newDataList <- list()
  all_defined_vars <- character()
  undefined_vars <- character()

  for(dataName in names(dataList)) {
    rel_col <- vector()
    control <- controlList[[dataName]]
    data_df <- dataList[[dataName]]
    orig_cols <- colnames(data_df)
    valueVar <- NA
    defined_vars <- colnames(data_df)
    for(formula in control$formulas) {
      targetName <- as.character(formula[[2]])
      rel_col <- c(rel_col, all.vars(formula)[!all.vars(formula) %in% rel_col])
      if(is.na(valueVar)) if(targetName %in% colnames(data_df) | (length(grep(paste0(targetName, "\\."), colnames(data_df), value = TRUE)) == 2)) valueVar <- targetName
      rhs_vars <- all.vars(formula[[3]])
      undefined_vars <- c(undefined_vars, rhs_vars[!rhs_vars %in% defined_vars])
      all_defined_vars <- c(all_defined_vars, undefined_vars[!undefined_vars %in% all_defined_vars])
      if(!targetName %in% c(colnames(data_df), valueVar)) defined_vars <- c(defined_vars, targetName)
    }
    censor_bounds <- control$censor_bounds
    #print(colnames(data_df))
    #print(valueVar)
    grep_vals <- grep(paste0(valueVar, "\\."), colnames(data_df), value = TRUE)
    if(length(grep_vals) > 2) stop("More than two columns matching target variable.")
    if(length(grep_vals ) == 2) {
      # Two values identified. Check if one is allways >= the other
      if(all((data_df[,grep_vals[1]]-data_df[,grep_vals[2]])>=0)) {
        upper_bound_var <- grep_vals[1]
        lower_bound_var <- grep_vals[2]
      } else if(all((data_df[,grep_vals[1]]-data_df[,grep_vals[2]])<=0)) {
        upper_bound_var <- grep_vals[2]
        lower_bound_var <- grep_vals[1]
      } else stop(paste0("Two columns matched the left-hand side: ", grep_vals[1], " and ", grep_vals[2], ". However, neither was allways >= the other, so they cannot be used as upper/lower bounds. Quitting."))

    } else {
      if(!valueVar %in% colnames(data_df)) {
        stop("No value columns specified. ")
      } else {
        upper_bound_var <- lower_bound_var <- valueVar
      }
    }
    controlList[[dataName]]$upper_bound_var <- upper_bound_var
    controlList[[dataName]]$lower_bound_var <- lower_bound_var
    controlList[[dataName]]$valueVar <- valueVar

    data_df[, c("internal_ub", "internal_lb")] <- data_df[, c(upper_bound_var, lower_bound_var)]


    formula[[2]] <- substitute(value)
    e_params <- list(...)

    if(length(latent_id_colname)) {
      if(latent_id_colname %in% colnames(data_df)) data_df$internal_classgroup <- data_df[, latent_id_colname]
    }
    else data_df$internal_classgroup <- paste0(dataName, "." , 1:NROW(data_df))

    data_df <- within(data_df, {
      both_na <- is.na(internal_ub) * is.na(internal_lb)

      internal_ub[internal_ub >= max(control$censor_bounds)] <- Inf
      internal_lb[internal_lb <= min(control$censor_bounds)] <- -Inf

      internal_ub[is.na(internal_ub)] <- Inf
      internal_lb[is.na(internal_lb)] <- -Inf

      internal_count <- 1

    })

    #print(data_df)

    if(sum(data_df$both_na)) {
      warning(paste0(sum(data_df$both_na," rows in which ", valueVar, " was NA were removed.")))
      data_df <- data_df[!data_df$both_na, -(NCOL(data_df))]
    }
    if(!is.na((weights_var <- control$weights_var)) & weights_var %in% colnames(data_df)) data_df$internal_count <- data_df[, weights_var]

    relevant_columns <- c(rel_col, colnames(data_df)[grep("internal_", colnames(data_df))])

    if(aggr) {
      new_df <- data_df[, colnames(data_df) %in% relevant_columns]
      new_df$internal_unique <- as.numeric(as.factor(apply(new_df, MARGIN = 1, FUN = "paste", collapse = "")))
      unique_df <- unique(new_df)
      unique_df$internal_count <- aggregate(new_df$internal_count, by = list(new_df$internal_unique), FUN = "sum")[unique_df$internal_unique, 2]
      newDataList[[dataName]] <- unique_df
    } else {
      newDataList[[dataName]] <- data_df

    }
    controlList[[dataName]]$orig_cols <- orig_cols

    startValues <- c(startValues, control$start_values[!names(control$start_values) %in% names(startValues)])
    fixedValues <- c(fixedValues, control$fixed_values[!names(control$fixed_values) %in% names(fixedValues)])

    all_defined_vars <- c(all_defined_vars, defined_vars[!defined_vars %in% all_defined_vars])
  }


  if(length(fixedValues)) {

    fixed_values[names(fixedValues)] <- fixedValues
    tmp <- rep(NA, NROW(fixed_values))
    names(tmp) <- names(fixed_values)
    tmp[rownames(fixed_df)] <- fixed_df[,2]
    fixed_df <- data.frame(Estimate = fixed_values, `Std. Error` = tmp)

  }
  startValues[names(fixed_values)[names(fixed_values) %in% names(startValues)]] <- fixed_values[names(fixed_values) %in% names(startValues)]


  # if(length(start_suggestions)) {
  #   if("xreg" %in% class(start_suggestions)) start_suggestions <- start_suggestions$coef
  #   if(is.numeric(start_suggestions)) {
  #     overlap <- names(startValues)[names(startValues) %in% names(start_suggestions)]
  #     startValues[overlap] <- start_suggestions[overlap]
  #   }
  # }

  undefined_vars <- unique(undefined_vars[!undefined_vars %in% c(names(startValues), latent_class_parameters)])
  if(length(undefined_vars)) {
    tmp <- paste(undefined_vars, collapse = ", ")
    used_from <- "xreg"
    if("run_from" %in% names(dot_args)) used_from <- dot_args[["run_from"]]
    if(!used_from == "hyreg") warning(paste0("Used formulas contained variables not found in start_values or column names (", tmp,"). These will be initiated with a starting value of 0.1. Defining start values either in xregControl or in the call to xreg is highly recommended."))
    tmp <- rep(0.1, length(undefined_vars))
    names(tmp) <- undefined_vars
    startValues <- c(startValues, tmp)
  }


  if(latent_classes > 1) {
    for(lc in latent_class_parameters) {
      tmp <- 0.1
      names(tmp) <- lc
      if(lc %in% names(startValues)) {
        tmp <- startValues[lc]
        startValues <- startValues[-which(names(startValues) == lc)]
        all_defined_vars <- all_defined_vars[-which(all_defined_vars== lc)]
      }
      for(i in 1:latent_classes) {
        lci <- paste(lc, i, sep = ".")
        if(!lci %in% names(startValues)) {
          startValues <- c(startValues, tmp)
          names(startValues)[length(startValues)] <- lci
        }
        if(!lci %in% all_defined_vars) {
          all_defined_vars <- c(all_defined_vars, lci)
        }
      }
    }
  }

  if(length((tmp <- names(startValues[!names(startValues) %in% all_defined_vars])))) {
    startValues <- startValues[!names(startValues) %in% tmp]
    tmp <- paste(tmp, collapse = ", ")
    warning(paste0("Start values defined for variables (", tmp, ") not found in column names or formulas. These have been removed."))
  }




  if("control" %in% names(e_params)){
    control <- e_params[["control"]]
  } else {
    boxconstr <- F
    if("method" %in% names(list(...))) {
      meth <- list(...)[["method"]]
      if(meth == "L-BFGS-B") boxconstr <- T
    }
    if(boxconstr) control <- list(maxit = 10000, factr = .Machine$double.eps)
    else control <- list(maxit = 10000, abstol = .Machine$double.eps, reltol = .Machine$double.eps)
    #control <- list(maxit = 5, abstol = .Machine$double.eps, reltol = .Machine$double.eps)
  }

  startValues[names(fixed_values)] <- fixed_values

  if(return_first) {
    ret <- do.call(lik_fun, c(startValues,
                              list(
                                controlList = controlList,
                                dataList = newDataList,
                                print_sum = F,
                                return_lik_df = return_lik_df,
                                return_first_df = return_first_df,
                                latent_classes = latent_classes,
                                latent_class_parameters = latent_class_parameters,
                                dot_args)))
    for(nm in names(ret)) {
      #ret[[nm]] <- ret[[nm]][, c(orig_cols, "Xb", "p")]
    }

    return(ret)
  }

  testmle <- do.call(x_mle, c(list(minuslogl = lik_fun,
                                   start = as.list(startValues),
                                   controlList = controlList,
                                   dataList = newDataList,
                                   control = control,
                                   print_sum = print_sum,
                                   fixed = as.list(fixed_values),
                                   latent_classes = latent_classes,
                                   latent_class_parameters = latent_class_parameters),
                              dot_args))

  res_types <- numeric()

  dfs <- do.call(lik_fun, c(c(testmle$coef,fixed_values),
                            list(
                              controlList = controlList,
                              dataList = newDataList,
                              print_sum = F,
                              return_lik_df = T,
                              return_first_df = T,
                              latent_classes = latent_classes,
                              latent_class_parameters = latent_class_parameters,
                              dot_args)))


  for(nm in names(dfs)) {

    res_types[nm] <- sum(controlList[[dataName]]$p_aggregation_fun(dfs[[nm]]))
  }
  res_types["sum"] <- sum(res_types)

  testmle$fixed_values <- fixed_df
  testmle$minima <- res_types
  pars <- rbind(cbind(testmle$fullcoef, type = rep("Fitted", NROW(testmle$fullcoef))),cbind(fixed_df, type = rep("Fixed", NROW(fixed_df))))
  testmle$pars <- pars
  testmle$pars_v <- pars[,1]
  names(testmle$pars_v) <- rownames(pars)
  res <- list(controlList = controlList,
              start_values = start_values,
              mle_obj = testmle,
              coef = testmle$coef,
              coefficients = testmle$coef,
              full_coef = testmle$fullcoef,
              valueVar = valueVar,
              fixed_values = fixed_df,
              minima = res_types,
              pars = pars,
              pars_v = testmle$pars_v)

  class(res) <- c("xreg", "list")


  return(res)

}