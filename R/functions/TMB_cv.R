# cross-validation for bespoke TMB models. Adapted from the sdmTMB_cv() function from
# the sdmTMB package (Anderson et al. 2022)

# testing
model_data=csyif_model1_v2_Data
data = csyif
time = "year"
model_parameters=csyif_model1_v2_Params
model_DLL= c("omega_s")
model_DLL = "csyif_model1_v2"



TMB_cv <- function(data,
                   time,
                   model_parameters,
                   model_random_effects,
                   model_DLL,
                   k_folds = 8,
                   parallel = TRUE,
                   model_obj,
                   data_observation_names){
  
  # function(
  #   formula, data, mesh_args, mesh = NULL, time = NULL,
  #   k_folds = 8, fold_ids = NULL,
  #   lfo = FALSE,
  #   lfo_forecast = 1,
  #   lfo_validations = 5,
  #   parallel = TRUE,
  #   use_initial_fit = FALSE,
  #   future_globals = NULL,
  #   ...) {
  
  # split data according to folds
      dd <- lapply(split(data, data[[time]]), function(x) {
        x$cv_fold <- sample(rep(seq(1L, k_folds), nrow(x)), size = nrow(x))
        x
      })
      data <- do.call(rbind, dd)
    fold_ids <- "cv_fold"
  
  if (k_folds > 1) {
    # data in kth fold get weight of 0:
    weights <- ifelse(data$cv_fold == 1L, 0, 1)
  } else {
    weights <- rep(1, nrow(data))
  }
  
  fit_func <- function(k) {
      # data in kth fold get weight of 0:
      weights <- ifelse(data$cv_fold == k, 0, 1)
      
      # set the model weights according to the folds
      cv_train_model_data <- model_data
      cv_train_model_data$weights_i <- weights
      
      ### fit a new TMB model
      cv_model_Obj = MakeADFun( data=cv_train_model_data,
                                       parameters=model_parameters,
                                       random= model_random_effects,
                                       DLL = model_DLL)
      
      # Optimize
      cv_model_Opt = nlminb( start=cv_model_Obj$par, obj=cv_model_Obj$fn, grad=cv_model_Obj$gr )
      # add getJointPrecision for index
      cv_model_Opt$SD = sdreport( cv_model_Obj, bias.correct=TRUE, getJointPrecision = TRUE )
      cv_model_v2_report = cv_model_Obj$report()
      
      ### Now, select the withheld data
      cv_data <- data[data$cv_fold == k, , drop = FALSE]
    
    # Predict model for the withheld data
    # predict for withheld data:
    predicted <- predict(object, newdata = cv_data, type = "response",
                         offset = if (!is.null(.offset)) cv_data[[.offset]] else rep(0, nrow(cv_data)))
    
    cv_data$cv_predicted <- predicted$est
    response <- get_response(object$formula[[1]])
    withheld_y <- predicted[[response]]
    withheld_mu <- cv_data$cv_predicted
    
    # FIXME: get LFO working with the TMB report() option below!
    # calculate log likelihood for each withheld observation:
    # trickery to get the log likelihood of the withheld data directly
    # from the TMB report():
    if (!lfo) {
      tmb_data <- object$tmb_data
      tmb_data$weights_i <- ifelse(tmb_data$weights_i == 1, 0, 1) # reversed
      new_tmb_obj <- TMB::MakeADFun(
        data = tmb_data,
        parameters = get_pars(object),
        map = object$tmb_map,
        random = object$tmb_random,
        DLL = "sdmTMB",
        silent = TRUE
      )
      lp <- object$tmb_obj$env$last.par.best
      r <- new_tmb_obj$report(lp)
      cv_loglik <- -1 * r$jnll_obs
      cv_data$cv_loglik <- cv_loglik[tmb_data$weights_i == 1]
    } else { # old method; doesn't work with delta models!
      cv_data$cv_loglik <- ll_sdmTMB(object, withheld_y, withheld_mu)
    }
    
    ## test
    # x2 <- ll_sdmTMB(object, withheld_y, withheld_mu)
    # identical(round(cv_data$cv_loglik, 6), round(x2, 6))
    # cv_data$cv_loglik <- ll_sdmTMB(object, withheld_y, withheld_mu)
    
    list(
      data = cv_data,
      model = object,
      pdHess = object$sd_report$pdHess,
      max_gradient = max(abs(object$gradients)),
      bad_eig = object$bad_eig
    )
  }
  
  if (requireNamespace("future.apply", quietly = TRUE) && parallel) {
    message(
      "Running fits with `future.apply()`.\n",
      "Set a parallel `future::plan()` to use parallel processing."
    )
    if (!is.null(future_globals)) {
      fg <- structure(TRUE, add = future_globals)
    } else {
      fg <- TRUE
    }
    if (lfo) {
      out <- future.apply::future_lapply(seq_len(lfo_validations), fit_func, future.seed = TRUE, future.globals = fg)
    } else {
      out <- future.apply::future_lapply(seq_len(k_folds), fit_func, future.seed = TRUE, future.globals = fg)
    }
  } else {
    message(
      "Running fits sequentially.\n",
      "Install the future and future.apply packages,\n",
      "set a parallel `future::plan()`, and set `parallel = TRUE` to use parallel processing."
    )
    if (lfo) {
      out <- lapply(seq_len(lfo_validations), fit_func)
    } else {
      out <- lapply(seq_len(k_folds), fit_func)
    }
  }
  
  models <- lapply(out, `[[`, "model")
  data <- lapply(out, `[[`, "data")
  fold_cv_ll <- vapply(data, function(.x) sum(.x$cv_loglik), FUN.VALUE = numeric(1L))
  data <- do.call(rbind, data)
  data <- data[order(data[["_sdm_order_"]]), , drop = FALSE]
  data[["_sdm_order_"]] <- NULL
  data[["_sdmTMB_time"]] <- NULL
  row.names(data) <- NULL
  pdHess <- vapply(out, `[[`, "pdHess", FUN.VALUE = logical(1L))
  max_grad <- vapply(out, `[[`, "max_gradient", FUN.VALUE = numeric(1L))
  converged <- all(pdHess)
  out <- list(
    data = data,
    models = models,
    fold_loglik = fold_cv_ll,
    sum_loglik = sum(data$cv_loglik),
    converged = converged,
    pdHess = pdHess,
    max_gradients = max_grad
  )
  `class<-`(out, "sdmTMB_cv")
}

log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

#' @export
#' @import methods
print.sdmTMB_cv <- function(x, ...) {
  nmods <- length(x$models)
  nconverged <- sum(x$pdHess)
  cat(paste0("Cross validation of sdmTMB models with ", nmods, " folds.\n"))
  cat("\n")
  cat("Summary of the first fold model fit:\n")
  cat("\n")
  print(x$models[[1]])
  cat("\n")
  cat("Access the rest of the models in a list element named `models`.\n")
  cat("E.g. `object$models[[2]]` for the 2nd fold model fit.\n")
  cat("\n")
  cat(paste0(nconverged, " out of ", nmods, " models are consistent with convergence.\n"))
  cat("Figure out which folds these are in the `converged` list element.\n")
  cat("\n")
  cat(paste0("Out-of-sample log likelihood for each fold: ", paste(round(x$fold_loglik, 2), collapse = ", "), ".\n"))
  cat("Access these values in the `fold_loglik` list element.\n")
  cat("\n")
  cat("Sum of out-of-sample log likelihoods:", round(x$sum_loglik, 2), "\n")
  cat("More positive values imply better out-of-sample prediction.\n")
  cat("Access this value in the `sum_loglik` list element.\n")
}