# Function to calculate conditional AIC. Modified from the sdmTMB function

cAIC <- function(model_data,
                 model_obj,  
                 # model_opt,
                 model_parameters,
                 random,
                 DLL,
                 data_observation_names) {
  
  
  tmb_data <- model_data
  
  ## Ensure profile = NULL
  # if (is.null(object$control$profile)) {
  #   obj <- object$tmb_obj
  # } else {
  #   obj <- TMB::MakeADFun(
  #     data = tmb_data,
  #     parameters = object$parlist,
  #     map = object$tmb_map,
  #     random = object$tmb_random,
  #     DLL = "sdmTMB",
  #     profile = NULL #<
  #   )
  # }
  
  
  
  
  ## Make obj_new
  # weights = 0 is equivalent to data = NA
  # essentially what you are doing with obj_new is downweighting the data
  # to basically turn it off entirely
  # tmb_data$weights_i[] <- 0
  # the weights argument in sdmTMB is used to multiply the individual observations of species density
  # to weight those points: https://github.com/pbs-assess/sdmTMB/blob/60061922abffb9b6cd190ba51d83bc9d72f30ceb/src/sdmTMB.cpp
  # by making them NA here we are basically turning them off, which achieves the same thing
  # Using our argument data_observation_names, loop through all observations of 
  # species density and make them NA
  # for (i in 1:length(data_observation_names)){
  #   tmb_data[[data_observation_names[i]]] <- NA
  # }
  tmb_data$weights_i[] <- 0
  
  # get params
  params_vector <- lapply(split(as.list(model_obj$env$last.par.best), names(as.list(model_obj$env$last.par.best))), unlist)
  params_vector <- params_vector[names(model_parameters)]
  # for any parameters that should be a matrix, reformat them as such
  parameter_type_vector <- rep(NA, length(params_vector))
  params_list <- list()
  for (i in 1:length(params_vector)){
    parameter_type_vector[i] <- class(model_parameters[[i]])[1]
    if (parameter_type_vector[i] == "matrix"){
      params_list[[i]] <- matrix(params_vector[[i]], 
                                    nrow = nrow(model_parameters[[i]]),
                                    ncol = ncol(model_parameters[[i]]))
    } else {
      params_list[[i]] <- params_vector[[i]]
    }
  }
  names(params_list) <- names(params_vector)
  
  
  
  obj_new <- TMB::MakeADFun(
    data = tmb_data,
    # parameters = model_obj$env$parList(par = model_obj$env$last.par.best),
    # parameters = model_parameters,
    # parameters = model_obj$env$last.par.best,
    parameters = params_list,
    # parameters = lapply(split(as.list(model_opt[["par"]]), names(as.list(model_opt[["par"]]))), unlist),
    # map = list(),
    random = random,
    DLL = DLL,
    profile = NULL
  )
  
  par <- model_obj$env$parList()
  parDataMode <- model_obj$env$last.par
  indx <- model_obj$env$lrandom()
  q <- sum(indx)
  p <- length(model_obj$par)
  
  ## use '-' for Hess because model returns negative loglikelihood
  if (is.null(model_obj$env$random)) {
    print(c("This model has no random effects.", "cAIC and EDF only apply to models with random effects."))
    return(invisible(NULL))
  }
  Hess_new <- -Matrix::Matrix(obj_new$env$f(parDataMode, order = 1, type = "ADGrad"), sparse = TRUE)
  Hess_new <- Hess_new[indx, indx] ## marginal precision matrix of REs
  
  ## Joint hessian etc
  Hess <- -Matrix::Matrix(model_obj$env$f(parDataMode, order = 1, type = "ADGrad"), sparse = TRUE)
  Hess <- Hess[indx, indx]
  negEDF <- Matrix::diag(Matrix::solve(Hess, Hess_new, sparse = FALSE))
  
  jnll <- model_obj$env$f(parDataMode)
  cnll <- jnll - obj_new$env$f(parDataMode)
  cAIC_out <- 2 * cnll + 2 * (p + q) - 2 * sum(negEDF)
  return(cAIC_out)
}
