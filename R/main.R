library(leaps)
library(nlme)
library(stringr)
library(gridExtra)

write_leaps_plot <- function(leaps_obj, name, out_dir) {
  out_path <- paste(out_dir, name, ".pdf", sep="")
  pdf(out_path)
  plot(leaps_obj$size, leaps_obj$adjr2, ylim=c(0,1))
  dev.off()
}

make_model_forms <- function(row_bool, res, predictors) {
  pred_list <- predictors[row_bool]
  pred_form <- paste(pred_list, collapse = '+')
  pred_str <- paste(res, '~', pred_form, sep='')[[1]]

  return(formula(pred_str))
}

get_mcp_form <- function(mod_form) {
  mcp_str <- as.character(mod_form)[[3]] # The third element should always be the mcp (first is tilde, second is the response)
  mcp_str <- str_extract(mcp_str, "[^+]+") # Extract the response variable
  mcp_form <- paste("~", mcp_str)
  return(formula(mcp_form))
}

extended_gls <- function(form, data, eta) {
  weights <- varPower(fixed = eta, form = get_mcp_form(form))
  print(get_mcp_form(form))
  model <- gls(form, data = data, weights = weights, method="REML")
  model$call$model <- form
  return(model)
}

extended_rand <- function(form, data, rand_form, eta) {
  weights <- varPower(fixed = eta, form = get_mcp_form(form))
  model <- lme(form, data = data, weights = weights, random = rand_form, method="REML")
  model$call$fixed <- form
  return(model)
}

fit_all_models <- function(responses, predictors, data, rand_str, sp_str, nbest=5, nvmax=5, eta=c(0, 0.5, 1), verbose=FALSE,
                           lps_plot_out_dir='.') {
  X <- data[, predictors]
  X <- X[,  -which(names(X) %in% c(rand_str))]
  predictors <- predictors[-which(predictors == rand_str)]
  print(predictors)
  
  # Convert strigns to needed formulae
  rand_form <- formula(paste(" ~ 1|", rand_str, sep=""))
  sp_form <- formula(paste("~", sp_str, "|STAND", sep = ""))
  
  # Define model containers
  all_mods <- list()

  for (res in responses) {
    if (verbose) {cat(paste("Fitting models for:", res, "\n"))}
    
    # Construct the leaps object
    y <- data[, res]
    lps <- leaps(x=X, y=y, nbest=nbest, method="adjr2")
    
    # Write leaps plot
    write_leaps_plot(lps, res, lps_plot_out_dir)

    # Get a list of model formulas for all models with <= nvmax predictors
    lps_sub <- lps$which[which(as.numeric(rownames(lps$which)) <= nvmax),]
    print(lps_sub)
    forms_list <- apply(lps_sub, 1, make_model_forms, res, predictors)

    for (i in 1:length(eta)) {
      weight <- as.character(eta[[i]])
      
      if (verbose) {cat(paste("---eta=", eta[[i]], "\n", sep=""))}
      
      if (verbose) {cat(paste("------Fitting fixed models.\n"))}
      gls_list <- lapply(forms_list, extended_gls, data=data, eta=eta[[i]])
      
      if (verbose) {cat(paste("------Fitting ranef models.\n"))}
      random_ind_list <-  lapply(forms_list, extended_rand, data=data, eta=eta[[i]],
                                rand_form=rand_form)
      
      #if (verbose) {cat(paste("------Fitting ranef+spatial models.\n"))}
      #random_sp_list <-   lapply(forms_list, extended_rand_sp, data=data, weights = varPower(form=weights_form, fixed=eta[[i]]),
      #                          rand_form=rand_form, sp_form)
      
      all_mods[[res]][[weight]][['fixed']] <- gls_list
      all_mods[[res]][[weight]][['random_ind']] <- random_ind_list
      #all_mods[[res]][[weight]][['random_sp']] <- random_sp_list
  }
  return(all_mods)
}

write_fixed_plots <- function(fixed_list, out_dir) {
  for (res in names(fixed_list)) {
    res_dir <- paste(out_dir, res, sep="")
    dir.create(res_dir)
    
    for (hskdcty in names(fixed_list[[res]])) {
      models <- fixed_list[[res]][[hskdcty]]$fixed
      model_plots <- list()
      
      i <- 1
      for (model in models) {
        print(model)
        model_plots[[i]] <- plot.lme(model, abline=c(0,99999), main = as.character(ceiling(i / 2)))
        model_plots[[i+1]] <- qqnorm(model, abline=c(0,1))
        i <- i + 2
      }
      
        out_path <- paste(res_dir, '/', hskdcty, '.pdf', sep="")
        pdf(out_path, width=7, height=100)
        do.call("grid.arrange", c(model_plots, nrow=25, ncol=2))
        dev.off()
    }
  }
}

get_var_cov_bf <- function(obj, individuals, type) {
  if (attr(obj, "class") == "lme") {
    sigma <- obj$sigma
    D <- as.matrix(obj$modelStruct$reStruct[[1]]) * sigma^2
    result <- list()
    groups  <-  obj$groups[[1]]
    ugroups <- unique(groups)
    
    if (missing(individuals)) individuals  <-  as.matrix(ugroups)[1,]
    if (is.numeric(individuals))
      individuals  <-  ugroups[individuals]
    for (individ in individuals)
    {
      indx <- (1:length(ugroups))[individ==ugroups]
      if (!length(indx))
        stop(gettextf("individual %s was not used in the fit",
                      sQuote(individ)), domain = NA)
      if (is.na(indx))
        stop(gettextf("individual %s was not used in the fit",
                      sQuote(individ)), domain = NA)
      ind <- groups == individ
      if(!is.null(obj$modelStruct$corStruct)) {
        # FIXME generalize to different group identifiers
        n_i <- sum(obj$data$STAND == as.numeric(individuals))
        if (n_i > 1) {
          V <- corMatrix(obj$modelStruct$corStruct)[[as.character(individ)]]
        } else {
          V <- as.matrix(1)
        }
      }
      
      else V <- diag(sum(ind))
      if(!is.null(obj$modelStruct$varStruct))
        sds <- 1/varWeights(obj$modelStruct$varStruct)[ind]
      else
        sds <- rep(1, sum(ind))
      sds <- obj$sigma * sds
      cond.var <- t(V * sds) * sds
      dimnames(cond.var)  <-  list(1:nrow(cond.var),1:ncol(cond.var))
      if (type=="conditional")
        result[[as.character(individ)]] <- cond.var
      else
      {
        Z <- model.matrix(obj$modelStruct$reStruc,
                          getData(obj))[ind, , drop = FALSE]
        result[[as.character(individ)]] <-
          cond.var + Z %*% D %*% t(Z)
      }
    }
    class(result)  <-  c(type,"VarCov")
    attr(result,"group.levels")  <-  names(obj$groups)
    return(result[[individuals]])
  } else {
    stand_inds <- which(obj$data$STAND == as.numeric(individuals))
    sigma_sq <- obj$sigma^2
    n_i <- length(stand_inds)
    covariates <- attr(obj$modelStruct$varStruct, "weights")[stand_inds]
    eta <- attr(obj$modelStruct$varStruct, "fixed")
    cov_mat <- diag(sigma_sq * abs(covariates)^ (2 * eta), nrow=n_i, ncol=n_i)
    return(cov_mat)
  }
}
