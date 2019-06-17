library(leaps)
library(nlme)

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

extended_gls <- function(form, data, weights) {
  model <- gls(form, data = data, weights = weights, method="ML")
  model$call$model <- form
  return(model)
}

extended_rand <- function(form, data, weights, rand_form) {
  model <- lme(form, data = data, weights = weights, random = rand_form, method="ML")
  model$call$fixed <- form
  return(model)
}

fit_all_models <- function(responses, predictors, gamma, data, rand_str, nbest=5, nvmax=5, eta=c(0, 0.5, 1)) {
  X <- data[, predictors]
  
  # Define model containers
  all_mods <- list()
  
  # Convert random string to a formula
  rand_form <- formula(rand_str)

  for (res in responses) {
    # Construct the leaps object
    y <- data[, res]
    lps <- leaps(x=X, y=y, nbest=nbest, method="adjr2")

    # Get a list of model formulas for all models with <= nvmax predictors
    lps_sub <- lps$which[which(as.numeric(rownames(lps$which)) <= nvmax),]
    forms_list <- apply(lps_sub, 1, make_model_forms, res, predictors)

    # Write the gls plots to file
    for (i in 1:length(eta)) {
      weight <- as.character(eta[[i]])
      # TODO generalize the form
      gls_list <- lapply(forms_list, extended_gls, data=data, weights = varPower(form=~p_95_mean, fixed=eta[[i]]))
      random_list <- lapply(forms_list, extended_rand, data=data, weights = varPower(form=~p_95_mean, fixed=eta[[i]]), rand_form=rand_form)
      
      all_mods[[res]][[weight]][['fixed']] <- gls_list
      all_mods[[res]][[weight]][['random_ind']] <- random_list
    }
  }
  return(all_mods)
}

write_fixed_plots <- function(fixed_list, out_dir) {
  for (res in names(fixed_list)) {
    res_dir <- paste(out_dir, res, sep="")
    dir.create(res_dir)
    
    for (hskdcty in names(fixed_list[[res]])) {
      models <- fixed_list[[res]][[hskdcty]]
      model_plots <- list()
      
      i <- 1
      for (model in models) {
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

get_model_obj <- function(list_mod, data, response) {
  # For some reason, I cannot call anova on the list elements, so this reconstructs the object
  predictors <- names(list_mod$coefficients[2:length(names(list_mod$coefficients))])
  model_frame <-  data[,c(predictors, response)]
  
  pred_form <- paste(predictors, collapse = ' + ')
  form <- paste(response, "~", pred_form)
}