###################################################################
# TMLE targeting step:                                            #
# estimate each treatment-specific mean across J treatment levels #
# https://github.com/cran/tmle/blob/master/R/tmle.R               #
###################################################################

## calculate ATTs and CIs for pairwise comparisons 
getTMLE <- function(QAW, g_preds_bounded, obs.treatment, Y, outcome.type=c("continuous","binomial"), target.gwt=TRUE){
  
  # create lists/vectors to store results
  weights <- vector(mode = "list", length = 1)
  
  # compute clever covariate
  
  if(target.gwt){ # propensity weights moved from denominator of clever covariate to regression weight
    # for ATE (target pop. is all units)
    weights[[1]] <-  rowSums(obs.treatment/g_preds_bounded)
    
    clever_covariates <- obs.treatment
  }else{ # propensity weights in denominator of clever covariate
    # for ATE
    weights[[1]]  <- rep(1, nrow(obs.treatment))
    
    clever_covariates <- obs.treatment/g_preds_bounded
  }
  
  # targeting step - refit outcome model using clever covariates
  if(outcome.type=="binomial"){
    updated_model_for_Y <- glm(Y~-1 + clever_covariates + offset(qlogis(QAW$QA)), weights=if(target.gwt){as.vector(scale(weights[[1]], center=FALSE))} else{weights[[1]]}, family=quasibinomial()) #  scale weights to avoid convergence problems 
  }else{
    updated_model_for_Y <- glm(Y~-1 + clever_covariates + offset(QAW$QA), weights=if(target.gwt){as.vector(scale(weights[[1]], center=FALSE))} else{weights[[1]]}, family=gaussian())  #  scale weights to avoid convergence problems
  }
  
  # for ATE
  if(target.gwt){
    if(outcome.type=="binomial"){
      Qstar <- data.frame(mapply(qlogis, QAW) + c(rowSums(coef(updated_model_for_Y)*clever_covariates), sapply(1:J, function(j){
        rep(coef(updated_model_for_Y)[j], length(Y))})))
    }else{
      Qstar <- data.frame(QAW + c(rowSums(coef(updated_model_for_Y)*clever_covariates), sapply(1:J, function(j){
        rep(coef(updated_model_for_Y)[j], length(Y))})))
    }
  }else{
    if(outcome.type=="binomial"){
      Qstar <- data.frame(mapply(qlogis, QAW) + c(rowSums(coef(updated_model_for_Y)*clever_covariates), coef(updated_model_for_Y)/g_preds_bounded))
    }else{
      Qstar <- data.frame(QAW + c(rowSums(coef(updated_model_for_Y)*clever_covariates), coef(updated_model_for_Y)/g_preds_bounded))
    }
  }
 
  return(Qstar)
}