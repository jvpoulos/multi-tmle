############################################################################################
# Static setting (T=1) simulations: Compare multinomial TMLE with binary TMLE             #
############################################################################################

# Setup parallel processing
doMPI <- TRUE
if(doMPI){
  library(doMPI)
  
  # Start cluster
  cl <- startMPIcluster(verbose=TRUE)
  
  # Register cluster
  registerDoMPI(cl)
  
  # Check cluster size
  print(paste0("cluster size: ", clusterSize(cl)))
  
} else{
  library(parallel)
  library(doParallel)
  library(foreach)
  
  cores <- parallel::detectCores()
  print(paste0("number of cores used: ", cores))
  
  cl <- parallel::makeCluster(cores, outfile="")
  
  doParallel::registerDoParallel(cl) # register cluster
}

# check directory
print(getwd())

# load utils
source('./src/data_generation.R')
source('./src/tmle_calculation.R')
source('./src/tmleContrast.R')
source('./src/misc_fns.R')

#######################
# Simulation function #
######################

staticSim <- function(r, J, n, gbound, ybound, n.folds, estimator, overlap.setting, gamma.setting, outcome.type, target.gwt, use.SL, scale.continuous){
  
  library(purrr)
  library(origami)
  library(sl3)
  library(nnet)
  library(ranger)
  library(xgboost)
  library(glmnet)
  library(MASS)
  library(lmtp)
  
  print(paste("This is simulation run number",r, "\n"))
  
  if(!J%in%c(3,6)){
    stop("J must be 3 or 6")
  }
  
  ## Generate data
  generated.data <- generateData(r, J, n, overlap.setting, gamma.setting, outcome.type, scale.continuous)
  true.ates <- generated.data$trueATE
  obs.treatment <- generated.data$observed.treatment
  
  C <- generated.data$C # contrast matrix
  
  Y <- generated.data$data$Y # extract outcome
  
  A <- as.factor(generated.data$data$Z) # extract treatment
  
  # define a dataframe of covariates
  L <- generated.data$data[!colnames(generated.data$data)%in%c("Y","Z")] 
  rm(generated.data)
  
  # store mean observed outcome under treatment, mean treatment, and true ATEs
  obs_outcome <- colSums(Y*obs.treatment)/colSums(obs.treatment)
  obs_treatment <- colMeans(obs.treatment)
  obs_covariates <- colMeans(L)
  
  ## Manual TMLE (ours)
  
  yinitial_tmle <-NULL # returns NULL if lmtp is estimated
  Ahat_tmle <- NULL
  yhat_tmle <- NULL
  ess_ate_tmle <- NULL
  Ahat_tmle_bin <- NULL
  yhat_tmle_bin <- NULL
  ess_ate_tmle_bin <- NULL
  
  bias_ate_tmle <- list()
  var_ate_tmle <- list()
  CP_ate_tmle <- list()
  CIW_ate_tmle <- list()
  
  bias_ate_tmle_bin <- list()
  var_ate_tmle_bin <- list()
  CP_ate_tmle_bin <- list()
  CIW_ate_tmle_bin <- list()
  
  # stack learners into a model
  learner_stack_A <- make_learner_stack(list("Lrnr_xgboost",nrounds=20, objective="multi:softprob", eval_metric="mlogloss",num_class=J), list("Lrnr_ranger",num.trees=100),list("Lrnr_ranger",num.trees=500), list("Lrnr_glmnet",nfolds = n.folds,alpha = 1, family = "multinomial"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.25, family = "multinomial"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.5, family = "multinomial"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.75, family = "multinomial"))  
  learner_stack_A_bin <- make_learner_stack(list("Lrnr_xgboost",nrounds=20, objective = "reg:logistic"), list("Lrnr_ranger",num.trees=100),list("Lrnr_ranger",num.trees=500), list("Lrnr_glmnet",nfolds = n.folds,alpha = 1, family = "binomial"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.25, family = "binomial"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.5, family = "binomial"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.75, family = "binomial"))  
  
  # metalearner defaults 
  if(outcome.type=="continuous"){
    metalearner_Y <- make_learner(Lrnr_solnp,learner_function=metalearner_linear,eval_function=loss_squared_error) # nonlinear optimization via augmented Lagrange
    learner_stack_Y <- make_learner_stack(list("Lrnr_xgboost",nrounds=20, objective = "reg:squarederror"), list("Lrnr_ranger",num.trees=100), list("Lrnr_ranger",num.trees=500), list("Lrnr_glmnet",nfolds = n.folds,alpha = 1, family = "gaussian"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.25, family = "gaussian"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.5, family = "gaussian"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.75, family = "gaussian")) 
  }else if (outcome.type =="binomial"){
    metalearner_Y <- make_learner(Lrnr_solnp,learner_function=metalearner_logistic_binomial,eval_function=loss_loglik_binomial)
    learner_stack_Y <- make_learner_stack(list("Lrnr_xgboost",nrounds=20, objective = "reg:logistic"), list("Lrnr_ranger",num.trees=100), list("Lrnr_ranger",num.trees=500), list("Lrnr_glmnet",nfolds = n.folds,alpha = 1, family = "binomial"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.25, family = "binomial"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.5, family = "binomial"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.75, family = "binomial")) 
  }
  metalearner_A <- make_learner(Lrnr_solnp,learner_function=metalearner_linear_multinomial, eval_function=loss_loglik_multinomial)
  metalearner_A_bin <- make_learner(Lrnr_solnp,learner_function=metalearner_logistic_binomial,eval_function=loss_loglik_binomial)
  
  if(estimator=="tmle"){
    
    ## fit initial outcome model
    if(use.SL){
      
      # define task and candidate learners
      initial_model_for_Y_task <- make_sl3_Task(cbind(Y,L,obs.treatment[,-1]), covariates = c(colnames(L),colnames(obs.treatment[,-1])), outcome = "Y", outcome_type=outcome.type, 
                                                folds = origami::make_folds(cbind(Y,L,obs.treatment[,-1]), fold_fun = folds_vfold, V = n.folds)) # A1 is reference
      
      # Super Learner algorithm
      
      initial_model_for_Y_sl <- make_learner(Lrnr_sl, # cross-validates base models
                                             learners = learner_stack_Y,
                                             metalearner = metalearner_Y,
                                             keep_extra=FALSE)
      initial_model_for_Y_sl_fit <- initial_model_for_Y_sl$train(initial_model_for_Y_task)
      initial_model_for_Y_preds <- initial_model_for_Y_sl_fit$predict(initial_model_for_Y_task) # predicted probs.
      
      for(j in 1:J){
        newdata <- cbind(Y,L,matrix(0, nrow=nrow(obs.treatment), ncol=ncol(obs.treatment), dimnames = dimnames(obs.treatment)))
        newdata[paste0("D",j)] <- 1
        assign(paste0("Q",j), initial_model_for_Y_sl_fit$predict(sl3_Task$new(newdata, covariates = c(colnames(L),colnames(obs.treatment)), outcome = "Y", outcome_type="binomial")))
      }
    } else{
      if(outcome.type=="continuous"){
        initial_model_for_Y <- glm(Y~.,data=cbind(Y,obs.treatment[,-1],L), family=gaussian(), control = glm.control(maxit = 100)) # A1 is reference
        initial_model_for_Y_preds <- predict(initial_model_for_Y)
      }else{
        initial_model_for_Y <- glm(Y~.,data=cbind(Y,obs.treatment[,-1],L), family=binomial(), control = glm.control(maxit = 100)) # A1 is reference
        initial_model_for_Y_preds <- predict(initial_model_for_Y,type="response") # predicted probs.
      }
      
      for(j in 1:J){
        newdata <- cbind(Y,L,matrix(0, nrow=nrow(obs.treatment), ncol=ncol(obs.treatment), dimnames = dimnames(obs.treatment)))
        newdata[paste0("D",j)] <- 1
        assign(paste0("Q",j), predict(initial_model_for_Y,newdata=newdata, type="response"))
      }
    }
    rm(newdata)
    
    Qs <- t(do.call(mapply, c(FUN = cbind, mget(paste0("Q", 1:J)))))
    colnames(Qs) <- paste0("Q", 1:J)
    
    QAW <- data.frame(apply(cbind(QA=initial_model_for_Y_preds,Qs), 2, boundProbs, bounds=ybound)) # bound predictions
    
    ##  fit initial treatment model
    
    if(use.SL){
      
      # multinomial
      
      initial_model_for_A_task <- make_sl3_Task(cbind(A,L), covariates = c(names(L)), outcome = "A", outcome_type="categorical", 
                                                folds = origami::make_folds(cbind(A,L), fold_fun = folds_vfold, V = n.folds))
      
      initial_model_for_A_sl <- make_learner(Lrnr_sl, # cross-validates base models
                                             learners = learner_stack_A,
                                             metalearner = metalearner_A,
                                             keep_extra=FALSE)
      initial_model_for_A_sl_fit <- initial_model_for_A_sl$train(initial_model_for_A_task)
      
      g_preds  <- initial_model_for_A_sl_fit$predict(initial_model_for_A_task)  # estimated propensity scores 
      
      g_preds  <-  data.frame(matrix(unlist(lapply(g_preds, unlist)), nrow=length(lapply(g_preds, unlist)), byrow=TRUE)) 
      
      colnames(g_preds) <- colnames(obs.treatment)
      
      # binomial
      
      initial_model_for_A_task_bin <- lapply(1:J, function(j) make_sl3_Task(cbind("A"=obs.treatment[,j],L), covariates = c(names(L)), outcome = "A", outcome_type="binomial",
                                                                            folds = origami::make_folds(cbind(A,L), fold_fun = folds_vfold, V = n.folds)))
      
      initial_model_for_A_sl_bin <- make_learner(Lrnr_sl, # cross-validates base models
                                                 learners = learner_stack_A_bin,
                                                 metalearner = metalearner_A_bin,
                                                 keep_extra=FALSE)
      initial_model_for_A_sl_fit_bin <- lapply(1:J, function(j) initial_model_for_A_sl_bin$train(initial_model_for_A_task_bin[[j]]))
      
      g_preds_bin  <- sapply(1:J, function(j) initial_model_for_A_sl_fit_bin[[j]]$predict(initial_model_for_A_task_bin[[j]]))  # estimated propensity scores 
      
      colnames(g_preds_bin) <- colnames(obs.treatment)
      
    } else{
      
      # multinomial logistic regression
      initial_model_for_A <- nnet::multinom(formula=A ~ ., data=cbind(A, L), maxit = 500, trace = FALSE) 
      g_preds <- fitted(initial_model_for_A)
      colnames(g_preds) <- colnames(obs.treatment)
      
      # binomial logistic regression
      initial_model_for_A_bin <- lapply(1:J, function(j) glm(formula=treatment ~ ., data=cbind("treatment"=obs.treatment[,j], L), family="binomial"))
      g_preds_bin <- sapply(1:J, function(j) predict(initial_model_for_A_bin[[j]], type="response"))
      colnames(g_preds_bin) <- colnames(obs.treatment)
    }
    
    g_preds_bounded <- boundProbs(g_preds, gbound) # (trimmed) estimated propensity scores
    g_preds_bin_bounded <- boundProbs(g_preds_bin, gbound)
    
    # Update equations, calculate point estimate and variance
    
    TMLE <- getTMLE(QAW, g_preds_bounded, obs.treatment, Y, outcome.type, target.gwt)
    TMLE_bin <- getTMLE(QAW, g_preds_bin_bounded, obs.treatment, Y, outcome.type, target.gwt)
    
    # calculate contrasts
    tmle_contrasts <- tmle_contrast(Qstar=TMLE, C, g_preds_bounded, obs.treatment, Y, outcome.type, alpha=0.05, multi.adjust=FALSE)
    tmle_contrasts_bin <- tmle_contrast(Qstar=TMLE_bin, C, g_preds_bin_bounded, obs.treatment, Y, outcome.type, alpha=0.05, multi.adjust=FALSE)
    
    # store results
    
    yinitial_tmle  <- QAW[-1] # estimated initial Y, treatment probs. 
    Ahat_tmle  <- g_preds
    Ahat_tmle_bin  <- g_preds_bin
    
    yhat_tmle <- colMeans(tmle_contrasts$Qstar)
    yhat_tmle_bin <- colMeans(tmle_contrasts_bin$Qstar)
    
    ess_ate_tmle <- (colSums(obs.treatment*(1/g_preds))**2)/(colSums(obs.treatment*(1/g_preds)**2)) 
    ess_ate_tmle_bin <- (colSums(obs.treatment*(1/g_preds_bin))**2)/(colSums(obs.treatment*(1/g_preds_bin)**2)) 
    
    bias_ate_tmle  <- tmle_contrasts$taus$ATE- true.ates # point and variance estimates for ATE (multinomial)
    var_ate_tmle <- sapply(tmle_contrasts$var,"[[",1)
    CP_ate_tmle <- as.numeric((sapply(tmle_contrasts$CI,"[[",1)[1,] < true.ates) & (sapply(tmle_contrasts$CI,"[[",1)[2,] > true.ates))
    CIW_ate_tmle  <- sapply(tmle_contrasts$CI,"[[",1)[2,] -sapply(tmle_contrasts$CI,"[[",1)[1,] 
    
    bias_ate_tmle_bin  <- tmle_contrasts_bin$taus$ATE- true.ates # point and variance estimates for ATE (multinomial)
    var_ate_tmle_bin <- sapply(tmle_contrasts_bin$var,"[[",1)
    CP_ate_tmle_bin <- as.numeric((sapply(tmle_contrasts_bin$CI,"[[",1)[1,] < true.ates) & (sapply(tmle_contrasts_bin$CI,"[[",1)[2,] > true.ates))
    CIW_ate_tmle_bin  <- sapply(tmle_contrasts_bin$CI,"[[",1)[2,] -sapply(tmle_contrasts_bin$CI,"[[",1)[1,] 
  }
  
  ## compare with lmtp
  
  yinitial_lmtp <- list()
  Ahat_lmtp <- list()
  ess_ate_lmtp <- list()
  
  results_ate_bias_lmtp <- list() # ATE results
  results_ate_var_lmtp <- list()
  results_ate_CP_lmtp <- list()
  results_ate_CIW_lmtp <- list()
  
  if(estimator=="lmtp"){
    
    lmtp_results <- list()
    
    # estimate outcomes under each J
    for(j in c(1:J)){      
      lmtp_results[[j]] <- lmtp_tmle(data = cbind(Y,A,L), 
                                     trt = "A",
                                     outcome = "Y", 
                                     baseline = colnames(cbind(Y,A,L)[grep("X",colnames(cbind(Y,A,L)))]), 
                                     shifted = data.frame(Y,"A"=factor(rep(j,times=length(A)), levels=levels(A)),L),
                                     intervention_type = "mtp",
                                     outcome_type = outcome.type, 
                                     folds = n.folds,
                                     .bound = if(outcome.type=="binomial") ybound[1] else 1e-5,
                                     .trim = gbound[2],
                                     .SL_folds = n.folds,
                                     learners_trt= if(use.SL) learner_stack_A_bin else sl3::make_learner(sl3::Lrnr_glm),
                                     learners_outcome= if(use.SL) learner_stack_Y else sl3::make_learner(sl3::Lrnr_glm))
    }
    
    # calculate contrasts
    lmtp_contrasts <- list()
    for(i in 1:nrow(C)){
      lmtp_contrasts[[i]] <- lmtp_contrast(lmtp_results[[as.numeric(colnames(C)[which(C[i,]==1)])]], ref=lmtp_results[[as.numeric(colnames(C)[which(C[i,]==-1)])]], type='additive') # treated - reference
    }
    
    # store results
    
    for(j in 1:J){
      yinitial_lmtp[[j]] <- lmtp_results[[j]]$outcome_reg
    }
    
    for(j in 1:J){
      Ahat_lmtp[[j]] <- (mean(obs.treatment[,j])*lmtp_results[[j]]$density_ratios)/(1+mean(obs.treatment[,j])*lmtp_results[[j]]$density_ratios) # convert density ratio to propensity score
    }
    
    for(j in 1:J){
      ess_ate_lmtp[[j]] <- (sum(obs.treatment[,j]*(1/Ahat_lmtp[[j]]))**2)/(sum(obs.treatment[,j]*(1/Ahat_lmtp[[j]])**2))
    }
    
    for(i in 1:nrow(C)){
      results_ate_bias_lmtp[[i]]  <- lmtp_contrasts[[i]]$vals$theta - true.ates[i] # LMTP ATE point and variance estimates
      results_ate_var_lmtp[[i]]  <- sqrt(lmtp_contrasts[[i]]$vals$std.error/n)
      results_ate_CP_lmtp[[i]]  <- as.numeric((lmtp_contrasts[[i]]$vals$conf.low < true.ates[i]) & (lmtp_contrasts[[i]]$vals$conf.high > true.ates[i]))
      results_ate_CIW_lmtp[[i]] <- lmtp_contrasts[[i]]$vals$conf.high-lmtp_contrasts[[i]]$vals$conf.low
    }
  }
  
  # cleanup
  gc()
  print(paste("Done with simulation run number",r, "\n"))
  
  # return results
  return(list("obs_outcome"=obs_outcome,"obs_treatment"=obs_treatment,"obs_covariates"=obs_covariates,"trueATE"=true.ates,
              "yinitial_tmle"=yinitial_tmle,"Ahat_tmle"=Ahat_tmle,"yhat_tmle"= yhat_tmle, "ess_ate_tmle"=ess_ate_tmle,
              "Ahat_tmle_bin"=Ahat_tmle_bin,"yhat_tmle_bin"= yhat_tmle_bin, "ess_ate_tmle_bin"=ess_ate_tmle_bin,
              "bias_ate_tmle"= unlist(bias_ate_tmle),"var_ate_tmle"= unlist(var_ate_tmle),"CP_ate_tmle"=unlist(CP_ate_tmle),"CIW_ate_tmle"=unlist(CIW_ate_tmle),
              "bias_ate_tmle_bin"= unlist(bias_ate_tmle_bin),"var_ate_tmle_bin"= unlist(var_ate_tmle_bin),"CP_ate_tmle_bin"=unlist(CP_ate_tmle_bin),"CIW_ate_tmle_bin"=unlist(CIW_ate_tmle_bin),
              "yinitial_lmtp"=unlist(yinitial_lmtp),"Ahat_lmtp"=unlist(Ahat_lmtp), "ess_ate_lmtp"=ess_ate_lmtp,
              "bias_ate_lmtp"=unlist(results_ate_bias_lmtp),"var_ate_lmtp"=unlist(results_ate_var_lmtp),"CP_ate_lmtp"=unlist(results_ate_CP_lmtp),"CIW_ate_lmtp"=unlist(results_ate_CIW_lmtp)))
}

#####################
# Set parameters    #
#####################

# define settings for simulation
settings <- expand.grid("J"=c(3,6),
                        "n"=c(6000),
                        "overlap.setting"=c("adequate","inadequate","rct"),
                        "gamma.setting"=c("zero","yang","low"))
settings$n <- ifelse(settings$J==6, 6000, 1500)

vary.n <- FALSE
if(vary.n){
  settings <- expand.grid("J"=c(6),
                          "n"=c(600,10000),
                          "overlap.setting"=c("adequate","inadequate","rct"),
                          "gamma.setting"=c("zero","yang","low"))
}

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE) # command line arguments
estimator <- as.character(args[1])
thisrun <- settings[as.numeric(args[2]),]
outcome.type <-  as.character(args[3]) # "continuous" or "binomial" 
use.SL <- as.logical(args[4])  # When TRUE, use Super Learner for initial Y model and treatment model estimation; if FALSE, use GLM

# define parameters

J <- as.numeric(thisrun[1]) # number of treatments
n <- as.numeric(thisrun[2]) # total sample size
overlap.setting <- thisrun[[3]] # "adequate"= adequate overlap scenario of Li and Li (2019); "inadequate" overlap scenario of Yang et al. (2016); "rct": kappa values all zero
gamma.setting <- thisrun[[4]] # zero"= gamma values are all zero, so that there is no treatment effect; "yang"= use values from Yang et al. (2016); "li" use values from Li and Li (2019)

R <- 1000 # number of simulation runs

gbound <- c(0.001,0.999) # define bounds to be used for the propensity score

if(outcome.type=="continuous"){
  ybound <- c(0,1000)
}else{
  ybound <- c(0.0001, 0.9999) # define bounds to be used for the Y predictions
}

n.folds <- 5

target.gwt <- TRUE # When TRUE, moves propensity weights from denominator of clever covariate to regression weight when fitting updated model for Y

scale.continuous <- TRUE # standardize continuous covariates

output_dir <- './outputs/'
simulation_version <- paste0(format(Sys.time(), "%Y%m%d"),"/")
if(!dir.exists(output_dir)){
  print(paste0('create folder for outputs at: ', output_dir))
  dir.create(output_dir)
}
output_dir = paste0(output_dir, simulation_version)
if(!dir.exists(output_dir)){
  print(paste0('create folder for outputs at: ', output_dir))
  dir.create(output_dir)
}

filename <- paste0(output_dir, 
                   "static_simulation_results_",
                   "estimator_",estimator,
                   "_R_", R,
                   "_n_", n,
                   "_J_", J,
                   "_n_folds_",n.folds,
                   "_overlap_setting_", overlap.setting,
                   "_gamma_setting_", gamma.setting,
                   "_outcome_type_", outcome.type,
                   "_use_SL_", use.SL,
                   "_scale_continuous_", scale.continuous,
                   "_target_gwt_", target.gwt,".rds")

#####################
# Run simulation #
#####################

print(paste0('simulation setting: ', "estimator = ", estimator, " R = ", R, ", n = ", n,", J = ",J, ", overlap.setting = ",overlap.setting, ", gamma.setting = ", gamma.setting, ", outcome.type = ", outcome.type, ", use.SL = ",use.SL, ", scale.continuous = ", scale.continuous))

sim.results <- foreach(r = 1:R, .combine='cbind', .verbose = TRUE) %dopar% {
  staticSim(r=r, J, n, gbound, ybound, n.folds, estimator, overlap.setting, gamma.setting, outcome.type, target.gwt, use.SL, scale.continuous)
}
sim.results

saveRDS(sim.results, filename)

if(doMPI){
  closeCluster(cl) # close down MPIcluster
  mpi.finalize()
}else{
  stopCluster(cl)
}