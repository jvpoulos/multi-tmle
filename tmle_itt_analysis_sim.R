##################################
# ITT analysis with TMLE         #
# (not longitudinal)             #
##################################

library(purrr)
library(SuperLearner)
library(origami)
library(weights)
library(dplyr)
library(reporttools)
library(ggplot2)
library(reshape2)
library(rpart)
library(ranger)
library(glmnet)
library(tictoc)
library(grid)
library(lmtp)
library(sl3)
options(sl3.verbose = FALSE)

# load utils
source('./src/tmle_calculation.R')
source('./src/tmleContrast.R')
source('./src/misc_fns.R')

# output directory
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

# command line args
args <- commandArgs(trailingOnly = TRUE) # command line arguments
outcome <- as.character(args[1]) # "combined" (death or diabetes), "diabetes", or "death"
outcome.type <-  as.character(args[2]) # binomial" 
condition <-  as.character(args[3]) # "none", "schizophrenia", or "mdd"
use.SL <- as.logical(args[4])  # When TRUE, use Super Learner for initial Y model and treatment model estimation; if FALSE, use GLM

scale.continuous <- TRUE # standardize continuous covariates

# load data
load("simdata_from_basevars.RData")

simdata_from_basevars <- simdata_from_basevars[simdata_from_basevars$month_number==0,] # keep index month

# create outcome

Y.combined <- ifelse((!is.na(simdata_from_basevars$days_to_death) & simdata_from_basevars$days_to_death <=1095) | (!is.na(simdata_from_basevars$days_to_diabetes) & simdata_from_basevars$days_to_diabetes <= 1095), 1, 0)
Y.death <- ifelse((!is.na(simdata_from_basevars$days_to_death) & simdata_from_basevars$days_to_death <=1095), 1, 0)
Y.diabetes <- ifelse((!is.na(simdata_from_basevars$days_to_diabetes) & simdata_from_basevars$days_to_diabetes <= 1095), 1, 0)

if(outcome=='combined'){ # outcome is death or diabetes within 36 months
  Y <- Y.combined
}

if(outcome=='death'){ # outcome is death within 36 months
  Y <- Y.death
}

if(outcome=='diabetes'){ # outcome is diabetes within 36 months
  Y <- Y.diabetes
}

# treatment

A <- factor(simdata_from_basevars$drug_group) # categorical
obs.treatment <- dummify(A) # dummy matrix

# covariates

L.unscaled <- cbind(dummify(factor(simdata_from_basevars$state_character),show.na = FALSE)[,-c(6)],
                    dummify(factor(simdata_from_basevars$race_defined_character),show.na = FALSE)[,-c(3)], 
                    dummify(factor(simdata_from_basevars$smi_condition))[,-c(1)], # omit bipolar
                    "year"=simdata_from_basevars$year, 
                    "female"=simdata_from_basevars$female, 
                    "payer_index_mdcr"=dummify(factor(simdata_from_basevars$payer_index))[,2],
                    "preperiod_ever_psych"=simdata_from_basevars$preperiod_ever_psych,
                    "preperiod_ever_metabolic"=simdata_from_basevars$preperiod_ever_metabolic,
                    "preperiod_ever_other"=simdata_from_basevars$preperiod_ever_other,
                    "preperiod_drug_use_days"=simdata_from_basevars$preperiod_drug_use_days,
                    "preperiod_ever_mt_gluc_or_lip"=simdata_from_basevars$preperiod_ever_mt_gluc_or_lip,
                    "preperiod_ever_rx_antidiab"=simdata_from_basevars$preperiod_ever_rx_antidiab,
                    "preperiod_ever_rx_cm_nondiab"=simdata_from_basevars$preperiod_ever_rx_cm_nondiab,
                    "preperiod_ever_rx_other"=simdata_from_basevars$preperiod_ever_rx_other,
                    "calculated_age"=simdata_from_basevars$calculated_age,
                    "preperiod_olanzeq_dose_total"=simdata_from_basevars$preperiod_olanzeq_dose_total,
                    "preperiod_er_mhsa"=simdata_from_basevars$preperiod_er_mhsa,
                    "preperiod_er_nonmhsa"=simdata_from_basevars$preperiod_er_nonmhsa,
                    "preperiod_er_injury"=simdata_from_basevars$preperiod_er_injury,
                    "preperiod_cond_mhsa"=simdata_from_basevars$preperiod_cond_mhsa,
                    "preperiod_cond_nonmhsa"=simdata_from_basevars$preperiod_cond_nonmhsa,
                    "preperiod_cond_injury"=simdata_from_basevars$preperiod_cond_injury,
                    "preperiod_los_mhsa"=simdata_from_basevars$preperiod_los_mhsa,
                    "preperiod_los_nonmhsa"=simdata_from_basevars$preperiod_los_nonmhsa,
                    "preperiod_los_injury"=simdata_from_basevars$preperiod_los_injury) 

colnames(L.unscaled)[which(colnames(L.unscaled)=="West Virginia")] <- "West_Virginia"

L <- L.unscaled

if(scale.continuous){
  continuous.vars <- c("calculated_age","preperiod_drug_use_days","preperiod_olanzeq_dose_total","preperiod_er_mhsa","preperiod_er_nonmhsa","preperiod_er_injury","preperiod_cond_mhsa",
                       "preperiod_cond_nonmhsa","preperiod_cond_injury","preperiod_los_mhsa","preperiod_los_nonmhsa","preperiod_los_injury")
  
  L[,continuous.vars] <- scale(L[,continuous.vars]) # scale continuous vars
}

## Summary statistics table

if(condition=="none"){
  summary.tables <- TRUE
  if(summary.tables){
    print(tableNominal(data.frame(L.unscaled,A,Y.combined,Y.death, Y.diabetes)[c( "California","Georgia","Iowa",                         
                                                                         "Mississippi","Oklahoma","West_Virginia","black","latino","white","mdd","schiz","year","female",
                                                                         "payer_index_mdcr","preperiod_ever_psych","preperiod_ever_metabolic","preperiod_ever_other","preperiod_ever_mt_gluc_or_lip",
                                                                         "preperiod_ever_rx_antidiab", "preperiod_ever_rx_cm_nondiab", "preperiod_ever_rx_other","Y.combined","Y.death", "Y.diabetes")], group=A, prec = 3, cumsum=FALSE, longtable = FALSE))
    
    print(tableContinuous(data.frame(L.unscaled,A)[c("calculated_age","preperiod_drug_use_days","preperiod_olanzeq_dose_total","preperiod_er_mhsa","preperiod_er_nonmhsa","preperiod_er_injury","preperiod_cond_mhsa",
                                                     "preperiod_cond_nonmhsa","preperiod_cond_injury","preperiod_los_mhsa","preperiod_los_nonmhsa","preperiod_los_injury")], group=A, prec = 3, cumsum=FALSE, stats= c("n", "min", "mean", "max", "s"), longtable = FALSE))
  }
}

rm(simdata_from_basevars)

## create condition dummy

if(condition=="schizophrenia"){
  x <- ifelse(L[,'schiz']==1,1,0)
} else if(condition=="mdd"){
  x <- ifelse(L[,'mdd']==1,1,0)
}else if(condition=="black"){
  x <- ifelse(L[,'black']==1,1,0)
}else if(condition=="latino"){
  x <- ifelse(L[,'latino']==1,1,0)
}else if(condition=="white"){
  x <- ifelse(L[,'white']==1,1,0)
}

## set parameters

J <- 6
n <- length(Y) # sample size

gbound <- c(0.001,0.999) # define bounds to be used for the propensity score

ybound <- c(0.0001, 0.9999) # define bounds to be used for the Y predictions

n.folds <- 5 # number of folds for SL

target.gwt <- TRUE # When TRUE, moves propensity weights from denominator of clever covariate to regression weight when fitting updated model for Y

est.binomial <- TRUE # When True, estimates binomial treatment model
est.LMTP <- FALSE # When True, estimates LMTP

# stack learners into a model
learner_stack_A <- make_learner_stack(list("Lrnr_ranger",num.trees=50),list("Lrnr_ranger",num.trees=100),list("Lrnr_ranger",num.trees=500), list("Lrnr_glmnet",nfolds = n.folds,alpha = 1, family = "multinomial"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0, family = "multinomial"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.25, family = "multinomial"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.5, family = "multinomial"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.75, family = "multinomial"))  
learner_stack_A_bin <- make_learner_stack(list("Lrnr_ranger",num.trees=50),list("Lrnr_ranger",num.trees=100),list("Lrnr_ranger",num.trees=500), list("Lrnr_glmnet",nfolds = n.folds,alpha = 1, family = "binomial"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0, family = "binomial"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.25, family = "binomial"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.5, family = "binomial"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.75, family = "binomial"))  

# metalearner defaults 
if(outcome.type=="continuous"){
  metalearner_Y <- make_learner(Lrnr_solnp,learner_function=metalearner_linear,eval_function=loss_squared_error) # nonlinear optimization via augmented Lagrange
  learner_stack_Y <- make_learner_stack(list("Lrnr_xgboost",nrounds=20, objective = "reg:linear"),list("Lrnr_ranger",num.trees=50), list("Lrnr_ranger",num.trees=100), list("Lrnr_ranger",num.trees=500), list("Lrnr_glmnet",nfolds = n.folds,alpha = 1, family = "gaussian"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0, family = "gaussian"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.25, family = "gaussian"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.5, family = "gaussian"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.75, family = "gaussian")) 
}else if (outcome.type =="binomial"){
  metalearner_Y <- make_learner(Lrnr_solnp,learner_function=metalearner_logistic_binomial,eval_function=loss_loglik_binomial)
  learner_stack_Y <- make_learner_stack(list("Lrnr_xgboost",nrounds=20, objective = "reg:logistic"),list("Lrnr_ranger",num.trees=50), list("Lrnr_ranger",num.trees=100), list("Lrnr_ranger",num.trees=500), list("Lrnr_glmnet",nfolds = n.folds,alpha = 1, family = "binomial"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0, family = "binomial"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.25, family = "binomial"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.5, family = "binomial"), list("Lrnr_glmnet",nfolds = n.folds,alpha = 0.75, family = "binomial")) 
}
metalearner_A <- make_learner(Lrnr_solnp,learner_function=metalearner_linear_multinomial, eval_function=loss_loglik_multinomial)
metalearner_A_bin <- make_learner(Lrnr_solnp,learner_function=metalearner_logistic_binomial,eval_function=loss_loglik_binomial)

#create contrast matrix (1 is treated, -1s are reference)
loc <- t(combn(J,2))
colnames(loc) <- c("reference","treated")

C <- matrix(0, nrow(loc), J)
for(i in 1:nrow(C)){
  C[i,loc[i,1]] = -1
  C[i,loc[i,2]] = 1
}
colnames(C) <-1:J

C <- C[which(C[,1]==-1),] # Aripripazole as reference

m <- nrow(C) # number of comparisons

## Manual TMLE (ours)

tic(paste0("Initialize manual TMLE","outcome = ", outcome, "; outcome.type = ", outcome.type, "condition = ", condition, "est.binomial = ", est.binomial, "use.SL =",use.SL))

## fit initial outcome model
if(use.SL){
  
  # define task and candidate learners
  initial_model_for_Y_task <- make_sl3_Task(data.frame(Y,L,obs.treatment[,-2]), covariates = c(colnames(L),colnames(obs.treatment[,-2])), outcome = "Y", outcome_type="binomial", 
                                            folds = origami::make_folds(data.frame(Y,L,obs.treatment[,-2]), fold_fun = folds_vfold, V = n.folds)) # HALOPERIDOL is reference
  
  # Super Learner algorithm
  
  initial_model_for_Y_sl <- make_learner(Lrnr_sl, # cross-validates base models
                                         learners = learner_stack_Y,
                                         metalearner = metalearner_Y)
  initial_model_for_Y_sl_fit <- initial_model_for_Y_sl$train(initial_model_for_Y_task)
  initial_model_for_Y_preds <- initial_model_for_Y_sl_fit$predict(initial_model_for_Y_task) # predicted probs.
  
  if(condition=="none"){
    saveRDS(initial_model_for_Y_sl_fit,paste0(output_dir,"tmle_",outcome,"_",outcome.type, "_", condition, "_", est.binomial, "_", est.LMTP, "_","initial_model_for_Y_sl_fit.rds")) # save model
  }
  
  if(outcome=="combined" & condition=="none" & outcome.type =="binomial"){
    Y.risks <- initial_model_for_Y_sl_fit$cv_risk(loss_loglik_binomial)
    print("Outcome model - SL CV risk")
    print(Y.risks)
  }
  
  for(j in 1:J){
    newdata <- data.frame(Y,L,matrix(0, nrow=nrow(obs.treatment), ncol=ncol(obs.treatment), dimnames = dimnames(obs.treatment)))
    newdata[paste0("D",j)] <- 1
    assign(paste0("Q",j), initial_model_for_Y_sl_fit$predict(sl3_Task$new(newdata, covariates = c(colnames(L),colnames(obs.treatment)), outcome = "Y", outcome_type="binomial")))
  }
} else{
  if(outcome.type=="continuous"){
    initial_model_for_Y <- glm(Y~.,data=data.frame(cbind(Y,obs.treatment[,-2],L)), family=gaussian(), control = glm.control(maxit = 100)) 
    initial_model_for_Y_preds <- predict(initial_model_for_Y)
  }else{
    initial_model_for_Y <- glm(Y~.,data=data.frame(cbind(Y,obs.treatment[,-2],L)), family=binomial(), control = glm.control(maxit = 100)) 
    initial_model_for_Y_preds <- predict(initial_model_for_Y,type="response") # predicted probs.
  }
  
  if(condition=="none"){
    saveRDS(initial_model_for_Y,paste0(output_dir,"tmle_",outcome,"_",outcome.type, "_", condition, "_", est.binomial, "_", est.LMTP, "_","initial_model_for_Y.rds")) # save model
  }
  
  for(j in 1:J){
    newdata <- data.frame(cbind(Y,L,matrix(0, nrow=nrow(obs.treatment), ncol=ncol(obs.treatment), dimnames = dimnames(obs.treatment))))
    newdata[paste0("D",j)] <- 1
    assign(paste0("Q",j), predict(initial_model_for_Y,newdata=data.frame(newdata), type="response"))
  }
}
rm(newdata)

Qs <- t(do.call(mapply, c(FUN = cbind, mget(paste0("Q", 1:J)))))
colnames(Qs) <- paste0("Q", 1:J)

QAW <- data.frame(apply(cbind(QA=initial_model_for_Y_preds,Qs), 2, boundProbs, bounds=ybound)) # bound predictions

##  fit initial treatment model

if(use.SL){
  
  # multinomial
  
  initial_model_for_A_task <- make_sl3_Task(data.frame(A,L), covariates = colnames(L), outcome = "A", outcome_type="categorical", 
                                            folds = origami::make_folds(data.frame(A,L), fold_fun = folds_vfold, V = n.folds))
  
  initial_model_for_A_sl <- make_learner(Lrnr_sl, # cross-validates base models
                                         learners = learner_stack_A,
                                         metalearner = metalearner_A)
  initial_model_for_A_sl_fit <- initial_model_for_A_sl$train(initial_model_for_A_task)
  
  g_preds  <- initial_model_for_A_sl_fit$predict(initial_model_for_A_task)  # estimated propensity scores 
  
  g_preds  <-  data.frame(matrix(unlist(lapply(g_preds, unlist)), nrow=length(lapply(g_preds, unlist)), byrow=TRUE)) 
  
  colnames(g_preds) <- colnames(obs.treatment)
  
  if(condition=="none"){
    saveRDS(initial_model_for_A_sl_fit,paste0(output_dir,"tmle_",outcome,"_",outcome.type, "_", condition, "_", est.binomial, "_", est.LMTP, "_","initial_model_for_A_sl_fit.rds")) # save model
  }
  
} else{
  
  # multinomial logistic regression
  initial_model_for_A <- nnet::multinom(formula=A ~ ., data=data.frame(cbind(A, L)), maxit = 500, trace = FALSE, model = TRUE) 
  g_preds <- fitted(initial_model_for_A)
  colnames(g_preds) <- colnames(obs.treatment)
  
  if(condition=="none"){
    saveRDS(initial_model_for_A,paste0(output_dir,"tmle_",outcome,"_",outcome.type, "_", condition, "_", est.binomial, "_", est.LMTP, "_","initial_model_for_A.rds")) # save model
  }
  
  # binomial logistic regression
  initial_model_for_A_bin <- lapply(1:J, function(j) glm(formula=treatment ~ ., data=data.frame(cbind("treatment"=obs.treatment[,j], L)), family="binomial"))
  g_preds_bin <- sapply(1:J, function(j) predict(initial_model_for_A_bin[[j]], type="response"))
  colnames(g_preds_bin) <- colnames(obs.treatment)
  
  if(condition=="none"){
    saveRDS(initial_model_for_A_bin,paste0(output_dir,"tmle_",outcome,"_",outcome.type, "_", condition, "_", est.binomial, "_", est.LMTP, "_","initial_model_for_A_bin.rds")) # save model
  }
}

if(outcome=="combined" & condition=="none" & use.SL){
  A.risks <- initial_model_for_A_sl_fit$cv_risk(loss_loglik_multinomial)
  print("Treatment model - multinomial - SL CV risk")
  print(A.risks)
}

Ahat_tmle_bin  <- c()
yhat_tmle_bin <- c()
ess_ate_tmle_bin <- c()
ATE_tmle_bin  <-  c()
ATE_CI_tmle_bin_upper  <-  c()
ATE_CI_tmle_bin_lower <-  c()
CATE_tmle_bin  <- c()
CATE_CI_tmle_bin_upper  <- c()
CATE_CI_tmle_bin_lower <- c()

if(est.binomial & use.SL){
  # binomial 
  
  tic(paste0("TMLE multiple binary treatment","outcome = ", outcome, "; outcome.type = ", outcome.type, "condition = ", condition))
  
  initial_model_for_A_task_bin <- lapply(1:J, function(j) make_sl3_Task(data.frame("A"=obs.treatment[,j],L), covariates = colnames(L), outcome = "A", outcome_type="binomial",
                                                                        folds = origami::make_folds(data.frame(A,L), fold_fun = folds_vfold, V = n.folds)))
  
  initial_model_for_A_sl_bin <- make_learner(Lrnr_sl, # cross-validates base models
                                             learners = learner_stack_A_bin,
                                             metalearner = metalearner_A_bin)
  initial_model_for_A_sl_fit_bin <- lapply(1:J, function(j) initial_model_for_A_sl_bin$train(initial_model_for_A_task_bin[[j]]))
  
  g_preds_bin  <- sapply(1:J, function(j) initial_model_for_A_sl_fit_bin[[j]]$predict(initial_model_for_A_task_bin[[j]]))  # estimated propensity scores 
  
  colnames(g_preds_bin) <- colnames(obs.treatment)
  
  if(condition=="none"){
    saveRDS(initial_model_for_A_sl_fit_bin,paste0(output_dir,"tmle_",outcome,"_",outcome.type, "_", condition, "_", est.binomial, "_", est.LMTP, "_","initial_model_for_A_sl_fit_bin.rds")) # save model
  }
  
  # perform checks on propensity scores
  
  if(any((g_preds_bin < 0) || (g_preds_bin > 1))){ #     # verify they are between 0 and 1
    stop("Binomial treatment probablities must be between 0 and 1")
  }
  
  g_preds_bin_bounded <- boundProbs(g_preds_bin, gbound) # (trimmed) estimated propensity scores
  
  # Update equations, calculate ATE and variance
  TMLE_bin <- getTMLE(QAW, g_preds_bin_bounded, obs.treatment, Y, outcome.type, target.gwt)
  
  # calculate contrasts
  if(condition=="none"){
    tmle_contrasts_bin <- tmle_contrast(Qstar=TMLE_bin, C, g_preds_bin_bounded, obs.treatment, Y, outcome.type, alpha=0.05, multi.adjust=FALSE)
  }else{
    tmle_contrasts_bin <- tmle_contrast(Qstar=TMLE_bin, C, g_preds_bin_bounded, obs.treatment, Y, outcome.type, alpha=0.05, multi.adjust=FALSE, x)
  }
  
  # store results
  Ahat_tmle_bin  <- g_preds_bin
  yhat_tmle_bin <- colMeans(tmle_contrasts_bin$Qstar)
  
  ess_ate_tmle_bin <- (colSums(obs.treatment*(1/g_preds_bin))**2)/(colSums(obs.treatment*(1/g_preds_bin)**2)) 
  
  ATE_tmle_bin  <- tmle_contrasts_bin$taus$ATE # point and variance estimates for ATE (binomial)
  ATE_CI_tmle_bin_upper  <- sapply(tmle_contrasts_bin$CI,"[[",1)[2,]
  ATE_CI_tmle_bin_lower <- sapply(tmle_contrasts_bin$CI,"[[",1)[1,]
  
  if(condition!="none"){
    CATE_tmle_bin  <- tmle_contrasts_bin$taus$CATE # point and variance estimates for CATE (binomial)
    CATE_CI_tmle_bin_upper  <- sapply(tmle_contrasts_bin$CI,"[[",2)[2,]
    CATE_CI_tmle_bin_lower <- sapply(tmle_contrasts_bin$CI,"[[",2)[1,] 
  }else if (condition=="none"){
    CATE_tmle_bin  <- NULL
    CATE_CI_tmle_bin_upper  <- NULL
    CATE_CI_tmle_bin_lower <- NULL
  }
}

# perform checks on propensity scores

if(sum(rowSums(g_preds))!=n){ #     # verify they sum to one
  print("Multinomial treatment probablities do not sum to 1")
}

if(any((g_preds < 0) || (g_preds > 1))){ #     # verify they are between 0 and 1
  stop("Multinomial treatment probablities must be between 0 and 1")
}

g_preds_bounded <- boundProbs(g_preds, gbound) # (trimmed) estimated propensity scores

# Update equations, calculate ATE and variance
TMLE <- getTMLE(QAW, g_preds_bounded, obs.treatment, Y, outcome.type, target.gwt)

# calculate contrasts
if(condition=="none"){
  tmle_contrasts <- tmle_contrast(Qstar=TMLE, C, g_preds_bounded, obs.treatment, Y, outcome.type, alpha=0.05, multi.adjust=FALSE)
}else{
  tmle_contrasts <- tmle_contrast(Qstar=TMLE, C, g_preds_bounded, obs.treatment, Y, outcome.type, alpha=0.05, multi.adjust=FALSE, x)
}

# store results

yinitial_tmle  <- QAW[-1] # estimated initial Y, treatment probs. 
Ahat_tmle  <- g_preds

yhat_tmle <- colMeans(tmle_contrasts$Qstar)

ess_ate_tmle <- (colSums(obs.treatment*(1/g_preds))**2)/(colSums(obs.treatment*(1/g_preds)**2)) 

ATE_tmle  <- tmle_contrasts$taus$ATE# point and variance estimates for ATE (multinomial)
ATE_CI_tmle_upper <- sapply(tmle_contrasts$CI,"[[",1)[2,]
ATE_CI_tmle_lower <- sapply(tmle_contrasts$CI,"[[",1)[1,]

if(condition!="none"){
  CATE_tmle <- tmle_contrasts$taus$CATE # point and variance estimates for CATE (binomial)
  CATE_CI_tmle_upper <- sapply(tmle_contrasts$CI,"[[",2)[2,]
  CATE_CI_tmle_lower <- sapply(tmle_contrasts$CI,"[[",2)[1,] 
}else if (condition=="none"){
  CATE_tmle  <- NULL
  CATE_CI_tmle_upper <- NULL
  CATE_CI_tmle_lower <- NULL
}

print(toc())

if(est.LMTP){

  lmtp_results <- list()
  yinitial_lmtp <- list()
  Ahat_lmtp <- list()
  ess_ate_lmtp <- list()
  
  results_ATE_lmtp <- list()
  results_ate_CI_lmtp <- list()
  
  # estimate outcomes under each J
  for(j in c(1:J)){
    lmtp_results[[j]] <- lmtp_tmle(data = data.frame("Y"=Y,"A"=A,L), 
                                   trt = "A",
                                   outcome = "Y", 
                                   baseline = colnames(L), 
                                   shifted = data.frame("Y"=Y,"A"=factor(rep(levels(A)[j],times=length(A))),L),
                                   intervention_type = "mtp",
                                   outcome_type = outcome.type, 
                                   folds = n.folds,
                                   .bound = if(outcome.type=="binomial") ybound[1] else 1e-5,
                                   .trim = gbound[2],
                                   .SL_folds = n.folds,
                                   learners_trt= learner_stack_A_bin,
                                   learners_outcome= learner_stack_Y)
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
    results_ATE_lmtp[[i]]  <- lmtp_contrasts[[i]]$vals$theta # point and variance estimates
    results_ate_CI_lmtp[[i]]  <- cbind("low"=lmtp_contrasts[[i]]$vals$conf.low, 
                                       "high"=lmtp_contrasts[[i]]$vals$conf.high)
  }
}

# save results
if(est.binomial==TRUE && est.LMTP==TRUE){
  saveRDS(list("yinitial_tmle"=yinitial_tmle,"Ahat_tmle"=Ahat_tmle,"yhat_tmle"= yhat_tmle, "ess_ate_tmle"=ess_ate_tmle,
               "Ahat_tmle_bin"=Ahat_tmle_bin,"yhat_tmle_bin"= yhat_tmle_bin, "ess_ate_tmle_bin"=ess_ate_tmle_bin,
               "ATE_tmle"= ATE_tmle,"ATE_CI_tmle_upper"=ATE_CI_tmle_upper, "ATE_CI_tmle_lower"=ATE_CI_tmle_lower,
               "CATE_tmle"= CATE_tmle,"CATE_CI_tmle_upper"=CATE_CI_tmle_upper, "CATE_CI_tmle_lower"=CATE_CI_tmle_lower,
               "ATE_tmle_bin"= ATE_tmle_bin,"ATE_CI_tmle_bin_upper"=ATE_CI_tmle_bin_upper, "ATE_CI_tmle_bin_lower"=ATE_CI_tmle_bin_lower,
               "CATE_tmle_bin"= CATE_tmle_bin,"CATE_CI_tmle_bin_upper"=CATE_CI_tmle_bin_upper, "CATE_CI_tmle_bin_lower"=CATE_CI_tmle_bin_lower,
               "yinitial_lmtp"=unlist(yinitial_lmtp), "Ahat_lmtp"=unlist(Ahat_lmtp), "ess_ate_lmtp"=unlist(ess_ate_lmtp),
               "ATE_lmtp"= unlist(results_ATE_lmtp),"ATE_CI_lmtp_upper"=unlist(lapply(results_ate_CI_lmtp, '[[', 2)), "ATE_CI_lmtp_lower"=unlist(lapply(results_ate_CI_lmtp, '[[', 1))),paste0(output_dir,"tmle_",outcome,"_",outcome.type, "_", condition, "_", est.binomial, "_", est.LMTP, "_","itt_analysis.rds"))
}else{
  saveRDS(list("yinitial_tmle"=yinitial_tmle,"Ahat_tmle"=Ahat_tmle,"yhat_tmle"= yhat_tmle, "ess_ate_tmle"=ess_ate_tmle,
               "Ahat_tmle_bin"=Ahat_tmle_bin,"yhat_tmle_bin"= yhat_tmle_bin, "ess_ate_tmle_bin"=ess_ate_tmle_bin,
               "ATE_tmle"= ATE_tmle,"ATE_CI_tmle_upper"=ATE_CI_tmle_upper, "ATE_CI_tmle_lower"=ATE_CI_tmle_lower,
               "CATE_tmle"= CATE_tmle,"CATE_CI_tmle_upper"=CATE_CI_tmle_upper, "CATE_CI_tmle_lower"=CATE_CI_tmle_lower,
               "ATE_tmle_bin"= ATE_tmle_bin,"ATE_CI_tmle_bin_upper"=ATE_CI_tmle_bin_upper, "ATE_CI_tmle_bin_lower"=ATE_CI_tmle_bin_lower,
               "CATE_tmle_bin"= CATE_tmle_bin,"CATE_CI_tmle_bin_upper"=CATE_CI_tmle_bin_upper, "CATE_CI_tmle_bin_lower"=CATE_CI_tmle_bin_lower),paste0(output_dir,"tmle_",outcome,"_",outcome.type, "_", condition, "_", est.binomial, "_","itt_analysis.rds"))
}

## Summary stats and results plots
source('./tmle_itt_analysis_eda.R')