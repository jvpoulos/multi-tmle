####################################################################
# Censoring analysis (ITT) for initial cohort of n=64120 patients  #
####################################################################

# Load required libraries
library(survival)
library(survminer)
library(weights)

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

use.simulated <- FALSE
# load data
if(use.simulated){
  load("simdata_from_basevars.RData")
  
  # make censoring variable name consistent with CMS dataset
  names(simdata_from_basevars)[names(simdata_from_basevars) == 'days_to_censored_or_outcome'] <- 'days_to_censored'
  
  fup3yr_episode_months_deid <- simdata_from_basevars
}else{
  load("/data/MedicaidAP_associate/poulos/fup3yr_episode_months_deid_fixmonthout.RData") # run on Argos
}


# covariates
L.unscaled <- cbind(dummify(factor(fup3yr_episode_months_deid $state_character),show.na = FALSE)[,-c(6)],
                    dummify(factor(fup3yr_episode_months_deid $race_defined_character),show.na = FALSE)[,-c(3)],
                    dummify(factor(fup3yr_episode_months_deid $smi_condition))[,-c(1)], # omit bipolar
                    "payer_index_mdcr"=dummify(factor(fup3yr_episode_months_deid $payer_index))[,2],
                    "year"=fup3yr_episode_months_deid $year,
                    "female"=fup3yr_episode_months_deid $female,
                    "preperiod_ever_psych"=fup3yr_episode_months_deid $preperiod_ever_psych,
                    "preperiod_ever_metabolic"=fup3yr_episode_months_deid $preperiod_ever_metabolic,
                    "preperiod_ever_other"=fup3yr_episode_months_deid $preperiod_ever_other,
                    "preperiod_drug_use_days"=fup3yr_episode_months_deid $preperiod_drug_use_days,
                    "preperiod_ever_mt_gluc_or_lip"=fup3yr_episode_months_deid $preperiod_ever_mt_gluc_or_lip,
                    "preperiod_ever_rx_antidiab"=fup3yr_episode_months_deid $preperiod_ever_rx_antidiab,
                    "preperiod_ever_rx_cm_nondiab"=fup3yr_episode_months_deid $preperiod_ever_rx_cm_nondiab,
                    "preperiod_ever_rx_other"=fup3yr_episode_months_deid $preperiod_ever_rx_other,
                    "calculated_age"=fup3yr_episode_months_deid $calculated_age,
                    "preperiod_er_mhsa"=fup3yr_episode_months_deid $preperiod_er_mhsa,
                    "preperiod_er_nonmhsa"=fup3yr_episode_months_deid $preperiod_er_nonmhsa,
                    "preperiod_er_injury"=fup3yr_episode_months_deid $preperiod_er_injury,
                    "preperiod_cond_mhsa"=fup3yr_episode_months_deid $preperiod_cond_mhsa,
                    "preperiod_cond_nonmhsa"=fup3yr_episode_months_deid $preperiod_cond_nonmhsa,
                    "preperiod_cond_injury"=fup3yr_episode_months_deid $preperiod_cond_injury,
                    "preperiod_los_mhsa"=fup3yr_episode_months_deid $preperiod_los_mhsa,
                    "preperiod_los_nonmhsa"=fup3yr_episode_months_deid $preperiod_los_nonmhsa,
                    "preperiod_los_injury"=fup3yr_episode_months_deid $preperiod_los_injury)

colnames(L.unscaled)[which(colnames(L.unscaled)=="West Virginia")] <- "West_Virginia"

L <- L.unscaled

scale.continuous <- TRUE
if(scale.continuous){
  continuous.vars <- c("calculated_age","preperiod_drug_use_days","preperiod_er_mhsa","preperiod_er_nonmhsa","preperiod_er_injury","preperiod_cond_mhsa",
                       "preperiod_cond_nonmhsa","preperiod_cond_injury","preperiod_los_mhsa","preperiod_los_nonmhsa","preperiod_los_injury")

  L[,continuous.vars] <- scale(L[,continuous.vars]) # scale continuous vars
}

use.SL <- FALSE
if(use.SL==FALSE){ # exclude multicolinear variables for GLM (VIF)
  glm.exlude <- c("California")
  L <- L[,!colnames(L)%in%glm.exlude]
}

# Censoring and median follow-up time (Table A2)

A <- factor(fup3yr_episode_months_deid$drug_group) # categorical

C <- dummify(factor(fup3yr_episode_months_deid$censoring_cause_4cats)) # Overall cause of censoring - 4 categories,

days <- fup3yr_episode_months_deid$days_to_end_of_followup # Days to end of followup (min (dod_both_sources,censored_date,days_to_index_plus_3yrs))

print(aggregate(C, by=list(A), FUN=mean))

print(aggregate(days, by=list(A), FUN=mean))

##Kaplan-Meier Survival Curves
## plot Kaplan-Meier survival curves for each treatment group

fup3yr_episode_months_deid$month_number <- fup3yr_episode_months_deid$month_number+1 # baseline starts at t=1

# Censoring dummy
fup3yr_episode_months_deid$censor <- ifelse(fup3yr_episode_months_deid$days_to_censored <= (1095/max(fup3yr_episode_months_deid$month_number))*fup3yr_episode_months_deid$month_number,1,0)

# Create the survival object
surv_obj <- Surv(time = fup3yr_episode_months_deid$month_number, event = fup3yr_episode_months_deid$censor)

# Plot Kaplan-Meier survival curves
g <- ggsurvplot(
  fit = survfit(surv_obj ~ factor(fup3yr_episode_months_deid$drug_group)),
  data = fup3yr_episode_months_deid,
  risk.table = "abs_pct",
  palette="lancet",
  combine = TRUE,
  ylab= "Censoring probability",
  xlab= "Month",
  legend="top",
  legend.title="Antipsychotic",
  legend.labs=c("Reference","A","B","C","D","E"),
  surv.scale="percent",
  xlim=c(1,37),
  tables.y.text=TRUE,
  risk.table.title="Number of patients at risk of censoring (%)",
  tables.theme = theme_cleantable()
)

png(paste0(output_dir,"censoring_KM_plot.png"),width = 780, height = 780)
print(g) 
dev.off() # Close the file

## Cox Proportional Hazards Models
## assess the impact of treatment assignment and other covariates on the hazard of being censored

# Fit the Cox proportional hazards model
cox_model <- coxph(surv_obj ~ fup3yr_episode_months_deid$drug_group + L)

# Summary of the model
summary(cox_model)

# Test the assumption of proportional hazards
cox_zph <- cox.zph(cox_model)
print(cox_zph)

# Plot the scaled Schoenfeld residuals to visualize the assumption
png(paste0(output_dir,"schoenfeld_plot.png"),width = 780, height = 780)
plot(cox_zph, xlab="Month")
dev.off() # Close the file