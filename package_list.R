packages <- c("devtools","ggplot2","nnet","tmle","MASS","purrr","tidyverse","data.table","SuperLearner","Rsolnp","reshape2","origami","tictoc","weights","grid")

super.learner <- TRUE
dependencies <- FALSE # data.table, stringi, HMisc dependencies might be needed for SuperLearner libraries
if(super.learner){
  packages <- c(packages, c("glmnet","ranger","rpart","nnls","xgboost"))
}
if(dependencies){
  remove.packages("data.table")                         # First remove the current version
  install.packages("data.table", type = "source", repos = "https://Rdatatable.gitlab.io/data.table")
  install.packages("stringi", dependencies=TRUE, INSTALL_opts = c('--no-lock'))

  install.packages(c("HMisc",packages),repos = "http://cran.us.r-project.org")
} else{
  install.packages(packages,repos = "http://cran.us.r-project.org")
}

# development packages
remotes::install_github("tlverse/sl3@v1.4.2") # Metalearner bug in v1.4.4
remotes::install_github("tlverse/tmle3")
devtools::install_github("nt-williams/lmtp@sl3")

devtools::install_github('osofr/gridisl', build_vignettes = FALSE)
devtools::install_github('osofr/stremr')
devtools::install_github('osofr/simcausal', build_vignettes = FALSE)

# doMPI
doMPI <- FALSE
if(doMPI){
  install.packages("Rmpi")
  install.packages("doMPI", dependencies=TRUE, repos = "http://cran.us.r-project.org")
}