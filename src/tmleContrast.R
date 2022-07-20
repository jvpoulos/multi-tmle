# function for calculating contrasts
tmle_contrast <- function(Qstar,C, g_preds_bounded, obs.treatment, Y, outcome.type, alpha=0.05, multi.adjust=FALSE, x=NULL){
  # x is binary covariate to condition on for CATE
  
  if(outcome.type=="binomial"){
    MuA <- mapply(plogis,Qstar$QA)
    Mu <- mapply(plogis,Qstar[-1])
  }else{
    MuA <- Qstar$QA
    Mu <- Qstar[-1]
  }
  
  taus <- vector(mode = "list", length = 2) 
  names(taus) <- c("ATE","CATE")
  
  taus[[1]] <- as.vector(C%*%colMeans(Mu))  # ATE_(treat, ref) = E[Y_i (treat) - Y_i (ref)]
  
  if(!is.null(x)){
    taus[[2]] <- as.vector(C%*%colMeans(Mu[which(x==1),]))  # CATE_(treat, ref) = E[Y_i (treat) - Y_i (ref) | x]
  }
  
  infcurv <- vector(mode = "list", length = nrow(C)) # structure is [[comparison]][[estimand]]
  var <- vector(mode = "list", length = nrow(C))
  CI <- vector(mode = "list", length = nrow(C))
  
  # loop through each pairwise comparison
  for(i in 1:nrow(C)){
    treat_index <- as.numeric(colnames(C)[which(C[i,]==1)]) # 1 is treated, -1s are reference
    ref_index <- as.numeric(colnames(C)[which(C[i,]==-1)]) 
    
    # calculate influence curve
    
    infcurv[[i]] <-  lapply(1:length(taus), function(t){((obs.treatment[,treat_index]/g_preds_bounded[,treat_index])-(obs.treatment[,ref_index]/g_preds_bounded[,ref_index]))*(Y-MuA) + Mu[,treat_index] - Mu[,ref_index] - taus[[t]][i]})

    if(!is.null(x)){
      var[[i]] <- lapply(1:length(taus), function(t){var(infcurv[[i]][[t]][which(x==1)])/sum(x)})
    }else{
      var[[i]] <- lapply(1:length(taus), function(t){var(infcurv[[i]][[t]])/length(Y)})
    }
    
    if(multi.adjust){
      CI[[i]] <- lapply(1:length(taus), function(t){c(taus[[t]][i] -qnorm(1-(alpha/(2*nrow(C)))) *sqrt(var[[i]][[t]]), taus[[t]][i] +qnorm(1-(alpha/(2*nrow(C)))) *sqrt(var[[i]][[t]]))})
    }else{
      CI[[i]] <- lapply(1:length(taus), function(t){c(taus[[t]][i] -qnorm(1-(alpha/2)) *sqrt(var[[i]][[t]]), taus[[t]][i] +qnorm(1-(alpha/2)) *sqrt(var[[i]][[t]]))})
    }
  }
  
  return(list("taus"=taus,
              "QA"=MuA,
              "Qstar"=Mu,
              "IC"=infcurv,
              "var" = var, 
              "CI"=CI))
}