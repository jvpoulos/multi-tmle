################################################
# Functions for generating simulation data     #
################################################

## data generation function for static simulation

generateData <- function(r, J, n, overlap.setting, gamma.setting, outcome.type, scale.continuous){
  # Code adapted from Li and Li (2019, AOAS) 
  # The simulation setup is similar to Yang et al. (2016, Biometrics)
  
  # X1-X3 are multivariate normal covariates
  # X4 is uniform variable
  # X5 is a Chi-squared variable
  # X6 is binary
  
  # set the seed
  set.seed(r, "L'Ecuyer-CMRG")
  
  vars <- c(2,1,1) # yang et al. values
  covars <- c(1,-1,-.5)
  mu <- c(0,0,0)
  
  Sigma <- diag(vars)
  Sigma[2,1] <- Sigma[1,2] <- covars[1]
  Sigma[3,1] <- Sigma[1,3] <- covars[2]
  Sigma[3,2] <- Sigma[2,3] <- covars[3]
  
  X13 <- mvrnorm(n, mu=mu, Sigma=Sigma, empirical = FALSE)
  X4 <- runif(n,-3,3)
  X5 <- rchisq(n, df=1)
  X6 <- rbinom(n, size=1, prob=.5)
  X16 <- cbind(X13, X4, X5, X6) # covariate set w.o. intercept
  X06 <- cbind(1, X13, X4, X5, X6) # covariate set w intercept
  
  # Assignment mechanism 
  
  if(overlap.setting=="adequate"){ # li and li values
    if(J==6){
      kappa <- c(0,0.1, 0.15, 0.2, 0.25, 0.3)
    } else if(J==3){
      kappa <- c(0,0.2, 0.1)
    }
  }else if(overlap.setting=="inadequate"){ # yang et al. values
    if(J==6){
      kappa <- c(0,0.4, 0.6, 0.8, 1.0, 1.2)
    } else if(J==3){
      kappa <- c(0,0.7, 0.4)
    }
  }else if(overlap.setting=="rct"){
    kappa <- rep(0, J)
  }

  ##########################
  # initialize tables
  ##########################
  beta_mat <- matrix(0, nrow=7, ncol=J)
  xb_mat <- matrix(0, nrow=n, ncol=J)
  exp_xb_mat <- xb_mat
  e_mat <- xb_mat
  e_denom <- 0
  
  # Simulate observed treatment 
  D <- matrix(NA, n, J) # store obs. treatment
  colnames(D) <- paste0("D", rep(1:J))
  
  w_mat <- matrix(0, nrow=7, ncol=J)   # yang et al. values
  if(J==6){
    w_mat[,1] <- c(0, 0, 0, 0, 0, 0, 0)
    w_mat[,2] <- c(0, 1, 1, 2, 1, 1, 1)
    w_mat[,3] <- c(0, 1, 1, 1, 1, 1, -5)
    w_mat[,4] <- c(0, 1, 1, 1, 1, 1, 5)
    w_mat[,5] <- c(0, 1, 1, 1, -2, 1, 1)
    w_mat[,6] <- c(1.5, 6, 1, 2,-1,-1,-1)
  } else if(J==3){
    w_mat[,1] <- c(0, 0, 0, 0, 0, 0, 0)
    w_mat[,2] <- c(0, 1, 1, 1, -1, 1, 1)
    w_mat[,3] <- c(0, 1, 1, 1, 1, 1, 1)
  }
  
  for(j in 1:J){
    beta_mat[,j] <- kappa[j]*w_mat[,j]
    xb_mat[,j] <- c(X06 %*% beta_mat[,j])
    exp_xb_mat[,j] <- exp(xb_mat[,j])
    e_denom <- e_denom + exp_xb_mat[,j]
  }
  
  # calculate e_mat
  for(j in 1:J){
    e_mat[, j] <- exp_xb_mat[,j] / e_denom
  }
  
  for(k in 1:n){
    D[k,] <- rmultinom(1, 1, prob = e_mat[k,])
  }
  
  Z <- D[,"D1"]
  for(j in 2:J){
    Z <- Z + j * D[,paste0("D",toString(j))]
  }
  
  # True potential outcome models
  
  if(gamma.setting == "yang"){ # values in yang et al. 
    if(J==6){
      gamma_mat <- matrix(data=c(-1.5, 1, 1, 1, 1, 1, 1, 
                                 -3, 2, 3, 1, 2, 2, 2, 
                                 3, 3, 1, 2, -1, -1, -4, 
                                 2.5, 4, 1, 2, -1, -1, -3, 
                                 2, 5, 1, 2, -1, -1, -2, 
                                 1.5, 6, 1, 2, -1, -1, -1),
                          nrow = 6, ncol = 7, byrow=TRUE)
    } else if(J==3){
      gamma_mat <- matrix(data=c(-1.5, 1, 1, 1, 1, 1, 1, 
                                 -3, 2, 3, 1, 2, 2, 2, 
                                 1.5, 3, 1, 2, -1, -1, -1),
                          nrow = 3, ncol = 7, byrow=TRUE)
    }

  } else if(gamma.setting == "li"){ # values in li and li
    if(J==6){
      gamma_mat <- matrix(data=c(-1.5, 1, 1, 1, 1, 1, 1, 
                                 -4, 2, 3, 1, 2, 2, 2, 
                                 4, 3, 1, 2, -1, -1, -4, 
                                 1, 4, 1, 2, -1, -1, -3, 
                                 3.5, 5, 1, 2, -1, -1, -2, 
                                 3.5, 6, 1, 2, -1, -1, -1),
                          nrow = 6, ncol = 7, byrow=TRUE)
    } else if(J==3){
      gamma_mat <- matrix(data=c(-1.5, 1, 1, 1, 1, 1, 1, 
                                 -4, 2, 3, 1, 2, 2, 2, 
                                 3, 3, 1, 2, -1, -1, -1),
                          nrow = 3, ncol = 7, byrow=TRUE)
    }

  }else if(gamma.setting == "low"){ # 
    if(J==6){
      gamma_mat <- matrix(data=c(-4, 1, -2, -1, 1, 1, 1, 
                                 -6, 1, -2, -1, 1, 1, 1, 
                                 -2, 1, -1, -1, -1, -1, -4, 
                                  1, 2, 1, 2, -1, -1, -3, 
                                 -2, 2, -1, 1, -2, -1, -3, 
                                 -3, 3, -1, 1, -2, -1, -2),
                          nrow = 6, ncol = 7, byrow=TRUE)
    } else if(J==3){
      gamma_mat <- matrix(data=c(-4, 1, -2, -1, 1, 1, 1, 
                                 -2, 1, -1, -1, -1, -1, -4, 
                                 -3, 3, -1, 1, -2, -1, -2),
                          nrow = 3, ncol = 7, byrow=TRUE)
    }
    
  }else if(gamma.setting == "zero"){
    gamma_mat <- matrix(0, nrow = J, ncol = 7, byrow=TRUE)
  } 

  EY <- list()
  
  if(outcome.type=="continuous"){
    u  <- rnorm(n)
    
    for(j in 1:J){
      EY[[j]] <- c(X06 %*% gamma_mat[j,])
    }
    
    EYs <- sapply(EY,unlist)
    Y.true <- EYs + u
    Y.obs <- rowSums(EYs*D) + u
  }else{
    Y <- list()
    for(j in 1:J){
      EY[[j]] <- exp(X06 %*% gamma_mat[j,])
      
      Y[[j]] <- as.integer(rbernoulli(n, p=EY[[j]]/(1+EY[[j]])))
      
      Y.true <- sapply(Y,unlist)
    }
    
    Y.obs <- rowSums(Y.true*D) # observed Y
  }
  
  # Store true ATEs and ATTs
  
  # create contrast matrix (1 is treated, -1s are reference)
  loc <- t(combn(J,2))
  colnames(loc) <- c("reference","treated")
  
  C <- matrix(0, nrow(loc), J)
  for(i in 1:nrow(C)){
    C[i,loc[i,1]] = -1
    C[i,loc[i,2]] = 1
  }
  colnames(C) <-1:J
  
  # estimate ATES
  true.ates <- as.vector(C%*%colMeans(Y.true)) # ATE_(treat, ref) = E[Y_i (treat) - Y_i (ref)]
  names(true.ates) <- rep(paste0("TrueATE", paste0(loc[,2]), paste0(loc[,1]))) # treated, reference


  true.atts <- c()
  for(i in 1:nrow(loc)){
    true.atts[i] <- mean(Y.true[,loc[i,][2]][which(D[,loc[i,][2]]==1)] - Y.true[,loc[i,][1]][which(D[,loc[i,][2]]==1)])  # ATT_(treat, ref) = E[Y_i (treat) - Y_i (ref) | A_i = treat]
  }
  names(true.atts) <- rep(paste0("TrueATT", paste0(loc[,2]), paste0(loc[,1]))) # treated, reference
  
  # We will now be agnostic to the DGP
  # and pretend as if we only have observed the analysis.data
  
  # Create the analysis datasets
  analysis.data <- data.frame(Y=Y.obs, Z=Z, X=X16)
  colnames(analysis.data) <- c("Y", "Z", "X1", "X2", "X3", "X4", "X5", "X6")
  
  if(scale.continuous){
    continuous.vars <- c("X1", "X2", "X3", "X4", "X5")
    
    analysis.data[,continuous.vars] <- scale(analysis.data[,continuous.vars]) # scale continuous vars
  }
  
  # return generated data
  return(list("data"=analysis.data,
              "observed.treatment"=D,
              "trueATE"=true.ates,
              "trueATT"=true.atts,
              "TrueYs"=colMeans(Y.true),
              "C"=C))
}