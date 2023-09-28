################################################
# Functions for generating simulation data     #
################################################

## data generation function for static simulation

generateData <- function(r, J, n, overlap.setting, gamma.setting, outcome.type, scale.continuous,covars40,covars100){
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
  
  if(covars40 | covars100){ # additional 34 covars
    # Additional code for 34 more covariates
    X7 <- rnorm(n)
    X8 <- runif(n, 0, 10)
    X9 <- rpois(n, lambda=3)
    X10 <- rexp_truncated(n, rate=0.1)
    X11 <- rgeom(n, prob=0.5)
    X12 <- rhyper(n, m=20, n=10, k=5)
    X13_new <- rlogis(n, location=0, scale=1)
    X14 <- rweibull(n, shape=2, scale=1)
    X15 <- rt(n, df=5)
    X16_new <- rwilcox(n, m=5, n=5)
    X17 <- rcauchy_truncated(n, location=0, scale=1)
    X18 <- rf(n, df1=5, df2=5)
    X19 <- rgamma(n, shape=1, scale=1)
    X20 <- VGAM::rlaplace(n, location=0, scale=1)
    X21 <- rlnorm_truncated(n, meanlog=0, sdlog=1)
    X22 <- rnbinom(n, size=10, prob=0.5)
    X23 <- rsignrank(n, n=10)
    X24 <- rbeta(n, shape1=1, shape2=1)
    X25 <- rchisq(n, df=2)
    X26 <- rbinom(n, size=10, prob=0.5)
    X27 <- rcat(n, prob=c(0.1, 0.9))
    X28 <- rdirichlet(n, alpha=c(1,2,3))[,1]
    X29 <- rdunif(n, min=1, max=10)
    X30 <- rfactorial(n, size=5)
    X31 <- rgeom(n, prob=0.5)
    X32 <- rhyper(n, m=20, n=20, k=10)
    X33 <- rlnorm_truncated(n, meanlog=0, sdlog=1)
    X34 <- rnbinom(n, size=10, prob=0.5)
    X35 <- rpois(n, lambda=3)
    X36 <- rt(n, df=5)
    X37 <- runif(n, 0, 10)
    X38 <- rweibull(n, shape=2, scale=1)
    X39 <- rchisq(n, df=1)
    X40 <- rexp_truncated(n, rate=0.1)
    
    # Combine 40 covariates
    X_all_40 <- cbind(X13, X4, X5, X6, X7, X8, X9, X10, X11, X12, X13_new, X14, X15, X16_new, X17, X18, X19, X20, X21, X22, X23, X24, X25, X26, X27, X28, X29, X30, X31, X32, X33, X34, X35, X36, X37, X38, X39, X40)
    
    # Covariate set with and without intercept
    X0_40 <- cbind(1, X_all_40)  # with intercept
    X_40 <- X_all_40  # without intercept
  }
  
  if(covars100){ # additional 60 covars using a variety of distributions including normal, uniform, Poisson, exponential, geometric, hypergeometric, and others. 
    # Additional code for 60 more covariates to make it 100
    X41 <- rgeom(n, prob=0.3)
    X42 <- rhyper(n, m=30, n=20, k=10)
    X43 <- rlogis(n, location=1, scale=2)
    X44 <- rweibull(n, shape=3, scale=2)
    X45 <- rt(n, df=10)
    X46 <- rwilcox(n, m=10, n=10)
    X47 <- rcauchy_truncated(n, location=1, scale=2)
    X48 <- rf(n, df1=10, df2=10)
    X49 <- rgamma(n, shape=2, scale=2)
    # Using custom function rfactorial
    X50 <- rfactorial(n, size=10)
    X51 <- rlnorm_truncated(n, meanlog=1, sdlog=2)
    X52 <- rnbinom(n, size=20, prob=0.4)
    X53 <- rsignrank(n, n=20)
    X54 <- rbeta(n, shape1=2, shape2=2)
    X55 <- rchisq(n, df=3)
    X56 <- rbinom(n, size=20, prob=0.4)
    X57 <- sample(1:10, n, replace = TRUE)  # Custom categorical variable
    X58 <- rexp_truncated(n, rate=0.2)
    X59 <- rpois(n, lambda=4)
    X60 <- rt(n, df=10)
    X61 <- runif(n, 0, 20)
    X62 <- rweibull(n, shape=3, scale=2)
    X63 <- rnorm(n, mean=10, sd=5)
    X64 <- rchisq(n, df=2)
    X65 <- rexp_truncated(n, rate=0.2)
    X66 <- rgeom(n, prob=0.4)
    X67 <- rhyper(n, m=30, n=30, k=15)
    X68 <- rlnorm_truncated(n, meanlog=1, sdlog=2)
    X69 <- rnbinom(n, size=20, prob=0.4)
    X70 <- rpois(n, lambda=5)
    X71 <- rt(n, df=15)
    X72 <- runif(n, 0, 20)
    X73 <- rweibull(n, shape=4, scale=2)
    X74 <- rchisq(n, df=3)
    X75 <- rexp_truncated(n, rate=0.3)
    X76 <- rgeom(n, prob=0.5)
    X77 <- rhyper(n, m=40, n=40, k=20)
    X78 <- rlnorm_truncated(n, meanlog=2, sdlog=2)
    X79 <- rnbinom(n, size=25, prob=0.5)
    X80 <- rpois(n, lambda=6)
    X81 <- rt(n, df=20)
    X82 <- runif(n, 0, 25)
    X83 <- rweibull(n, shape=5, scale=2)
    X84 <- rnorm(n, mean=15, sd=5)
    X85 <- rchisq(n, df=4)
    X86 <- rexp_truncated(n, rate=0.4)
    X87 <- rgeom(n, prob=0.6)
    X88 <- rhyper(n, m=50, n=50, k=25)
    X89 <- rlnorm_truncated(n, meanlog=2, sdlog=3)
    X90 <- rnbinom(n, size=30, prob=0.6)
    X91 <- rpois(n, lambda=7)
    X92 <- rt(n, df=25)
    X93 <- runif(n, 0, 30)
    X94 <- rweibull(n, shape=6, scale=3)
    X95 <- rnorm(n, mean=20, sd=6)
    X96 <- rchisq(n, df=5)
    X97 <- rexp_truncated(n, rate=0.5)
    X98 <- rgeom(n, prob=0.7)
    X99 <- rhyper(n, m=60, n=60, k=30)
    X100 <- rlnorm_truncated(n, meanlog=3, sdlog=3)
    
    # Combine 100 covariates
    X_all_100 <- cbind(X13, X4, X5, X6, X7, X8, X9, X10, X11, X12, X13_new, X14, X15, X16_new, X17, X18, X19, X20, X21, X22, X23, X24, X25, X26, X27, X28, X29, X30, X31, X32, X33, X34, X35, X36, X37, X38, X39, X40, X41, X42, X43, X44, X45, X46, X47, X48, X49, X50, X51, X52, X53, X54, X55, X56, X57, X58, X59, X60, X61, X62, X63, X64, X65, X66, X67, X68, X69, X70, X71, X72, X73, X74, X75, X76, X77, X78, X79, X80, X81, X82, X83, X84, X85, X86, X87, X88, X89, X90, X91, X92, X93, X94, X95, X96, X97, X98, X99, X100)
    
    # Covariate set with and without intercept
    X0_100 <- cbind(1, X_all_100)  # with intercept
    X_100 <- X_all_100  # without intercept
  }
  
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
  
  # Decide which covariate set to use
  if (covars40) {
    chosen_covars <- X0_40
    num_covars <- 41  # 40 covariates + intercept
  } else if (covars100) {
    chosen_covars <- X0_100
    num_covars <- 101  # 100 covariates + intercept
  } else {
    chosen_covars <- X06  # Default to the original 6 covariates and an intercept
    num_covars <- 7
  }
  
  # Adjust the matrix sizes to accommodate the number of covariates
  beta_mat <- matrix(0, nrow=num_covars, ncol=J)
  w_mat <- matrix(0, nrow=num_covars, ncol=J)  # Initialize to zeros or any other value
  xb_mat <- matrix(0, nrow=n, ncol=J)
  exp_xb_mat <- xb_mat
  e_mat <- xb_mat
  e_denom <- 0
  
  if(J==6){ 
    if (num_covars == 7) { # yang et al. values
      w_mat[,1] <- c(0, 0, 0, 0, 0, 0, 0)
      w_mat[,2] <- c(0, 1, 1, 2, 1, 1, 1)
      w_mat[,3] <- c(0, 1, 1, 1, 1, 1, -5)
      w_mat[,4] <- c(0, 1, 1, 1, 1, 1, 5)
      w_mat[,5] <- c(0, 1, 1, 1, -2, 1, 1)
      w_mat[,6] <- c(1.5, 6, 1, 2,-1,-1,-1)
    }else if (num_covars == 41) {
      w_mat[,1] <- rep(0, 41)
      w_mat[,2] <- c(0, rep(c(1, 1, 2, 1, 1, 1), 6), 1, 1, 1, 1)
      w_mat[,3] <- c(0, rep(c(1, -5), 20))
      w_mat[,4] <- c(0, rep(c(1, 5), 20))
      w_mat[,5] <- c(0, rep(c(1, -2), 20))
      w_mat[,6] <- c(1.5, rep(c(6, 1, 2, -1, -1, -1), 6), 6, 1, 2, -1)
    } else if (num_covars == 101) {
      w_mat[,1] <- rep(0, 101)
      w_mat[,2] <- c(0, rep(c(1, 1, 2, 1, 1, 1), 15), rep(1, 10))
      w_mat[,3] <- c(0, rep(c(1, -5), 50))
      w_mat[,4] <- c(0, rep(c(1, 5), 50))
      w_mat[,5] <- c(0, rep(c(1, -2), 50))
      w_mat[,6] <- c(1.5, rep(c(6, 1, 2, -1, -1, -1), 15), rep(6, 1), rep(1, 9))
    }
  } else if(J==3){
    w_mat[,1] <- c(0, 0, 0, 0, 0, 0, 0)
    w_mat[,2] <- c(0, 1, 1, 1, -1, 1, 1)
    w_mat[,3] <- c(0, 1, 1, 1, 1, 1, 1)
  }
  
  # Simulate observed treatment 
  D <- matrix(NA, n, J) # store obs. treatment
  colnames(D) <- paste0("D", rep(1:J))
  
  for(j in 1:J){
    beta_mat[,j] <- kappa[j]*w_mat[,j]
    xb_mat[,j] <- c(chosen_covars %*% beta_mat[,j])
    exp_xb_mat[,j] <- exp(xb_mat[,j])
    e_denom <- e_denom + exp_xb_mat[,j]
  }
  
  # calculate e_mat
  for(j in 1:J){
    e_mat[, j] <- exp_xb_mat[,j] / e_denom
  }
  
  e_mat[k, is.na(e_mat[k, ])] <- 1 # replace any NA with 1
  
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
      base_values <- c(-1.5, 1, 1, 1, 1, 1, 1, 
                                 -3, 2, 3, 1, 2, 2, 2, 
                                 3, 3, 1, 2, -1, -1, -4, 
                                 2.5, 4, 1, 2, -1, -1, -3, 
                                 2, 5, 1, 2, -1, -1, -2, 
                                 1.5, 6, 1, 2, -1, -1, -1)
      
      # Create an empty matrix of dimensions J x n_covars
      gamma_mat <- matrix(0, nrow = J, ncol = num_covars)
      
      # Populate the gamma_mat by repeating the pattern for each group
      for (j in 1:J) {
        start_idx <- (j - 1) * 7 + 1
        end_idx <- j * 7
        row_pattern <- base_values[start_idx:end_idx]
        
        # Repeat or truncate the pattern to fit n_covars
        repeated_values <- rep_len(row_pattern, num_covars)
        
        # Assign to the j-th row of gamma_mat
        gamma_mat[j, ] <- repeated_values
      }
    } else if(J==3){
      gamma_mat <- matrix(data=c(-1.5, 1, 1, 1, 1, 1, 1, 
                                 -3, 2, 3, 1, 2, 2, 2, 
                                 1.5, 3, 1, 2, -1, -1, -1),
                          nrow = 3, ncol = 7, byrow=TRUE)
    }

  } else if(gamma.setting == "li"){ # values in li and li
    if(J==6){
      base_values <- c(-1.5, 1, 1, 1, 1, 1, 1, 
                                 -4, 2, 3, 1, 2, 2, 2, 
                                 4, 3, 1, 2, -1, -1, -4, 
                                 1, 4, 1, 2, -1, -1, -3, 
                                 3.5, 5, 1, 2, -1, -1, -2, 
                                 3.5, 6, 1, 2, -1, -1, -1)
      # Create an empty matrix of dimensions J x n_covars
      gamma_mat <- matrix(0, nrow = J, ncol = num_covars)
      
      # Populate the gamma_mat by repeating the pattern for each group
      for (j in 1:J) {
        start_idx <- (j - 1) * 7 + 1
        end_idx <- j * 7
        row_pattern <- base_values[start_idx:end_idx]
        
        # Repeat or truncate the pattern to fit n_covars
        repeated_values <- rep_len(row_pattern, num_covars)
        
        # Assign to the j-th row of gamma_mat
        gamma_mat[j, ] <- repeated_values
      }
    } else if(J==3){
      gamma_mat <- matrix(data=c(-1.5, 1, 1, 1, 1, 1, 1, 
                                 -4, 2, 3, 1, 2, 2, 2, 
                                 3, 3, 1, 2, -1, -1, -1),
                          nrow = 3, ncol = 7, byrow=TRUE)
    }

  }else if(gamma.setting == "low"){ # 
    if(J==6){
      base_values <- c(-4, 1, -2, -1, 1, 1, 1, 
                                 -6, 1, -2, -1, 1, 1, 1, 
                                 -2, 1, -1, -1, -1, -1, -4, 
                                  1, 2, 1, 2, -1, -1, -3, 
                                 -2, 2, -1, 1, -2, -1, -3, 
                                 -3, 3, -1, 1, -2, -1, -2)
      # Create an empty matrix of dimensions J x n_covars
      gamma_mat <- matrix(0, nrow = J, ncol = num_covars)
      
      # Populate the gamma_mat by repeating the pattern for each group
      for (j in 1:J) {
        start_idx <- (j - 1) * 7 + 1
        end_idx <- j * 7
        row_pattern <- base_values[start_idx:end_idx]
        
        # Repeat or truncate the pattern to fit n_covars
        repeated_values <- rep_len(row_pattern, num_covars)
        
        # Assign to the j-th row of gamma_mat
        gamma_mat[j, ] <- repeated_values
      }
    } else if(J==3){
      gamma_mat <- matrix(data=c(-4, 1, -2, -1, 1, 1, 1, 
                                 -2, 1, -1, -1, -1, -1, -4, 
                                 -3, 3, -1, 1, -2, -1, -2),
                          nrow = 3, ncol = 7, byrow=TRUE)
    }
    
  }else if(gamma.setting == "zero"){
    gamma_mat <- matrix(0, nrow = J, ncol = num_covars, byrow=TRUE)
  } 

  EY <- list()
  
  if(outcome.type=="continuous"){
    u  <- rnorm(n)
    
    for(j in 1:J){
      EY[[j]] <- c(chosen_covars %*% gamma_mat[j,] + D[,j]) # incl. treatment dummies
    }
    
    EYs <- sapply(EY,unlist)
    Y.true <- EYs + u
    Y.obs <- rowSums(EYs*D) + u
  }else{
    Y <- list()
    for(j in 1:J){
      EY[[j]] <- exp(chosen_covars %*% gamma_mat[j,] + D[,j]) # incl. treatment dummies
      
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

  # We will now be agnostic to the DGP
  # and pretend as if we only have observed the analysis.data
  
  # Create the analysis datasets
  if(covars40){
    analysis.data <- data.frame(Y=Y.obs, Z=Z, X=X_40)
    colnames(analysis.data) <- c("Y", "Z", paste0("X",seq(1,40)))
  }else if(covars100){
    analysis.data <- data.frame(Y=Y.obs, Z=Z, X=X_100)
    colnames(analysis.data) <- c("Y", "Z", paste0("X",seq(1,100)))
  }else{
    analysis.data <- data.frame(Y=Y.obs, Z=Z, X=X16)
    colnames(analysis.data) <- c("Y", "Z", "X1", "X2", "X3", "X4", "X5", "X6")
  }

  if(scale.continuous){
    all_vars <- colnames(analysis.data)[-c(1:2)]
    
    # Specify the names of non-continuous variables
    if(covars100){
      non_continuous_vars <- c("X6", "X50", "X57", "X76", "X86", "X98", 
                               "X42", "X67", "X77", "X88", "X99", 
                               "X52", "X69", "X79", "X90", 
                               "X53", "X63", 
                               "X91", "X59", "X70", "X80")
    }else{
      non_continuous_vars <- c("X6")
    }

    
    # Identify the names of continuous variables
    continuous_vars <- setdiff(all_vars, non_continuous_vars)
    
    analysis.data[,continuous_vars] <- scale(analysis.data[,continuous_vars]) # scale continuous vars
  }
  
  # return generated data
  return(list("data"=analysis.data,
              "observed.treatment"=D,
              "trueATE"=true.ates,
              "TrueYs"=colMeans(Y.true),
              "C"=C))
}