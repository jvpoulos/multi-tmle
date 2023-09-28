######################
# Misc. functions    #
######################

# Custom rfactorial function to generate random numbers from a factorial distribution
rfactorial <- function(n, size) {
  # n: number of random numbers to generate
  # size: the upper limit for the factorial distribution, from 1 to size
  if (size <= 0 || n <= 0) {
    stop("Both n and size should be positive integers.")
  }
  
  # Generate n random numbers from the set 1:size with equal probability
  return(sample(1:size, n, replace = TRUE))
}

# custom truncated log normal distribution
rlnorm_truncated <- function(n, meanlog = 0, sdlog = 1, upper_limit = 10) {
  x <- rlnorm(n, meanlog, sdlog)
  x[x > upper_limit] <- upper_limit
  return(x)
}

# custom truncated exponential distribution
rexp_truncated <- function(n, rate = 1, upper_limit = 10) {
  x <- rexp(n, rate)
  x[x > upper_limit] <- upper_limit
  return(x)
}

# truncated Cauchy Distribution
rcauchy_truncated <- function(n, location = 0, scale = 1, lower_limit = -10, upper_limit = 10) {
  x <- rcauchy(n, location, scale)
  x[x > upper_limit] <- upper_limit
  x[x < lower_limit] <- lower_limit
  return(x)
}

# function to bound probabilities to be used when making predictions
boundProbs <- function(x,bounds=c(0.001,0.999)){
  x[x>max(bounds)] <- max(bounds)
  x[x<min(bounds)] <- min(bounds)
  return(x)
}

# proper characters
proper <- function(s) sub("(.)", ("\\U\\1"), tolower(s), pe=TRUE)

# Summary figure for estimates
ForestPlot <- function(d, xlab, ylab){
  # Forest plot for summary figure
  p <- ggplot(d, aes(x=x, y=y, ymin=y.lo, ymax=y.hi,colour=forcats::fct_rev(Analysis))) + 
    geom_pointrange(size=1, position = position_dodge(width = -0.5), alpha=0.7) + 
    coord_flip() +
    geom_hline(aes(yintercept=0), lty=2) +
    ylab(xlab) +
    xlab(ylab) #switch because of the coord_flip() above
  return(p)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}