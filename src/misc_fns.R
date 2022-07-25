######################
# Misc. functions    #
######################

# function to bound probabilities to be used when making predictions
boundProbs <- function(x,bounds=c(0.001,0.999)){
  x[x>max(bounds)] <- max(bounds)
  x[x<min(bounds)] <- min(bounds)
  return(x)
}

# proper characters
proper <- function(s) sub("(.)", ("\\U\\1"), tolower(s), pe=TRUE)

# definition: left-truncated normal distribution

rnorm_trunc <- function(n, mean, sd, minval = 0, maxval = 10000,
                        min.low = 0, max.low = 50, min.high = 5000, max.high = 10000)
{                                                                                      
  out <- rnorm(n = n, mean = mean, sd = sd)
  minval <- minval[1]; min1 <- min.low[1]; max1 <- max.low[1]
  maxval <- maxval[1]; min2 <- min.high[1]; max2 <- max.high[1]
  leq.zero <- length(out[out <= minval])
  geq.max <- length(out[out >= maxval])
  out[out <= minval] <- runif(n = leq.zero, min = min1, max = max1)
  out[out >= maxval] <- runif(n = geq.max, min = min2, max = max2)
  out
}

# definition: negative binomial
neg_binom <- function(n, mu){                                                                                      
  rnbinom(n=n, size=1, mu=mu)
}

# definition: multinomial
multinom <- function(n, prob){
  rmultinom(n=n, size=1, prob=prob)
}

# Summary figure for estimates
ForestPlot <- function(d, xlab, ylab){
  # Forest plot for summary figure
  p <- ggplot(d, aes(x=x, y=y, ymin=y.lo, ymax=y.hi,colour=forcats::fct_rev(Analysis))) + 
    geom_pointrange(size=1, position = position_dodge(width = -0.5)) + 
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