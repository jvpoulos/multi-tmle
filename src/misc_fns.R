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