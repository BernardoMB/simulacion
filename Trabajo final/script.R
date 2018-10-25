library(ggplot2)
library(latex2exp)

# PDF Compound Poisson Distribution
dpoiscomp <- function(x, lambda_po=4,lambda_exp=2) {
  # Dont take into account terms of the sum that are smaller than the following tolerance
  tolerance <- 0.0000000001
  prob <- 0
  # The rest of terms
  n <- 1
  repeat {
    next_term <- exp(-lambda_po+n*log(lambda_po)-lfactorial(n)) * dgamma(x, shape=n, rate=lambda_exp)
    if (next_term < tolerance) {
      break
    }
    prob <- prob + next_term
    n <- n + 1
  }
  prob
}

# Get values to be ploted
getValues.pdf.pois.comp <- function(x, lambda_po=4, lambda_exp=2) sapply(x, FUN=dpoiscomp, lambda_po=lambda_po,lambda_exp=lambda_exp)

# PDF Gamma Trasladada
getValues.pdf.gamma.tras <- function(x, desp=1, alpha=4, beta=2) {
  dgamma(x-desp, shape = alpha, rate=beta)
}

# Create a plot
p <- ggplot()
p

# Stablish background and leyend position
theme <- theme(
  # Background
  #panel.grid.major = element_blank(),
  #panel.grid.minor = element_blank(),
  #panel.background = element_blank(),

  # Legend position
  legend.position = c(.915, .955),

  # Legend boxes
  legend.key = element_rect(colour = 'black', fill = 'white', size = .5, linetype=0),
  legend.key.width = unit(3.5, "lines"),

  # Center plot title
  plot.title = element_text(hjust = 0.5)
)
p + theme

# Plot axes to make it look nicer
xAxis <- geom_hline(yintercept = 0, size = .1)
yAxis <-geom_vline(xintercept = 0, size = .1)
p + theme + xAxis + yAxis

# Plot PDFs
# Labels
labels <- labs(x = TeX('$s$'), y = TeX('$f_S(s)$'))
# Title
title <- ggtitle("Densidad de la severidad total \n vs. \n aproximaci?n")
# Empty plot
p + theme + xAxis + yAxis + labels + title
# Number of different function to be ploted
types <- c("1","2")
# Create plot
p1 <- ggplot(data.frame(type=types), aes(colour=type, linetype=type, size=type))
# Create layers
pdf1 <- geom_path(
  data=data.frame(type="1"),
  stat="function",
  fun = getValues.pdf.pois.comp,
  args=list(lambda_po=0.1, lambda_exp=1/10000),
  aes(colour=type, linetype=type, size=type)
)
pdf2 <- stat_function(
  data=data.frame(type="2"),
  fun = getValues.pdf.gamma.tras,
  args = list(desp=-1000/3, alpha=4/45, beta=1/15000),
  aes(colour=type, linetype=type, size=type)
)
# Colors and line types
color <- FALSE
if (color) {
  colors <- c("blue","red")
  sizes <- c(1,1)
  linetypes <-c(1,1)
} else {
  colors <- c("black", "black")
  sizes <- c(1,1)
  linetypes <-c(1,5)
}
# Legends
an <- labs(TeX('Severidad total'), TeX('Aproximaci?n'))
# Adjust scales on axes
xlim <- xlim(0,80000)
ylim <- ylim(0,.00001)
# Plot
p1 + theme + xAxis + yAxis + labels + title + xlim + ylim + pdf1 + pdf2 +
  scale_colour_manual(
    #name="Densidades",
    name=NULL,
    values=colors,
    labels = an,
    breaks = types
  ) +
  scale_linetype_manual(
    #name="Densidades",
    name=NULL,
    values = linetypes,
    labels = an,
    breaks = types
  ) +
  scale_size_manual(
    #name="Densidades",
    name=NULL,
    values = sizes,
    labels = an,
    breaks = types
  )








