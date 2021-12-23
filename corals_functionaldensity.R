rm(list = ls())
library(fda)
library(robCompositions)

#import data
oneyeardf <- read.csv("C:/Users/624225/Box/PhD/coralfda/oneyeardb.csv", 
                      row.names=1)
axisscores <- read.csv("C:/Users/624225/Box/PhD/coralfda/axisscores.csv", 
                       row.names=1)
axisscores <- axisscores[order(axisscores$Site), ]


nclasses <- 10
sites <- sort(unique(oneyeardf$Site))
nsites <- length(sites)

breaks <- seq(from = min(oneyeardf$logArea), to = max(oneyeardf$logArea), 
              length.out = nclasses + 1)
#breaks <- breaks[3:9]
classwidth <- diff(breaks)[1] # assuming all classes have the same width
#create empty matrix for things to go in
dens.raw <- matrix(nrow = nsites, ncol = nclasses)

# get histogram density for just the subset of logArea values from each site
for(i in 1:nsites){
  sitedens <- hist(oneyeardf$logArea[oneyeardf$Site==sites[i]],
                   breaks = breaks, plot = FALSE)
  # for the rows of the matrix
  sitedens$density[sitedens$density == 0] <- 2/3 * 1/(sum(sitedens$counts)) 
                      # zero count imputation from Machalova et al 2021 p.1053
  dens.raw[i,] <- sitedens$density / sum(sitedens$density)
}

# mid points (t is what people called independent variables in FDA)
t.raw <- sitedens$mids # assuming same classes for all sites

# clr transformation
clr.raw <- cenLR(dens.raw)$x.clr

#compositional spline
nknots <- 5
knots <- seq(from = min(oneyeardf$logArea), to = max(oneyeardf$logArea),
             length.out = nknots) # five test values equally spaced

weights <- rep(1,nclasses)

order <- 3

par(mfrow = c(1,1))
# make matrix for coef
coef <- matrix(nrow = nsites, ncol = nknots)
for (i in 1:nsites){
  cspline <- compositionalSpline(t = t.raw, clrf = clr.raw[i,], 
              knots = knots, w = weights, order = order, der = 1, alpha = 0.9, 
              spline.plot = TRUE, basis.plot = TRUE)
  coef[i,] <- cspline$ZB_coef
}

# ZB-spline basis evaluated on the grid "t.fine"
t.fine <- seq(from = min(knots), to = max(knots),length.out = 1000)
ZB_base <- ZBsplineBasis(t = t.fine, knots = knots, order = order)$ZBsplineBasis

# Compute estimation of B matrix by using command lm
splinemodel <- lm(coef ~ axisscores$areax + axisscores$areay)
B <-  coef(splinemodel)
#beta.fd <- fd(t(B), ZB_base)

comp.spline.clr <- ZB_base %*% t(B)

plot(t.fine, comp.spline.clr[,1], type = "l")
lines(t.fine, comp.spline.clr[,2], col = "red")
lines(t.fine, comp.spline.clr[,3], col = "blue")
abline(a = 0, b = 0, lty = "dashed")

# interpretation of betas 
y_pred.l <- matrix(nrow = length(t.fine), ncol = nsites)
for (i in 1:nsites){
  y_pred.l[,i] <- comp.spline.clr[,1] + comp.spline.clr[,2] * axisscores$areax[i]
                 + comp.spline.clr[,3] * axisscores$areay[i]
}

# y_pred = NULL
# for (i in 1:length(weight.l)){
#   y_pred = cbind(y_pred, clr2density(t.fine,t_step,y_pred.l[,i]))
# }

plot(t.fine, y_pred.l[,1], type = "l")
for (i in 2:nsites){
  lines(t.fine, y_pred.l[,i])
}

# smoothed clr observations
smoothedobservations <- ZB_base %*% t(coef)

par(mfrow=c(3,6))
for (i in 1:nsites){
  plot(t.fine, smoothedobservations[,i], type = "l", main = sites[i])
  lines(t.fine, y_pred.l[,i], col="blue")
}