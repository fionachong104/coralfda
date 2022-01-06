rm(list = ls())
library(fda)
library(robCompositions)

# make polygon
makepolygon95 <- function(y, t.fine){ # y is a 2 dimensional array, size of the array is number of values of t.fine 
                                    # and number of bootstrap iterations, t.fine is x values for polygon
  q025 <- apply(y, 1, "quantile", 0.025)
  q975 <- apply(y, 1, "quantile", 0.975)
  polygon(c(t.fine, rev(t.fine)), c(q025,rev(q975)), col = adjustcolor("green", alpha.f = 0.1), 
          border = NA)
}

#construct approximate 95% confidence band using asymptotic Gaussian approximation
#Arguments:
#splinemodel: an lm object used to estimate coefficients of ZB-splines
#Z: the ZB-spline basis evaluated at points in t.fine
#i: indices in vcov(splinemodel) of the variable we want (check the row names of vcov(splinemodel to get these))
#t.fine: points at which to evaluate function
#f: estimated values of function
#Value: draws a polygon representing approximate 95% confidence band around f at the points in t.fine
make_asymp_polygon <- function(splinemodel, Z, i, t.fine, f){
  Omegafull <- vcov(splinemodel)
  Omega <- Omegafull[i, i]
  V <- Z %*% Omega %*% t(Z)
  se <- sqrt(diag(V))
  clow <- f - 1.96 * se
  chigh <- f + 1.96 * se
  polygon(c(t.fine, rev(t.fine)), c(clow, rev(chigh)), col = adjustcolor("red", alpha.f = 0.1), 
          border = NA)
}

#import data
oneyeardf <- read.csv("oneyeardb.csv", 
                      row.names=1)
oneyeardf <- oneyeardf[is.na(oneyeardf$ROI.LabelCode), ] # to remove corals that were out of frame 
axisscores <- read.csv("axisscores.csv", 
                       row.names=1)
axisscores <- axisscores[order(axisscores$Site), ]



sites <- sort(unique(oneyeardf$Site))
nsites <- length(sites)

# check how many histogram bins are used, as per Sturges rule

for (i in sites){
  nc <- nclass.Sturges(oneyeardf$logArea[oneyeardf$Site==i])
  print(c(i, nc))
}

nclasses <- 9 # this might need changing later as data changes

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

nalpha <- 30
alphas <- seq(from = 0.01, to = 1, length.out = nalpha)
gcv <- array(dim = c(nsites, nalpha))
for (i in 1:nsites){
  for(j in 1:nalpha){#try out different values of smoothing parameter
    cspline <- compositionalSpline(t = t.raw, clrf = clr.raw[i,], 
                                   knots = knots, w = weights, order = order, der = 1, alpha = alphas[j],                                    spline.plot = FALSE, basis.plot = FALSE)
    gcv[i, j] <- cspline$GCV #generalized cross-validation score
  }
}
plot(range(alphas), range(gcv), type = "n", xlab = expression(alpha), ylab = "GCV")
for(i in 1:nsites){
  lines(alphas, gcv[i, ])
}

for(i in 1:nsites){
  cspline <- compositionalSpline(t = t.raw, clrf = clr.raw[i,], 
              knots = knots, w = weights, order = order, der = 1, alpha = 0.9, 
              spline.plot = TRUE, basis.plot = TRUE)
  coef[i,] <- cspline$ZB_coef
}

# ZB-spline basis evaluated on the grid "t.fine"
nt.fine <- 1000
t.fine <- seq(from = min(knots), to = max(knots),length.out = nt.fine)
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

# Bootstrap
# functional residuals
residua  <- smoothedobservations - y_pred.l 

# compute bootstrap response Y_boot, R bootstrap repetitions
R <- 1000  

betaboot <- array(dim = c(3, nt.fine, R))

# generate new dataset, fit model to new dataset, keeping coefs
for (i in 1:R){
  j <- sample(1:nsites, replace = TRUE) #resample set of residuals 
  yboot <- t(y_pred.l + residua[, j]) # generate new dataset, fit model, keeping coefs
  betaboot[, , i] <- coef(lm(yboot ~ axisscores$areax + axisscores$areay)) 
}
par(mfrow = c(1,1))
plot(range(t.fine),range(betaboot[1, , ]), type = "n", xlab = "log coral areas", ylab = "clr of intercept" )
makepolygon95(y = betaboot[1, , ], t.fine = t.fine)
lines(t.fine, comp.spline.clr[, 1])
abline(a = 0, b = 0, lty = "dashed")
make_asymp_polygon(splinemodel = splinemodel, Z = ZB_base, i = c(1, 4, 7, 10, 13), t.fine = t.fine, f = comp.spline.clr[, 1])

plot(range(t.fine),range(betaboot[2, , ]), type = "n", xlab = "log coral areas", ylab = "clr of first axis scores" )
makepolygon95(y = betaboot[2, , ], t.fine = t.fine)
lines(t.fine, comp.spline.clr[, 2])
abline(a = 0, b = 0, lty = "dashed")
make_asymp_polygon(splinemodel = splinemodel, Z = ZB_base, i = c(2, 5, 8, 11, 14), t.fine = t.fine, f = comp.spline.clr[, 2])

# 
# matlines(t.fine, betaboot[2,,], col = "grey", type = "l", lty = "solid")

plot(range(t.fine),range(betaboot[3, , ]), type = "n", xlab = "log coral areas", ylab = "clr of second axis scores" )
makepolygon95(y = betaboot[3, , ], t.fine = t.fine)
lines(t.fine, comp.spline.clr[, 3])
abline(a = 0, b = 0, lty = "dashed")
make_asymp_polygon(splinemodel = splinemodel, Z = ZB_base, i = c(3, 6, 9, 12, 15), t.fine = t.fine, f = comp.spline.clr[, 3])


# matplot(t.fine, betaboot[1,,], col = "grey", type = "l", lty = "solid", 
#         xlab = "log coral areas", ylab = "clr of intercept")
# lines(t.fine, comp.spline.clr[,1])
# abline(a = 0, b = 0, lty = "dashed")
# 
# matplot(t.fine, betaboot[2,,], col = "grey", type = "l", lty = "solid", 
#         xlab = "log coral areas", ylab = "clr of first axis scores")
# lines(t.fine, comp.spline.clr[,2], col = "red")
# abline(a = 0, b = 0, lty = "dashed")
# 
# matplot(t.fine, betaboot[3,,], col = "grey", type = "l", lty = "solid", 
#         xlab = "log coral areas", ylab = "clr of second axis scores")
# lines(t.fine, comp.spline.clr[,3], col = "blue")
# abline(a = 0, b = 0, lty = "dashed")

