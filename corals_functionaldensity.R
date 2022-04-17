rm(list = ls())
library(fda)
library(robCompositions)
library(RColorBrewer)

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

# Numerical integration via trapezoidal formula (on a regular grid)
#copied from code suppied with Talska et al. 2018
#Arguments: 
#step: step size for grid
#y: function to be integrated, evaluated at the grid points
#Value: the value of the required integral
trapzc <- function(step, y) 
{
  return(step * (0.5 * y[1] + sum(y[2:(length(y) - 1)]) + 0.5 * y[length(y)]))
}

#functional F statistic as in Ramsay et al. 2009, Functional Data Analysis with R and Matlab, p. 168, eq. 10.20
#Arguments:
#smoothedobservations (array, number of values of t.fine x number of observations): smoothed clr observations
#y_pred.l (array, number of values of t.fine x number of observations): predicted clr observations
#Value: functional F statistic vector (length is number of values of t.fine)
Ft <- function(smoothedobservations, y_pred.l){
  Fnum <- apply(y_pred.l, 1, "var")
  n <- dim(smoothedobservations)[2]
  Fdenom <- 1 / n * rowSums((smoothedobservations - y_pred.l)^2)
  return(Fnum / Fdenom)
}

#functional F permutation test as in Ramsay et al. 2009, Functional Data Analysis with R and Matlab, p. 168
#Arguments:
#smoothedobservations (array, number of values of t.fine x number of observations): smoothed clr observations
#y_pred.l (array, number of values of t.fine x number of observations): predicted clr observations
#nperm: number of permutations (default 1e3)
#Falpha: size of test (default 0.05)
#Value: list containing
#Fperm: array (nperm x number of values of t.fine): distribution of pointwise functional F statistics when observations permuted
#Fobs: vector (length number of values of t.fine): observed pointwise functional F statistics
#Fcrit vector (length number of values of t.fine): (1 - alpha) quantile of distribution of pointwise functional F statistics when observations permuted
#Fmaxperm: vector (length nperm) of max values of F(t) in each permutation
#Fmaxobs: scalar, max observed value of F(t)
#Fmaxcrit: (1 - alpha)-quantile of distribution of max functional F statistics when observations permuted
functionalF <- function(smoothedobservations, y_pred.l, nperm = 1e3, Falpha = 0.05, coef, axisscores, ZB_base, t.fine){
  N <- dim(smoothedobservations)[2]
  nt <- dim(smoothedobservations)[1]
  Fobs <- fitZBmodel(coef = coef, axisscores = axisscores, ZB_base = ZB_base, nsites = N, t.fine = t.fine)$F #observed vector of pointwise functional F statistics
  Fperm <- array(dim = c(nperm, nt))
  for(i in 1:nperm){
    iperm <- sample(1:N, replace = FALSE)
    permutedmodel <- fitZBmodel(coef = coef[iperm, ], axisscores = axisscores, ZB_base = ZB_base, nsites = N, t.fine = t.fine) #permute the rows of coefficients, equivalent to permuting the observations
    Fperm[i, ] <- permutedmodel$F
  }
  Fcrit <- apply(Fperm, 2, "quantile", 1 - Falpha)
  Fmaxobs <- max(Fobs) #observed max(F(t))
  Fmaxperm <- apply(Fperm, 1, "max") #max(F(t)) in each permutation
  Fmaxcrit <- quantile(Fmaxperm, 1 - Falpha)
  return(list(Fperm = Fperm, Fobs = Fobs, Fcrit = Fcrit, Fmaxperm = Fmaxperm, Fmaxobs = Fmaxobs, Fmaxcrit = Fmaxcrit))
}

#fit a ZB-spline model
#Arguments:
#coef (matrix, sites x dimension of spline space) of ZB-spline coefficients of observations
#axisscores: data frame, variables PC1, PC2 (second and third columns are explanatory variables)
#ZB_base: matrix (points in t.fine x splines): basis for ZB-splines evaluated at points in t.fine
#nsites: number of sites
#t.fine: vector, values of log colony area at which to evaluate splines
#Value:
#list containing:
#y_pred.l: matrix (points in t.fine x sites) of clr predictions
#smoothedobservations: matrix (points in t.fine x sites) of clr smoothed observations
#F: vector of pointwise F statistics
#comp.spline.clr: matrix (points in t.fine x 3) of coefficient functions (intercept, PC1, PC2)
#B: matrix (3 x dimension of spline space) of regression coefficients
#splinemodel: fitted lm object
fitZBmodel <- function(coef, axisscores, ZB_base, nsites, t.fine){
  splinemodel <- lm(coef ~ axisscores$PC1 + axisscores$PC2)
  B <-  coef(splinemodel)  
  comp.spline.clr <- ZB_base %*% t(B)
  X <- as.matrix(cbind(rep(1, nsites), axisscores[, 1:2]))
  y_pred.l <- comp.spline.clr %*% t(X) #clr predictions
  smoothedobservations <- ZB_base %*% t(coef) # smoothed clr observations
  F <- Ft(smoothedobservations = smoothedobservations, y_pred.l = y_pred.l)
  return(list(y_pred.l = y_pred.l, smoothedobservations = smoothedobservations, F = F, comp.spline.clr = comp.spline.clr, B = B, splinemodel = splinemodel))
}

# Inverse of clr transformation from Talska et al 2018
# Input: z = grid of point defining the abscissa 
#        z_step = step of the grid of the abscissa
#        clr = grid evaluation of the clr transformed density
# Output: grid evaluation of the density
clr2density <- function(z, z_step, clr) 
{
  if(is.fd(clr))
    return(exp(eval.fd(z,clr))/trapzc(z_step,exp(eval.fd(z,clr))))
  if(!is.fd(clr))
    return(exp(clr)/trapzc(z_step,exp(clr)))
}

# PC1 predictions plot
# Input: axisscores = PCA axis scores
        # fittedsplinemodel = object returned by fitZBmodel
        # nt.fine = number of values on the x-axis (log coral area)
        # t.fine = grid of log coral areas
        # t_step = interval between the values in t.fine
pc1predictions <-  function(axisscores, fittedsplinemodel, nt.fine, t.fine, t_step){
  npc1 <- 10
  pc1grid <- seq(from = min(axisscores$PC1), to = max(axisscores$PC1), length.out = npc1)
  Xgrid <- as.matrix(cbind(rep(1, npc1), pc1grid, rep(0, npc1))) #third column 0 : mean of PC2
  pc1gridclr <- fittedsplinemodel$comp.spline.clr %*% t(Xgrid) #clr predictions on a grid of equally spaced PC1 scores from min to max, with PC2 = 0 (mean value)
  pc1gridpred <- array(dim = c(npc1, nt.fine))
  for(i in 1:npc1){
    pc1gridpred[i,] <- clr2density(t.fine,t_step,pc1gridclr[,i])
  }
  pc1colors <- brewer.pal(10, "RdBu")
  matplot(t.fine, t(pc1gridpred), type = "l", lty = "solid", xlab = "log coral area", ylab = "probability density", col = pc1colors)
}

# residual plot coloured by PC1 scores (blue (more positive)-red (more negative))
#  input: t.fine = grid of log coral areas
# residua = residuals from fitted model (array, rows are log coral area values, columns are sites) 
# nsites = number of sites
# sites = vector of names of sites
# axisscores = PCA axis scores
residualplot <- function(t.fine, residua, nsites, sites, axisscores){
  colorfunc <- colorRampPalette(brewer.pal(8, "RdBu"))
  matplot(t.fine, residua, type = "l", lty = "solid", xlab = "log coral area", ylab = "clr residuals", col = colorfunc(8)[findInterval(axisscores$PC1,seq(from = min(axisscores$PC1),to = max(axisscores$PC1),length.out = 8))])
  for(i in 1:nsites){
    text(t.fine[1],
      residua[1,i],sites[i],pos = 4, cex = 0.5)
  }
  for(i in 1:nsites){
   text(t.fine[length(t.fine)],
    residua[length(t.fine),i],sites[i],pos = 2, cex = 0.5)
  }
}

#compute pointwise and global functional R-squared. Based on Talska et al. 2018 (p. 79 and code in supporting info)
#arguments:
#fittedsplinemodel: object returned by fitZBmodel()
#t.fine: grid of values of log coral area at which to evaluate functions (vector)
#t_step: step size in t.fine
#Value: list containing 
#R.t pointwise R2
#R2global
#also plots pointwise and global R2
computeR2 <- function(fittedsplinemodel, t.fine, t_step){
  #compute pointwise and global R^2
  mean.l <- apply(fittedsplinemodel$smoothedobservations, 1, "mean") #note this equals mean of predictions, provided we have an intercept in the model (Talska code used mean of predictions)
  SSE <- rowSums((fittedsplinemodel$smoothedobservations - fittedsplinemodel$y_pred.l)^2)
  SST <- rowSums((fittedsplinemodel$smoothedobservations - mean.l)^2)
  SSF <- rowSums((fittedsplinemodel$y_pred.l - mean.l)^2)
  R.t <- SSF / SST #pointwise R^2
  plot(t.fine, R.t, xlab = "log coral area", ylab = expression(paste("pointwise", ~italic(R)^2)), ylim = c(0, 1), type = "l")
  
  SST.norm <- 0
  SSF.norm <- 0
  for(i in 1:nsites){
    SST.norm <- SST.norm + trapzc(t_step, (fittedsplinemodel$smoothedobservations[, i] - mean.l)^2)
    SSF.norm <- SSF.norm + trapzc(t_step, (fittedsplinemodel$y_pred.l[, i] - mean.l)^2)
  }
  R2global <- SSF.norm / SST.norm
  print(paste("global R^2:", R2global, sep = " "))
  legend("topright", bty = "n", legend = bquote(paste(italic(R)[global]^2==.(round(R2global, 2)))))
  
  return(list(R.t = R.t, R2global = R2global))
}

#smooth a set of log areas using compositionalSpline: choose number of histogram bins using Sturges' rule, choose smoothing parameter alpha by generalized cross-validation
#Arguments:
#site: name of site
#oneyeardf: data frame containing variables logArea and Site
#nalpha: number of values of smoothing parameter alpha (between 0.01 and 1) to try on regular grid (default 30)
#knots: number of knots in compositional splines
#order: order of splines
#Value: list containing
#nc: number of bins selected by Sturges' rule
#breaks: histogram breaks selected by Sturges' rule
#dens.raw: density in each bin (after zero count imputation)
#t.raw: bin midpoints
#clr.raw: clr-transformed density in each bin
#gcv: generalized cross-validation score for each value of alpha
#alphas: vector of smoothing parameters considered
#alpha: smoothing parameter value that minimizes gcv score
#coef: coefficients of ZB-splines, estimated using selected value of alpha
smoothhistogram <- function(site, oneyeardf, nalpha = 30, knots, order){
  nc <- nclass.Sturges(oneyeardf$logArea[oneyeardf$Site == site]) #Apply Sturges' rule to data for site i
  breaks <- seq(from = min(oneyeardf$logArea), to = max(oneyeardf$logArea), length.out = nc + 1)
  classwidth <- diff(breaks)[1] # assuming all classes have the same width, which is the case with nclass.Sturges()
  sitedens <- hist(oneyeardf$logArea[oneyeardf$Site == site], breaks = breaks, plot = FALSE) # get histogram density for just the subset of logArea values from each site
  sitedens$density[sitedens$density == 0] <- 2 / 3 * 1 / (sum(sitedens$counts)) # zero count imputation from Machalova et al 2021 p.1053
  dens.raw <- sitedens$density / sum(sitedens$density)
  t.raw <- matrix(sitedens$mids, nrow = 1) # mid points (t is what people called independent variables in FDA)
  clr.raw <- cenLR(matrix(dens.raw, nrow = 1))$x.clr # clr transformation
  alphas <- seq(from = 0.01, to = 1, length.out = nalpha)
  weights <- rep(1, nc)
  gcv <- array(dim = c(1, nalpha))
  for(j in 1:nalpha){#try out different values of smoothing parameter
    cspline <- compositionalSpline(t = t.raw, clrf = as.numeric(clr.raw), knots = knots, w = weights, order = order, der = 1, alpha = alphas[j], spline.plot = FALSE, basis.plot = FALSE)
    gcv[j] <- cspline$GCV #generalized cross-validation score
  }
  alpha <- alphas[which.min(gcv)] #choose the alpha from our grid that minimizes GCV score
  cspline <- compositionalSpline(t = t.raw, clrf = as.numeric(clr.raw), knots = knots, w = weights, order = order, der = 1, alpha = alpha, spline.plot = FALSE, basis.plot = FALSE) #refit with selected alpha
  return(list(nc = nc, breaks = breaks, dens.raw = dens.raw, t.raw = t.raw, clr.raw = clr.raw, gcv = gcv, alphas = alphas, alpha = alpha, coef = cspline$ZB_coef))  
}

#plot raw and smoothed histograms and predicted densities on clr scale
#Arguments:
#fittesplinemodel: object returned by fitZBmodel()
#t.fine: grid of log areas on which clr densities evaluated
#sites: list of site names
#shists: list of objects returned by smoothhistogram() for each site
#Value: plot with a panel for each site, clr density on y axis, log area on x axis, points are raw clr densities, solid lines are smoothed clr densities, dashed lines are predicted smoothed clr densities
plotfit <- function(fittedsplinemodel, t.fine, sites, shists){
  par(mfrow=c(4,5))
  par(mar = c(4, 5, 2, 2))
  yl <- range(range(fittedsplinemodel$smoothedobservations), range(fittedsplinemodel$y_pred.l))
  for (i in 1:nsites){
    plot(t.fine, fittedsplinemodel$smoothedobservations[, i], type = "l", main = paste(letters[i], ": ", sites[i], sep = ""), ylim = yl, xlab = "log coral area", ylab = "clr density") #smoothed raw data
    lines(t.fine, fittedsplinemodel$y_pred.l[, i], lty = "dashed") #predicted
    points(shists[[i]]$t.raw, shists[[i]]$clr.raw, pch = 16, col = adjustcolor("black", 0.4)) #raw clr densities in bins
    if(i == 1){
      legend("bottomleft", bty = "n", pch = c(16, NA, NA), col = c(adjustcolor("black", 0.4), "black", "black"), lty = c(NA, "solid", "dashed"), legend = c("raw", "smoothed", "predicted"))
    }
    axlims <- par("usr")
    text(axlims[2] - 0.2 * (axlims[2] - axlims[1]), axlims[3] + 0.1 * (axlims[4] - axlims[3]), bquote(alpha == .(round(shists[[i]]$alpha, 2))))
  }
}

#import data
oneyeardf <- read.csv("oneyeardb.csv", 
                      row.names=1)
oneyeardf <- oneyeardf[is.na(oneyeardf$ROI.LabelCode), ] # to remove corals that were out of frame 
axisscores <- read.csv("axisscores.csv", 
                       row.names=1)
axisscores <- axisscores[order(rownames(axisscores)), ]

# pcascores <- read.csv("pcascores.csv") # calling PCA scores from all environmental variables
# pcascores <- pcascores[order(pcascores$Site), ]


sites <- sort(unique(oneyeardf$Site))
nsites <- length(sites)

#compositional spline
g <- 2
nknots <- g + 2
knots <- seq(from = min(oneyeardf$logArea), to = max(oneyeardf$logArea),
             length.out = nknots) #g + 2 equally-spaced values
k <- 3 #2 is quadratic, 3 is cubic, etc
order <- k + 1

par(mfrow = c(1, 1))
# make matrix for coef
coef <- matrix(nrow = nsites, ncol = g + k) #Machalova et al 2021, Theorem 1: dimension of the ZB-spline vector space is g + k

nalpha <- 30
gcv <- array(dim = c(nsites, nalpha))
alphas <- seq(from = 0.1, to = 1, length.out = nalpha)
selectedalphas <- rep(NA, nsites) #smoothing parameter for each site
shists <- vector(mode = "list", length = nsites) #list of lists: smoothed histogram data for each site
for (i in 1:nsites){
  shists[[i]] <- smoothhistogram(site = sites[i], oneyeardf = oneyeardf, nalpha = nalpha, knots = knots, order = order)
  gcv[i, ] <- shists[[i]]$gcv
  coef[i, ] <- shists[[i]]$coef
  selectedalphas[i] <- shists[[i]]$alpha
}
plot(range(alphas), range(gcv), type = "n", xlab = expression(alpha), ylab = "GCV")
for(i in 1:nsites){
  lines(alphas, gcv[i, ])
}

# ZB-spline basis evaluated on the grid "t.fine"
nt.fine <- 1000
t.fine <- seq(from = min(knots), to = max(knots), length.out = nt.fine)
t_step <- diff(t.fine)[1]
ZB_base <- ZBsplineBasis(t = t.fine, knots = knots, order = order)$ZBsplineBasis

fittedsplinemodel <- fitZBmodel(coef = coef, axisscores = axisscores, ZB_base = ZB_base, nsites = nsites, t.fine = t.fine)
plotfit(fittedsplinemodel = fittedsplinemodel, t.fine = t.fine, sites = sites, shists = shists)

plot(t.fine, fittedsplinemodel$comp.spline.clr[, 1], type = "l")
lines(t.fine, fittedsplinemodel$comp.spline.clr[, 2], col = "red")
lines(t.fine, fittedsplinemodel$comp.spline.clr[, 3], col = "blue")
abline(a = 0, b = 0, lty = "dashed")

plot(t.fine, fittedsplinemodel$y_pred.l[, 1], type = "l")
for (i in 2:nsites){
  lines(t.fine, fittedsplinemodel$y_pred.l[, i])
}


bootstrap <- FALSE

  # Bootstrap
  # functional residuals
  residua  <- fittedsplinemodel$smoothedobservations - fittedsplinemodel$y_pred.l 
  
  # compute bootstrap response Y_boot, R bootstrap repetitions
  R <- 1000  
  
  betaboot <- array(dim = c(3, nt.fine, R))
  
  # generate new dataset, fit model to new dataset, keeping coefs
  for (i in 1:R){
    j <- sample(1:nsites, replace = TRUE) #resample set of residuals 
    yboot <- t(fittedsplinemodel$y_pred.l + residua[, j]) # generate new dataset, fit model, keeping coefs
    betaboot[, , i] <- coef(lm(yboot ~ axisscores$PC1 + axisscores$PC2)) 
  }


par(mfrow = c(1,1))
plot(range(t.fine),range(betaboot[1, , ]), type = "n", xlab = "log coral area", ylab = "clr of intercept" )
if(bootstrap){
  makepolygon95(y = betaboot[1, , ], t.fine = t.fine)
}
lines(t.fine, fittedsplinemodel$comp.spline.clr[, 1])
abline(a = 0, b = 0, lty = "dashed")
coefindices <- seq(from = 1, to = dim(vcov(fittedsplinemodel$splinemodel))[1], by = 3) #every third row/column in covariance matrix of parameters is intercept, because we have intercept and two explanatory variables
make_asymp_polygon(splinemodel = fittedsplinemodel$splinemodel, Z = ZB_base, i = coefindices, t.fine = t.fine, f = fittedsplinemodel$comp.spline.clr[, 1])

plot(range(t.fine),range(betaboot[2, , ]), type = "n", xlab = "log coral area", ylab = "clr of first axis scores" )
if(bootstrap){
  makepolygon95(y = betaboot[2, , ], t.fine = t.fine)
}
lines(t.fine, fittedsplinemodel$comp.spline.clr[, 2])
abline(a = 0, b = 0, lty = "dashed")
make_asymp_polygon(splinemodel = fittedsplinemodel$splinemodel, Z = ZB_base, i = coefindices + 1, t.fine = t.fine, f = fittedsplinemodel$comp.spline.clr[, 2])

plot(range(t.fine),range(betaboot[3, , ]), type = "n", xlab = "log coral area", ylab = "clr of second axis scores" )
if(bootstrap){
  makepolygon95(y = betaboot[3, , ], t.fine = t.fine)
}
lines(t.fine, fittedsplinemodel$comp.spline.clr[, 3])
abline(a = 0, b = 0, lty = "dashed")
make_asymp_polygon(splinemodel = fittedsplinemodel$splinemodel, Z = ZB_base, i = coefindices + 2, t.fine = t.fine, f = fittedsplinemodel$comp.spline.clr[, 3])

Rsquared <- computeR2(fittedsplinemodel = fittedsplinemodel, t.fine = t.fine, t_step = t_step) #compute and plot pointwise and global R-squared
  
#permutation F-test
Falpha <- 0.05
nperm <- 1e2 #WANT A BIGGER NUMBER E.G. 1E4 FOR FINAL VERSION
Ftest <- functionalF(smoothedobservations = fittedsplinemodel$smoothedobservations, y_pred.l = fittedsplinemodel$y_pred.l, nperm = nperm, Falpha = Falpha, coef = coef, ZB_base = ZB_base, axisscores = axisscores, t.fine = t.fine)
par(mfrow = c(1,1))
plot(t.fine, Ftest$Fobs, type = "l", xlab = "log coral area", ylab = expression(paste("pointwise", ~italic(F))), ylim = c(0, max(c(Ftest$Fobs, Ftest$Fmaxcrit))))
lines(t.fine, Ftest$Fcrit, lty = "dotted")
abline(h = Ftest$Fmaxcrit, lty = "dashed")
legend("topright", bty = "n", lty = c("solid", "dotted", "dashed"), legend = c("observed", as.expression(bquote(paste("pointwise ", .(Falpha), " critical value"))), as.expression(bquote(paste("maximum ", .(Falpha), " critical value")))))

# figure showing predicted size distributions with increasing PC1
pc1predictions(axisscores = axisscores, fittedsplinemodel = fittedsplinemodel, nt.fine = nt.fine, t.fine = t.fine, t_step = t_step)

#residuals
residualplot(t.fine = t.fine, residua = residua, nsites = nsites, sites = sites, axisscores = axisscores)
