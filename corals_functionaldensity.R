rm(list = ls())
library(fda)
library(robCompositions)
library(RColorBrewer)

set.seed(12345)

# make polygon representing 95% pointwise confidence band by bootstrap
makepolygon95 <- function(y, t.fine){ # y is a 2 dimensional array, size of the array is number of values of t.fine 
                                    # and number of bootstrap iterations, t.fine is x values for polygon
  q025 <- apply(y, 1, "quantile", 0.025)
  q975 <- apply(y, 1, "quantile", 0.975)
  polygon(c(t.fine, rev(t.fine)), c(q025,rev(q975)), col = adjustcolor("gray", alpha.f = 0.5), 
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
#Note: we probably should not trust this, even though in our model it agrees well with the bootstrap version
make_asymp_polygon <- function(splinemodel, Z, i, t.fine, f){
  Omegafull <- vcov(splinemodel)
  Omega <- Omegafull[i, i]
  V <- Z %*% Omega %*% t(Z)
  se <- sqrt(diag(V))
  clow <- f - 1.96 * se
  chigh <- f + 1.96 * se
  polygon(c(t.fine, rev(t.fine)), c(clow, rev(chigh)), col = adjustcolor("grey", alpha.f = 0.5), 
          border = NA)
}

#plot a single coefficient function
#Arguments:
#t.fine: points at which to evaluate function
#xlab: x axis label
#ylab: y axis label
#bootstrap: use bootstrap for confidence band? Default TRUE
#betaboot: array of bootstrap fits (number of parameters x nt.fine x R, where nt.fine is number of points in t.fine and R is number of bootstrap replicates). Can be NULL if bootstrap not done
#whichparm: which parameter do we want? Ordered as in vcov
#fittedsplinemodel: object returned by fitZBmodel()
#ZB_base: matrix (points in t.fine x splines): basis for ZB-splines evaluated at points in t.fine
#coefindices: indices of parameter of interest in vcov: only used if not doing bootstrap, can be NULL otherwise (default)
#textlabel (default ""): label to add to top left corner of plot
#Value: plots the coefficient function, with a confidence band
plotcoefficientfunction <- function(t.fine, xlab, ylab, bootstrap = TRUE, betaboot, whichparm, fittedsplinemodel, ZB_base, coefindices = NULL, textlabel = ""){
  par(mfrow = c(1,1))
  plot(range(t.fine),range(betaboot[whichparm, , ]), type = "n", xlab = xlab, ylab = ylab, cex.lab = 1.5, cex.axis = 1.5 )
  if(bootstrap){
    makepolygon95(y = betaboot[whichparm, , ], t.fine = t.fine)
  } else {
    make_asymp_polygon(splinemodel = fittedsplinemodel$splinemodel, Z = ZB_base, i = coefindices, t.fine = t.fine, f = fittedsplinemodel$comp.spline.clr[, whichparm])
  }
  lines(t.fine, fittedsplinemodel$comp.spline.clr[, whichparm])
  abline(a = 0, b = 0, lty = "dashed")
  axls <- par("usr")
  text(axls[1] + 0.05 * (axls[2] - axls[1]), axls[4] - 0.05 * (axls[4] - axls[3]), textlabel, pos = 4, cex = 1.5)
}

# Numerical integration via trapezoidal formula (on a regular grid)
#copied from code supplied with Talska et al. 2018
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
##coef (matrix, sites x dimension of spline space) of ZB-spline coefficients of observations
#explanatory: data frame or matrix of explanatory variables (rows are sites, columns are variables)
#ZB_base: matrix (points in t.fine x splines): basis for ZB-splines evaluated at points in t.fine
#t.fine: points at which to evaluate function
#Value: list containing
#Fperm: array (nperm x number of values of t.fine): distribution of pointwise functional F statistics when observations permuted
#Fobs: vector (length number of values of t.fine): observed pointwise functional F statistics
#Fcrit vector (length number of values of t.fine): (1 - alpha) quantile of distribution of pointwise functional F statistics when observations permuted
#Fmaxperm: vector (length nperm) of max values of F(t) in each permutation
#Fmaxobs: scalar, max observed value of F(t)
#Fmaxcrit: (1 - alpha)-quantile of distribution of max functional F statistics when observations permuted
#FmaxP: randomization P-value based on Fmax (they call it permutation, but we don't look at all possible permutations)
functionalF <- function(smoothedobservations, y_pred.l, nperm = 1e3, Falpha = 0.05, coef, explanatory, ZB_base, t.fine){
  N <- dim(smoothedobservations)[2]
  nt <- dim(smoothedobservations)[1]
  Fobs <- fitZBmodel(coef = coef, explanatory = explanatory, ZB_base = ZB_base, nsites = N, t.fine = t.fine)$F #observed vector of pointwise functional F statistics
  Fperm <- array(dim = c(nperm, nt))
  for(i in 1:nperm){
    iperm <- sample(1:N, replace = FALSE)
    permutedmodel <- fitZBmodel(coef = coef[iperm, ], explanatory = explanatory, ZB_base = ZB_base, nsites = N, t.fine = t.fine) #permute the rows of coefficients, equivalent to permuting the observations
    Fperm[i, ] <- permutedmodel$F
  }
  Fcrit <- apply(rbind(Fperm, Fobs), 2, "quantile", 1 - Falpha)
  Fmaxobs <- max(Fobs) #observed max(F(t))
  Fmaxperm <- apply(Fperm, 1, "max") #max(F(t)) in each permutation
  Fmaxcrit <- quantile(c(Fmaxperm, Fmaxobs), 1 - Falpha)
  FmaxP <- (sum(Fmaxperm >= Fmaxobs) + 1) / (nperm + 1) #include the observed max(F(t))
  return(list(Fperm = Fperm, Fobs = Fobs, Fcrit = Fcrit, Fmaxperm = Fmaxperm, Fmaxobs = Fmaxobs, Fmaxcrit = Fmaxcrit, FmaxP = FmaxP))
}

#fit a ZB-spline model
#Arguments:
#coef (matrix, sites x dimension of spline space) of ZB-spline coefficients of observations
#explanatory: data frame or matrix of explanatory variables at each site
#ZB_base: matrix (points in t.fine x splines): basis for ZB-splines evaluated at points in t.fine
#nsites: number of sites
#t.fine: vector, values of log colony area at which to evaluate splines
#Value:
#list containing:
#y_pred.l: matrix (points in t.fine x sites) of clr predictions
#smoothedobservations: matrix (points in t.fine x sites) of clr smoothed observations
#F: vector of pointwise F statistics
#comp.spline.clr: matrix (points in t.fine x 3) of coefficient functions (intercept, PC1, PC2)
#B: matrix (# explanatory x dimension of spline space) of regression coefficients
#splinemodel: fitted lm object
fitZBmodel <- function(coef, explanatory, ZB_base, nsites, t.fine){
    explanatory <- as.matrix(explanatory)
    splinemodel <- lm(coef ~ explanatory)
    B <-  coef(splinemodel)  
    comp.spline.clr <- ZB_base %*% t(B)
    X <- cbind(1, explanatory)
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
        # fname: name of pdf file to which we will write (if not NA, default)
        # oldpar: default par settings to restore (helps make figure predictable if drawing in RStudio plot panel, where size not known)
pc1predictions <-  function(axisscores, fittedsplinemodel, nt.fine, t.fine, t_step, fname = NA, oldpar){
  npc1 <- 10
  pc1grid <- seq(from = min(axisscores$PC1), to = max(axisscores$PC1), length.out = npc1)
  Xgrid <- as.matrix(cbind(rep(1, npc1), pc1grid, rep(0, npc1))) #third column 0 : mean of PC2
  pc1gridclr <- fittedsplinemodel$comp.spline.clr %*% t(Xgrid) #clr predictions on a grid of equally spaced PC1 scores from min to max, with PC2 = 0 (mean value)
  pc1gridpred <- array(dim = c(npc1, nt.fine))
  for(i in 1:npc1){
    pc1gridpred[i,] <- clr2density(t.fine,t_step,pc1gridclr[,i])
    
  }
  pc1colors <- brewer.pal(10, "RdBu")
  if(!is.na(fname)){
    pdf(file = fname, width = 7, height = 7)
  }
  par(oldpar, no.readonly = TRUE) #reset to original to get predictable behaviour
  par(mar = c(5, 6, 2, 2))
  matplot(t.fine, t(pc1gridpred), type = "l", lty = "solid", xlab = expression(paste("log(coral area/"*cm^2*")")), ylab = "Probability density", col = pc1colors, cex.lab = 1.5, cex.axis = 1.5)
  par(fig = c(grconvertX(c(7, 8), from = "user", to = "ndc"), grconvertY(c(0.15, 0.3), from = "user", to = "ndc")), mar = c(1, 1, 2, 1), new = TRUE) #coordinates in which we'll plot the color bar image
  image(y = 1:10, z = t(1:10), col = pc1colors, axes = FALSE, xlab = NA, ylab = NA, main = "PC1") #use image to plot the color bar
  axis(4, at = 1:10, las = 2, labels = round(pc1grid, 2), col = NA, col.ticks = NA, cex.axis = 1) #label the colors, rotate the labels and make the axis itself invisible
  if(!is.na(fname)){
    dev.off()
  }
  par(oldpar, no.readonly = TRUE) #reset to original to get predictable behaviour
}

# residual plot coloured by PC1 scores (blue (more positive)-red (more negative))
#  input: t.fine = grid of log coral areas
# residua = residuals from fitted model (array, rows are log coral area values, columns are sites) 
# nsites = number of sites
# sites = vector of names of sites
# axisscores = PCA axis scores
residualplot <- function(t.fine, residua, nsites, sites, axisscores){
  colorfunc <- colorRampPalette(brewer.pal(8, "RdBu"))
  matplot(t.fine, residua, type = "l", lty = "solid", xlab = expression(paste("log(coral area/"*cm^2*")")), ylab = "clr residuals", cex.lab = 1.5, cex.axis = 1.5, col = colorfunc(8)[findInterval(axisscores$PC1,seq(from = min(axisscores$PC1),to = max(axisscores$PC1),length.out = 8))])
  for(i in 1:nsites){
    text(t.fine[1],
      residua[1,i],sites[i],pos = 4, cex = 0.9)
  }
  for(i in 1:nsites){
   text(t.fine[length(t.fine)],
    residua[length(t.fine),i],sites[i],pos = 2, cex = 0.9)
  }
}

#compute pointwise and global functional R-squared. Based on Talska et al. 2018 (p. 79 and code in supporting info)
#arguments:
#fittedsplinemodel: object returned by fitZBmodel()
#t.fine: grid of values of log coral area at which to evaluate functions (vector)
#t_step: step size in t.fine
#nsites: number of sites
#textlabel: extra label (default "") giving more information about what was plotted
#Value: list containing 
#R.t pointwise R2
#R2global
#also plots pointwise and global R2
computeR2 <- function(fittedsplinemodel, t.fine, t_step, nsites, textlabel = ""){
  #compute pointwise and global R^2
  mean.l <- apply(fittedsplinemodel$smoothedobservations, 1, "mean") #note this equals mean of predictions, provided we have an intercept in the model (Talska code used mean of predictions)
  SSE <- rowSums((fittedsplinemodel$smoothedobservations - fittedsplinemodel$y_pred.l)^2)
  SST <- rowSums((fittedsplinemodel$smoothedobservations - mean.l)^2)
  SSF <- rowSums((fittedsplinemodel$y_pred.l - mean.l)^2)
  R.t <- SSF / SST #pointwise R^2
  plot(t.fine, R.t, xlab = expression(paste("log(coral area/"*cm^2*")")), ylab = expression(paste("Pointwise", ~italic(R)^2)), ylim = c(0, 1), type = "l", cex.lab = 1.5, cex.axis = 1.5)
  
  SST.norm <- 0
  SSF.norm <- 0
  for(i in 1:nsites){
    SST.norm <- SST.norm + trapzc(t_step, (fittedsplinemodel$smoothedobservations[, i] - mean.l)^2)
    SSF.norm <- SSF.norm + trapzc(t_step, (fittedsplinemodel$y_pred.l[, i] - mean.l)^2)
  }
  R2global <- SSF.norm / SST.norm
  print(paste("global R^2:", R2global, sep = " "))
  legend("topright", bty = "n", cex = 1.5, legend = bquote(paste(italic(R)[global]^2==.(round(R2global, 2)))))
  axls <- par("usr")
  text(axls[1] + 0.05 * (axls[2] - axls[1]), axls[4] - 0.05 * (axls[4] - axls[3]), textlabel, pos = 4)
  
  return(list(R.t = R.t, R2global = R2global))
}

#smooth a set of log areas using compositionalSpline: by default, choose number of histogram bins using Sturges' rule, choose smoothing parameter alpha by generalized cross-validation
#Arguments:
#site: name of site
#oneyeardf: data frame containing variables logArea and Site
#nalpha: number of values of smoothing parameter alpha (between 0.01 and 1) to try on regular grid (default 30)
#knots: sequence of knots in compositional splines
#order: order of splines
#nc (default NULL): number of bins (if NULL, use Sturges' rule)
#Value: list containing
#nc: number of bins (selected by Sturges' rule if argument nc NULL)
#breaks: histogram breaks
#dens.raw: density in each bin (after zero count imputation)
#t.raw: bin midpoints
#clr.raw: clr-transformed density in each bin
#gcv: generalized cross-validation score for each value of alpha
#alphas: vector of smoothing parameters considered
#alpha: smoothing parameter value that minimizes gcv score
#coef: coefficients of ZB-splines, estimated using selected value of alpha
smoothhistogram <- function(site, oneyeardf, nalpha = 30, knots, order, nc = NULL){
  if(is.null(nc)){
    nc <- nclass.Sturges(oneyeardf$logArea[oneyeardf$Site == site]) #Apply Sturges' rule to data for site i
  }
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
#fittedsplinemodel: object returned by fitZBmodel()
#t.fine: grid of log areas on which clr densities evaluated
#sites: list of site names
#shists: list of objects returned by smoothhistogram() for each site
#oneyeardf: data frame including variable Site
#Value: plot with a panel for each site, clr density on y axis, log area on x axis, points are raw clr densities, solid lines are smoothed clr densities, dashed lines are predicted smoothed clr densities
plotfit <- function(fittedsplinemodel, t.fine, sites, shists, oneyeardf){
  par(mfrow = c(4,5)) # 4,4 for reduced sites
  par(mar = c(4, 5, 2, 2))
  yl <- range(range(fittedsplinemodel$smoothedobservations), range(fittedsplinemodel$y_pred.l))
  for (i in 1:nsites){
    plot(t.fine, fittedsplinemodel$smoothedobservations[, i], type = "l", main = paste(letters[i], ": ", sites[i], sep = ""), ylim = yl, xlab = expression(paste("log(coral area/"*cm^2*")")), ylab = "clr density", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5) #smoothed raw data
    lines(t.fine, fittedsplinemodel$y_pred.l[, i], lty = "dashed") #predicted
    points(shists[[i]]$t.raw, shists[[i]]$clr.raw, pch = 16, col = adjustcolor("black", 0.4)) #raw clr densities in bins
    if(i == 1){
      legend("bottomright", bty = "n", pch = c(16, NA, NA), col = c(adjustcolor("black", 0.4), "black", "black"), lty = c(NA, "solid", "dashed"), legend = c("raw", "smoothed", "predicted"), cex = 1.2)
    }
    axlims <- par("usr")
    text(axlims[4] - 0.2 * (axlims[2] - axlims[1]), axlims[3] + 0.2 * (axlims[4] - axlims[3]), bquote(alpha == .(round(shists[[i]]$alpha, 2))), cex = 1.5)
    ni <- dim(oneyeardf[oneyeardf$Site == sites[i], ])[1]
    text(axlims[4] - 0.2 * (axlims[2] - axlims[1]), axlims[3] + 0.1 * (axlims[4] - axlims[3]), bquote(italic(n)[i] == .(ni)), cex = 1.5)
    
  }
}

#refit model excluding sites with few colonies
#Arguments:
#nthreshold: min number of colonies at a site for inclusion
#oneyeardf: data frame containing variables logArea and Site
#g: number of knots is g + 2
#k: order of splines: #2 is quadratic, 3 is cubic, etc
#nalpha: number of values of smoothing parameter alpha (between 0.01 and 1) to try on regular grid (default 30)
#sites: vector of names of sites
#axisscores: data frame, variables PC1, PC2 (second and third columns are explanatory variables)
#knots: sequence of knots in compositional splines
#order: order of splines
#ZB_base: matrix (points in t.fine x splines): basis for ZB-splines evaluated at points in t.fine
#nt.fine: number of points at which to evaluate function
#t.fine: points at which to evaluate function
#bootstrap (logical, default TRUE): bootstrap the fitted model to get estimate of uncertainty
#R (default 1000: number of bootstrap replicates)
#Value: just plots the coefficient functions
refitwithoutsmallsites <- function(nthreshold, oneyeardf, g, k, nalpha, sites, axisscores, knots, order, ZB_base, nt.fine, t.fine, bootstrap = TRUE, R = 1000){
  ninclude <- table(oneyeardf$Site) > nthreshold
  siteorder <- match(sites, names(ninclude)) #order of sites is not alphabetical
  ninclude <- ninclude[siteorder]
  nsitestrim <- sum(ninclude)
  sitestrim <- sites[ninclude]
  axisscorestrim <- axisscores[ninclude, ]
  order <- k + 1
  coeftrim <- matrix(nrow = nsitestrim, ncol = g + k) #Machalova et al 2021, Theorem 1: dimension of the ZB-spline vector space is g + k
  gcvtrim <- array(dim = c(nsitestrim, nalpha))
  alphastrim <- seq(from = 0.1, to = 1, length.out = nalpha)
  selectedalphastrim <- rep(NA, nsitestrim) #smoothing parameter for each site
  shiststrim <- vector(mode = "list", length = nsitestrim) #list of lists: smoothed histogram data for each site
  for (i in 1:nsitestrim){
    shiststrim[[i]] <- smoothhistogram(site = sitestrim[i], oneyeardf = oneyeardf, nalpha = nalpha, knots = knots, order = order)
    gcvtrim[i, ] <- shiststrim[[i]]$gcv
    coeftrim[i, ] <- shiststrim[[i]]$coef
    selectedalphastrim[i] <- shiststrim[[i]]$alpha
  }
  
  # ZB-spline basis evaluated on the grid "t.fine"
  fittedsplinemodeltrim <- fitZBmodel(coef = coeftrim, explanatory = axisscorestrim, ZB_base = ZB_base, nsites = nsitestrim, t.fine = t.fine)
  betaboot <- bootstrapsplinemodel(fittedsplinemodel = fittedsplinemodeltrim, R = R, nt.fine = nt.fine, nsites = nsitestrim, explanatory = axisscorestrim)
  coefindices <- seq(from = 1, to = dim(vcov(fittedsplinemodeltrim$splinemodel))[1], by = 3) #every third row/column in covariance matrix of parameters is intercept, because we have intercept and two explanatory variables
  plotcoefficientfunction(t.fine = t.fine, xlab = expression(paste("log(coral area/"*cm^2*")")), ylab = "clr(density)", bootstrap = bootstrap, betaboot = betaboot, whichparm = 1, fittedsplinemodel = fittedsplinemodeltrim, ZB_base = ZB_base, coefindices = coefindices, textlabel = bquote(beta[0]*", sites with more than"~.(nthreshold)~"colonies"))
  plotcoefficientfunction(t.fine = t.fine, xlab = expression(paste("log(coral area/"*cm^2*")")), ylab = "clr(density)", bootstrap = bootstrap, betaboot = betaboot, whichparm = 2, fittedsplinemodel = fittedsplinemodeltrim, ZB_base = ZB_base, coefindices = coefindices + 1, textlabel = bquote(beta[1]*", sites with more than"~.(nthreshold)~"colonies"))
  plotcoefficientfunction(t.fine = t.fine, xlab = expression(paste("log(coral area/"*cm^2*")")), ylab = "clr(density)", bootstrap = bootstrap, betaboot = betaboot, whichparm = 3, fittedsplinemodel = fittedsplinemodeltrim, ZB_base = ZB_base, coefindices = coefindices + 2, textlabel = bquote(beta[2]*", sites with more than"~.(nthreshold)~"colonies"))
  print(paste("Excluding sites with fewer than", nthreshold, "colonies:"))
  Rsquaredtrim <- computeR2(fittedsplinemodel = fittedsplinemodeltrim, t.fine = t.fine, t_step = diff(t.fine)[1], nsites = nsitestrim, textlabel = bquote("Only sites with more than"~.(nthreshold)~"colonies")) #compute and plot pointwise and global R-squared
}

#categorize each site as pre- or post-bleaching
#Arguments:
#axisscores: PC1 and PC2 of environmental variables at each site
#oneyeardf: data frame including Site and Year variables
#Value: vector, with sites ordered as in rows of axisscores, with TRUE if site has Year either 2016b or 2018, FALSE otherwise
makebleachvariable <- function(axisscores, oneyeardf){
  oneyeardf$postbleach <- oneyeardf$Year == "2016b" | oneyeardf$Year == "2018"
  bleachtable <- table(oneyeardf$Site, oneyeardf$postbleach)
  sitepostbleach <- bleachtable[, 2] > 0 #all observations at site were in same year
  siteorder <- match(sites, rownames(bleachtable)) #order of sites is not alphabetical
  sitepostbleach <- sitepostbleach[siteorder]
  return(sitepostbleach)
}

#create array of bootstrapped coefficients
#Arguments:
#fittedsplinemodel: object returned by fitZBmodel()
#R: number of bootstrap iterations
#nt.fine: number of points at which to evaluate functions
#nsites: number of sites
#explanatory: data frame or matrix of explanatory variables (rows are sites, columns are variables)
#Value:
#array (number of explanatory variables + 1 x nt.fine x R) of bootstrapped coefficients
bootstrapsplinemodel <- function(fittedsplinemodel, R, nt.fine, nsites, explanatory){
  resids  <- fittedsplinemodel$smoothedobservations - fittedsplinemodel$y_pred.l
  nvars <- dim(explanatory)[2]
  betaboot <- array(dim = c(nvars + 1, nt.fine, R))
  for (i in 1:R){# generate new dataset, fit model to new dataset, keeping coefs
    j <- sample(1:nsites, replace = TRUE) #resample set of residuals 
    yboot <- t(fittedsplinemodel$y_pred.l + resids[, j]) # generate new dataset, fit model, keeping coefs
    betaboot[, , i] <- coef(lm(yboot ~ as.matrix(explanatory))) 
  }
  return(betaboot)
}

#permutation F-test
#Arguments:
#Falpha: size of test, e.g. 0.05
#nperm: number of permutations
#fittedsplinemodel: object returned by fit_ZBmodel()
#coef (matrix, sites x dimension of spline space) of ZB-spline coefficients of observations
#ZB_base: matrix (points in t.fine x splines): basis for ZB-splines evaluated at points in t.fine
#explanatory: data frame or matrix of explanatory variables (rows are sites, columns are variables)
#t.fine: points at which to evaluate function
#textlabel: extra label for plot (default "")
#Value:
#Plots the pointwise functional F
#Ftest, from functionalF(), a list containing
#Fperm: array (nperm x number of values of t.fine): distribution of pointwise functional F statistics when observations permuted
#Fobs: vector (length number of values of t.fine): observed pointwise functional F statistics
#Fcrit vector (length number of values of t.fine): (1 - alpha) quantile of distribution of pointwise functional F statistics when observations permuted
#Fmaxperm: vector (length nperm) of max values of F(t) in each permutation
#Fmaxobs: scalar, max observed value of F(t)
#Fmaxcrit: (1 - alpha)-quantile of distribution of max functional F statistics when observations permuted
#FmaxP: randomization P-value based on Fmax (they call it permutation, but we don't look at all possible permutations)
doFtest <- function(Falpha, nperm, fittedsplinemodel, coef, ZB_base, explanatory, t.fine, textlabel = ""){
  Ftest <- functionalF(smoothedobservations = fittedsplinemodel$smoothedobservations, y_pred.l = fittedsplinemodel$y_pred.l, nperm = nperm, Falpha = Falpha, coef = coef, ZB_base = ZB_base, explanatory = explanatory, t.fine = t.fine)
  par(mfrow = c(1,1))
  plot(t.fine, Ftest$Fobs, type = "l", xlab = expression(paste("log(coral area/"*cm^2*")")), ylab = expression(paste("Pointwise", ~italic(F))), ylim = c(0, max(c(Ftest$Fobs, Ftest$Fmaxcrit))), cex.lab = 1.5, cex.axis = 1.5)
  lines(t.fine, Ftest$Fcrit, lty = "dotted")
  abline(h = Ftest$Fmaxcrit, lty = "dashed")
  legend(x = 6, y = 0.9, bty = "n", lty = c("solid", "dotted", "dashed"), legend = c("observed", as.expression(bquote(paste("pointwise ", .(Falpha), " critical value"))), as.expression(bquote(paste("maximum ", .(Falpha), " critical value")))))
  print(paste("Randomization Fmax test: observed Fmax = ", Ftest$Fmaxobs, "P =", Ftest$FmaxP, "from", nperm, "randomizations"))
  legend(x = 0, y = 0.9, bty = "n", cex = 1.5, legend = bquote(paste(italic(P)==.(round(Ftest$FmaxP, 2)))))
  axls <- par("usr")
  text(axls[1] + 0.05 * (axls[2] - axls[1]), axls[4] - 0.05 * (axls[4] - axls[3]), textlabel, pos = 4)
  return(Ftest)
}

oldpar <- par(no.readonly = TRUE) #default par settings (restore them to get predictable behaviour)

#import data
oneyeardf <- read.csv("oneyeardf.csv", row.names=1)
oneyeardf <- oneyeardf[is.na(oneyeardf$ROI.LabelCode), ] # to remove corals that were out of frame 
axisscores <- read.csv("axisscores.csv", row.names=1)
axisscores <- axisscores[order(axisscores$PC1), ]

sites <- unique(row.names(axisscores))
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

fittedsplinemodel <- fitZBmodel(coef = coef, explanatory = axisscores, ZB_base = ZB_base, nsites = nsites, t.fine = t.fine)
plotfit(fittedsplinemodel = fittedsplinemodel, t.fine = t.fine, sites = sites, shists = shists, oneyeardf = oneyeardf)

bootstrap <- TRUE #CI by bootstrap? Otherwise by asymptotic method on which we probably can't rely, although results similar
R <- 1e3
betaboot <- bootstrapsplinemodel(fittedsplinemodel = fittedsplinemodel, R = R, nt.fine = nt.fine, nsites = nsites, explanatory = axisscores)
residua  <- fittedsplinemodel$smoothedobservations - fittedsplinemodel$y_pred.l 

par(mfrow = c(1,1))
coefindices <- seq(from = 1, to = dim(vcov(fittedsplinemodel$splinemodel))[1], by = 3) #every third row/column in covariance matrix of parameters is intercept, because we have intercept and two explanatory variables (but won't actually use these: don't have the asymptotic theory for them)
plotcoefficientfunction(t.fine = t.fine, xlab = expression(paste("log(coral area/"*cm^2*")")), ylab = "clr(density)", bootstrap = bootstrap, betaboot = betaboot, whichparm = 1, fittedsplinemodel = fittedsplinemodel, ZB_base = ZB_base, coefindices = coefindices)
plotcoefficientfunction(t.fine = t.fine, xlab = expression(paste("log(coral area/"*cm^2*")")), ylab = "clr(density)", bootstrap = bootstrap, betaboot = betaboot, whichparm = 2, fittedsplinemodel = fittedsplinemodel, ZB_base = ZB_base, coefindices = coefindices + 1)
plotcoefficientfunction(t.fine = t.fine, xlab = expression(paste("log(coral area/"*cm^2*")")), ylab = "clr(density)", bootstrap = bootstrap, betaboot = betaboot, whichparm = 3, fittedsplinemodel = fittedsplinemodel, ZB_base = ZB_base, coefindices = coefindices + 2)
Rsquared <- computeR2(fittedsplinemodel = fittedsplinemodel, t.fine = t.fine, t_step = t_step, nsites = nsites) #compute and plot pointwise and global R-squared
  
#permutation F-test
Falpha <- 0.05
nperm <- 1e4 - 1
Ftest <- doFtest(Falpha = Falpha, nperm = nperm, fittedsplinemodel = fittedsplinemodel, coef = coef, ZB_base = ZB_base, explanatory = axisscores, t.fine = t.fine)

# figure showing predicted size distributions with increasing PC1
pc1predictions(axisscores = axisscores, fittedsplinemodel = fittedsplinemodel, nt.fine = nt.fine, t.fine = t.fine, t_step = t_step, fname = NA, oldpar = oldpar)

#residuals
residualplot(t.fine = t.fine, residua = residua, nsites = nsites, sites = sites, axisscores = axisscores)

#refit without sites that have only a small number of colonies
nthreshold <- 400
refitwithoutsmallsites(nthreshold = nthreshold, oneyeardf = oneyeardf, g = g, k = k, nalpha = nalpha, sites = sites, axisscores = axisscores, knots = knots, order = order, ZB_base = ZB_base, nt.fine = nt.fine, t.fine = t.fine, bootstrap = bootstrap, R = R)

#consider including a pre- vs post-2016 bleaching variable
bleach <- makebleachvariable(axisscores = axisscores, oneyeardf = oneyeardf)
bleachexplanatory <- cbind(axisscores, bleach)
fsmbleach <- fitZBmodel(coef = coef, explanatory = bleachexplanatory, ZB_base = ZB_base, nsites = nsites, t.fine = t.fine)
plotfit(fittedsplinemodel = fsmbleach, t.fine = t.fine, sites = sites, shists = shists, oneyeardf = oneyeardf)
betabootbleach <- bootstrapsplinemodel(fittedsplinemodel = fsmbleach, R = R, nt.fine = nt.fine, nsites = nsites, explanatory = bleachexplanatory)
residbleach  <- fsmbleach$smoothedobservations - fsmbleach$y_pred.l 

par(mfrow = c(1,1))
plotcoefficientfunction(t.fine = t.fine, xlab = expression(paste("log(coral area/"*cm^2*")")), ylab = "clr(density)", bootstrap = bootstrap, betaboot = betabootbleach, whichparm = 1, fittedsplinemodel = fsmbleach, ZB_base = ZB_base, coefindices = NULL, textlabel = bquote(beta[0]~"with bleaching variable included")) #intercept
plotcoefficientfunction(t.fine = t.fine, xlab = expression(paste("log(coral area/"*cm^2*")")), ylab = "clr(density)", bootstrap = bootstrap, betaboot = betabootbleach, whichparm = 2, fittedsplinemodel = fsmbleach, ZB_base = ZB_base, coefindices = NULL, textlabel = bquote(beta[1]~"with bleaching variable included")) #PC1
plotcoefficientfunction(t.fine = t.fine, xlab = expression(paste("log(coral area/"*cm^2*")")), ylab = "clr(density)", bootstrap = bootstrap, betaboot = betabootbleach, whichparm = 3, fittedsplinemodel = fsmbleach, ZB_base = ZB_base, coefindices = NULL, textlabel = bquote(beta[2]~"with bleaching variable included")) #PC2
plotcoefficientfunction(t.fine = t.fine, xlab = expression(paste("log(coral area/"*cm^2*")")), ylab = "clr(density)", bootstrap = bootstrap, betaboot = betabootbleach, whichparm = 4, fittedsplinemodel = fsmbleach, ZB_base = ZB_base, coefindices = NULL, textlabel = bquote(beta[3]*", bleaching variable")) #bleach variable
print("With bleaching variable included:")
Rsquaredbleach <- computeR2(fittedsplinemodel = fsmbleach, t.fine = t.fine, t_step = t_step, nsites = nsites, textlabel = "with bleaching variable included") #compute and plot pointwise and global R-squared
Ftestbleach <- doFtest(Falpha = Falpha, nperm = nperm, fittedsplinemodel = fsmbleach, coef = coef, ZB_base = ZB_base, explanatory = bleachexplanatory, t.fine = t.fine, textlabel = "with bleaching variable included")

#check effects of number of bins: fit with min and max from range selected by Sturges' rule over all sites
nbins <- unlist(lapply(shists, "[[", "nc")) #extracts number of bins at each site

coefmin <- matrix(nrow = nsites, ncol = g + k) #Machalova et al 2021, Theorem 1: dimension of the ZB-spline vector space is g + k
gcvmin <- array(dim = c(nsites, nalpha))
alphasmin <- seq(from = 0.1, to = 1, length.out = nalpha)
selectedalphasmin <- rep(NA, nsites) #smoothing parameter for each site
shistsmin <- vector(mode = "list", length = nsites) #list of lists: smoothed histogram data for each site
coefmax <- matrix(nrow = nsites, ncol = g + k) #Machalova et al 2021, Theorem 1: dimension of the ZB-spline vector space is g + k
gcvmax <- array(dim = c(nsites, nalpha))
alphasmax <- seq(from = 0.1, to = 1, length.out = nalpha)
selectedalphasmax <- rep(NA, nsites) #smoothing parameter for each site
shistsmax <- vector(mode = "list", length = nsites) #list of lists: smoothed histogram data for each site
for (i in 1:nsites){
  shistsmin[[i]] <- smoothhistogram(site = sites[i], oneyeardf = oneyeardf, nalpha = nalpha, knots = knots, order = order, nc = min(nbins))
  gcvmin[i, ] <- shistsmin[[i]]$gcv
  coefmin[i, ] <- shistsmin[[i]]$coef
  selectedalphasmin[i] <- shistsmin[[i]]$alpha
  shistsmax[[i]] <- smoothhistogram(site = sites[i], oneyeardf = oneyeardf, nalpha = nalpha, knots = knots, order = order, nc = max(nbins))
  gcvmax[i, ] <- shistsmax[[i]]$gcv
  coefmax[i, ] <- shistsmax[[i]]$coef
  selectedalphasmax[i] <- shistsmax[[i]]$alpha
}
fittedsplinemodelmin <- fitZBmodel(coef = coefmin, explanatory = axisscores, ZB_base = ZB_base, nsites = nsites, t.fine = t.fine)
plotfit(fittedsplinemodel = fittedsplinemodelmin, t.fine = t.fine, sites = sites, shists = shistsmin, oneyeardf = oneyeardf)
fittedsplinemodelmax <- fitZBmodel(coef = coefmax, explanatory = axisscores, ZB_base = ZB_base, nsites = nsites, t.fine = t.fine)
plotfit(fittedsplinemodel = fittedsplinemodelmax, t.fine = t.fine, sites = sites, shists = shistsmax, oneyeardf = oneyeardf)
betaboot <- bootstrapsplinemodel(fittedsplinemodel = fittedsplinemodel, R = R, nt.fine = nt.fine, nsites = nsites, explanatory = axisscores)
residua  <- fittedsplinemodel$smoothedobservations - fittedsplinemodel$y_pred.l 

par(mfrow = c(1,1))
betabootmin <- bootstrapsplinemodel(fittedsplinemodel = fittedsplinemodelmin, R = R, nt.fine = nt.fine, nsites = nsites, explanatory = axisscores)
betabootmax <- bootstrapsplinemodel(fittedsplinemodel = fittedsplinemodelmax, R = R, nt.fine = nt.fine, nsites = nsites, explanatory = axisscores)
plotcoefficientfunction(t.fine = t.fine, xlab = expression(paste("log(coral area/"*cm^2*")")), ylab = "clr(density)", bootstrap = bootstrap, betaboot = betabootmin, whichparm = 1, fittedsplinemodel = fittedsplinemodelmin, ZB_base = ZB_base, coefindices = NULL, textlabel = bquote(beta[0]~"with"~.(min(nbins))~"bins at all sites"))
plotcoefficientfunction(t.fine = t.fine, xlab = expression(paste("log(coral area/"*cm^2*")")), ylab = "clr(density)", bootstrap = bootstrap, betaboot = betabootmin, whichparm = 2, fittedsplinemodel = fittedsplinemodelmin, ZB_base = ZB_base, coefindices = NULL, textlabel = bquote(beta[1]~"with"~.(min(nbins))~"bins at all sites"))
plotcoefficientfunction(t.fine = t.fine, xlab = expression(paste("log(coral area/"*cm^2*")")), ylab = "clr(density)", bootstrap = bootstrap, betaboot = betabootmin, whichparm = 3, fittedsplinemodel = fittedsplinemodelmin, ZB_base = ZB_base, coefindices = NULL, textlabel = bquote(beta[2]~"with"~.(min(nbins))~"bins at all sites"))
plotcoefficientfunction(t.fine = t.fine, xlab = expression(paste("log(coral area/"*cm^2*")")), ylab = "clr(density)", bootstrap = bootstrap, betaboot = betabootmax, whichparm = 1, fittedsplinemodel = fittedsplinemodelmax, ZB_base = ZB_base, coefindices = NULL, textlabel = bquote(beta[0]~"with"~.(max(nbins))~"bins at all sites"))
plotcoefficientfunction(t.fine = t.fine, xlab = expression(paste("log(coral area/"*cm^2*")")), ylab = "clr(density)", bootstrap = bootstrap, betaboot = betabootmax, whichparm = 2, fittedsplinemodel = fittedsplinemodelmax, ZB_base = ZB_base, coefindices = NULL, textlabel = bquote(beta[1]~"with"~.(max(nbins))~"bins at all sites"))
plotcoefficientfunction(t.fine = t.fine, xlab = expression(paste("log(coral area/"*cm^2*")")), ylab = "clr(density)", bootstrap = bootstrap, betaboot = betabootmax, whichparm = 3, fittedsplinemodel = fittedsplinemodelmax, ZB_base = ZB_base, coefindices = NULL, textlabel = bquote(beta[2]~"with"~.(max(nbins))~"bins at all sites"))
