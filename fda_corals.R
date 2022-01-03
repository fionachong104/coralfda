#can we get Regina precipitation example from Ramsay et al to work?
rm(list = ls())
library(fda)

#get pls first axis scores from areascores df in RMD, alphabetically
#areascores$areax[order(areascores$Site)]

#import data
oneyeardf <- read.csv("C:/Users/624225/Box/PhD/coralfda/oneyeardb.csv", 
                      row.names=1)
axisscores <- read.csv("C:/Users/624225/Box/PhD/coralfda/axisscores.csv", 
                       row.names=1)

axisscores <- as.numeric(axisscores[,2])

#area or logarea
SortedArea <- sort(oneyeardf$Area)

 #from description on p. 75
#N <- length(RegPrec) #this is a guess: could it be where we are wrong? But the knot values in the next line match their figure 5.6, I think
#Wknots = RegPrec[round(N * seq(1 / N, 1, len = 11), 0)]
#Wnbasis = length(Wknots) + 2
#Wbasis = create.bspline.basis(range(RegPrec), 13, 4, Wknots)
#Wlambda = 1e-1
#WfdPar = fdPar(Wbasis, 2, Wlambda)
#densityList = density.fd(RegPrec, WfdPar) #same error about number of items to replace is not a multiple of replacement length as in our coral attempt
#Wfd = densityList$Wfdobj
#C = densityList$C
#Zfine = seq(RegPrec[1],RegPrec[N],len=201)
#Wfine = eval.fd(Zfine, Wfd)
#Pfine = exp(Wfine)/C
nbasis <- 4
#this is adapted from the example in ?density.fd(), but using the Regina precipitation data as in the book
basisobj <- create.bspline.basis(range(SortedArea), nbasis)
#  set up initial value for Wfdobj
Wfd0 <- fd(matrix(0, nbasis, 1), basisobj)
WfdParobj <- fdPar(Wfd0)
#  estimate density for each site
denslist <- density.fd(SortedArea, WfdParobj)
responsedensity <- denslist[[1]]
responsedensity$coefs <- array(dim=c(nbasis,18))
sites <- unique(oneyeardf$Site)
nsites <- length(sites)
for (i in 1:nsites){
  sitesarea <- sort(oneyeardf$Area[oneyeardf$Site == sites[i]])
  sitesdenslist <- density.fd(sitesarea, WfdParobj)
  responsedensity$coefs[,i] <- sitesdenslist[[1]]$coefs
  print(sitesdenslist[[1]]$basis)
}
# randomly generated site scores, need to use real scores
# testscores <- rnorm(n = 18)
SortedArea.f <- fRegress(responsedensity ~ axisscores)

scoresfit <- SortedArea.f$yhatfdobj
plot(scoresfit)

#effect of increasing pls test score (getting colder)
plot(SortedArea.f$betaestlist[[2]])

# think about smoothing - chapter 5 of book. 
# fourier spline function is for periodic data e.g. annual temp.
# chapter 10.1.3 on smoothing parameters

# 
# #  plot density
# xval <- seq(from = min(SortedArea), to = max(SortedArea), length.out = 100)
# wval <- eval.fd(xval, denslist$Wfdobj)
# pval <- exp(wval) / denslist$C
# plot(xval, pval, type="l")
# points(FR_SortedArea, rep(0, length(FR_SortedArea)))


