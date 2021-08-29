
#setwd("Choose your working directory")

library(R2WinBUGS)
library(spdep)
library(maps)
library(rgdal)
library(reshape)

#read data objects
load("DetNonDet_Data.RData")
load("COVSdata.RData")

ni = 750000
nb = 350000
nthin = 250
nc = 3


###############################################
### All data + effort + site COVS + CAR

data <- list(muni1 = dat$IDENT, muni2 = dat2$IDENT, muni3 = dat3$IDENT, muni4 = dat4$IDENT, 
              Y1 = dat$A_VINACEA, Y2 = dat2$A_VINACEA, Y3 = dat3$AVINACEA, Y4 = dat4$A_VINACEA, 
              nMuni = length(muni), nObs1 = nrow(dat), nObs2 = nrow(dat2), 
              nObs3 = nrow(dat3), nObs4 = nrow(dat4), 
              TObs = dat$EFFORT_MIN/60, SSee = dat2$NSPECIES, TObs2 = dat2$DURATION.M/60, 
              RLen = dat2$EFFORT.DIS, NPho = dat3$NPIC, NAud = dat3$NSONG,
              NAud2 = dat4$NSONGS, VegCover = VegCover, ArauCover = ArauCover, 
              Alt = Altitude, nCell = nrow(hex.centroids), cell.id = cell.id, 
              nMuni = length(muni), adj = adj, num = num, sumNeigh = sumNeigh) 

inits = function() {list(z = rep(1, data$nMuni))}#, alpha=rep(0, nrow(hex.centroids)))}#, beta = c(0.3), BETA = rep(0.01,7), alpha = rep(0,data$nMuni))} 
params <- c("beta","psi","BETA","z", "spacesigma","alpha", "muP1", "muP2", "muP3", "muP4")

mod = "path.../FullModelWinbugs_Code.txt"
fullModel <- bugs(data = data, inits = inits, parameters.to.save = params, model.file = mod,
             n.chains = nc, n.iter = ni,n.burnin = nb, n.thin = nthin, debug = TRUE)

save(fullModel, file = "outFULLModel.RData")


