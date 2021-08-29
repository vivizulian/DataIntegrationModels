
setwd("choose working directory")

library(R2WinBUGS)
library(spdep)
library(maps)
library(rgdal)
library(reshape)

load("DevianceDataCalculation.RData")
load("COVSdata.RData")

ni = 200000
nb = 150000
nthin = 100
nc = 3

#############################################
# 20% of the data Out - M1

dataM1 <- list(muni1 = datIn$IDENT, muni2 = datOut$IDENT, muni3 = dat2In$IDENT, muni4 = dat2Out$IDENT, 
              muni5 = dat3In$IDENT, muni6 = dat3Out$IDENT, muni7 = dat4In$IDENT, muni8 = dat4Out$IDENT,
              Y1 = datIn$A_VINACEA, Y3 = dat2In$A_VINACEA, Y5 = dat3In$AVINACEA, Y7 = dat4In$A_VINACEA, 
              nMuni = length(muni), nObs1 = nrow(datIn), nObs2 = nrow(datOut), nObs3 = nrow(dat2In), 
              nObs4 = nrow(dat2Out), nObs5 = nrow(dat3In), nObs6 = nrow(dat3Out), nObs7 = nrow(dat4In), 
              nObs8 = nrow(dat4Out), Ytruth2 = datOut$A_VINACEA, Ytruth4 = dat2Out$A_VINACEA, 
              Ytruth6 = dat3Out$AVINACEA, Ytruth8 = dat4Out$A_VINACEA,
              TObs = datIn$EFFORT_MIN/60, TObs2 = datOut$EFFORT_MIN/60, SSee = dat2In$NSPECIES, 
              SSee2 = dat2Out$NSPECIES, TObs3 = dat2In$DURATION.M/60, TObs4 = dat2Out$DURATION.M/60, 
              RLen = dat2In$EFFORT.DIS, RLen2 = dat2Out$EFFORT.DIS, NPho = dat3In$NPIC, NPho2 = dat3Out$NPIC, 
              NAud = dat3In$NSONG, NAud2 = dat3Out$NSONG, NAud3 = dat4In$NSONGS, NAud4 = dat4Out$NSONGS, 
              VegCover = VegCover, ArauCover = ArauCover, Alt = Altitude, nCell = nrow(hex.centroids), 
              cell.id = cell.id, adj = adj, num = num, sumNeigh = sumNeigh)

inits = function() {list(z = rep(1, dataM1$nMuni))}#, alpha=rep(0, nrow(hex.centroids)))}#, beta = c(0.3), BETA = rep(0.01,7), alpha = rep(0,data$nMuni))} 
params <- c("beta","psi","z", "BETA", "muP1", "muP2", "muP3", "muP4", "muP5", "muP6", "muP7", "muP8", 
            "Y2", "Y4", "Y6", "Y8", "spacesigma", "alpha")

mod4 = "path...ValidationDataInOut.txt"
outDataInOutM1 <- bugs(data = dataM1, inits = inits, parameters.to.save = params, model.file = mod4,
                           n.chains = nc, n.iter = ni,n.burnin = nb, n.thin = nthin, debug = FALSE)


#### Compute fit scores ####
# Data set
Ytruth2 <- datOut$A_VINACEA 
Yhat2 <- outDataInOutM1$mean$Y2

Ytruth4 <- dat2Out$A_VINACEA
Yhat4 <- outDataInOutM1$mean$Y4

Ytruth6 <- dat3Out$AVINACEA
Yhat6 <- outDataInOutM1$mean$Y6

Ytruth8 <- dat4Out$A_VINACEA
Yhat8 <- outDataInOutM1$mean$Y8

## Deviance:
likhood2 <- (Yhat2^Ytruth2)*((1-Yhat2)^(1-Ytruth2))
DEV2 <- -(2*(sum(log(likhood2))))

likhood4 <- (Yhat4^Ytruth4)*((1-Yhat4)^(1-Ytruth4))
DEV4 <- -(2*(sum(log(likhood4))))

likhood6 <- (Yhat6^Ytruth6)*((1-Yhat6)^(1-Ytruth6))
DEV6 <- -(2*(sum(log(likhood6))))

likhood8 <- (Yhat8^Ytruth8)*((1-Yhat8)^(1-Ytruth8))
DEV8 <- -(2*(sum(log(likhood8))))

DEVtotal <- DEV2 + DEV4 + DEV6 + DEV8


#############################################
### Model 2 = Taking out the effort equation and 20% of the munis data:

data <- list(muni1 = datIn$IDENT, muni2 = datOut$IDENT, muni3 = dat2In$IDENT, muni4 = dat2Out$IDENT, 
             muni5 = dat3In$IDENT, muni6 = dat3Out$IDENT, muni7 = dat4In$IDENT, muni8 = dat4Out$IDENT,
             Y1 = datIn$A_VINACEA, Y3 = dat2In$A_VINACEA, Y5 = dat3In$AVINACEA, Y7 = dat4In$A_VINACEA, 
             nMuni = length(muni), nObs1 = nrow(datIn), nObs2 = nrow(datOut), nObs3 = nrow(dat2In), 
             nObs4 = nrow(dat2Out), nObs5 = nrow(dat3In), nObs6 = nrow(dat3Out), nObs7 = nrow(dat4In), 
             nObs8 = nrow(dat4Out), Ytruth2 = datOut$A_VINACEA, Ytruth4 = dat2Out$A_VINACEA, 
             Ytruth6 = dat3Out$AVINACEA, Ytruth8 = dat4Out$A_VINACEA, VegCover = VegCover, 
             ArauCover = ArauCover, Alt = Altitude, nCell = nrow(hex.centroids), cell.id = cell.id, 
             adj = adj, num = num, sumNeigh = sumNeigh)
#TObs = dat$EFFORT_MIN/60, SSee = dat2$NSPECIES, TObs2 = dat2$DURATION.M/60, 
#RLen = dat2$EFFORT.DIS, NPho = dat3$NPIC, NAud = dat3$NSONG,
#NAud2 = dat4$NSONGS, 

inits = function() {list(z = rep(1, data$nMuni))}#, alpha=rep(0, nrow(hex.centroids)))}#, beta = c(0.3), BETA = rep(0.01,7), alpha = rep(0,data$nMuni))} 
params <- c("beta","psi","z", "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", 
            "Y2", "Y4", "Y6", "Y8", "BRSC2", "BRSC4", "BRSC6", "BRSC8", "spacesigma", "alpha")

mod2 = "path...ValidationEffortOut.txt"
outDataEffortM1 <- bugs(data = data, inits = inits, parameters.to.save = params, model.file = mod2,
                           n.chains = nc, n.iter = ni,n.burnin = nb, n.thin = nthin, debug = FALSE)


#### Compute fit scores ####
# Data set
Ytruth2 <- datOut$A_VINACEA 
Yhat2 <- outDataEffortM1$mean$Y2

Ytruth4 <- dat2Out$A_VINACEA
Yhat4 <- outDataEffortM1$mean$Y4

Ytruth6 <- dat3Out$AVINACEA
Yhat6 <- outDataEffortM1$mean$Y6

Ytruth8 <- dat4Out$A_VINACEA
Yhat8 <- outDataEffortM1$mean$Y8

## Deviance:
likhood2 <- (Yhat2^Ytruth2)*((1-Yhat2)^(1-Ytruth2))
DEV2 <- -(2*(sum(log(likhood2))))

likhood4 <- (Yhat4^Ytruth4)*((1-Yhat4)^(1-Ytruth4))
DEV4 <- -(2*(sum(log(likhood4))))

likhood6 <- (Yhat6^Ytruth6)*((1-Yhat6)^(1-Ytruth6))
DEV6 <- -(2*(sum(log(likhood6))))

likhood8 <- (Yhat8^Ytruth8)*((1-Yhat8)^(1-Ytruth8))
DEV8 <- -(2*(sum(log(likhood8))))

DEVtotal <- DEV2 + DEV4 + DEV6 + DEV8

###############################################
### 3. Taking out the CAR and the same 20% of the munis data:

data3 <- list(muni1 = datIn$IDENT, muni2 = datOut$IDENT, muni3 = dat2In$IDENT, muni4 = dat2Out$IDENT, 
              muni5 = dat3In$IDENT, muni6 = dat3Out$IDENT, muni7 = dat4In$IDENT, muni8 = dat4Out$IDENT,
              Y1 = datIn$A_VINACEA, Y3 = dat2In$A_VINACEA, Y5 = dat3In$AVINACEA, Y7 = dat4In$A_VINACEA, 
              nMuni = length(muni), nObs1 = nrow(datIn), nObs2 = nrow(datOut), nObs3 = nrow(dat2In), 
              nObs4 = nrow(dat2Out), nObs5 = nrow(dat3In), nObs6 = nrow(dat3Out), nObs7 = nrow(dat4In), 
              nObs8 = nrow(dat4Out), Ytruth2 = datOut$A_VINACEA, Ytruth4 = dat2Out$A_VINACEA, 
              Ytruth6 = dat3Out$AVINACEA, Ytruth8 = dat4Out$A_VINACEA,
              TObs = datIn$EFFORT_MIN/60, TObs2 = datOut$EFFORT_MIN/60, SSee = dat2In$NSPECIES, 
              SSee2 = dat2Out$NSPECIES, TObs3 = dat2In$DURATION.M/60, TObs4 = dat2Out$DURATION.M/60, 
              RLen = dat2In$EFFORT.DIS, RLen2 = dat2Out$EFFORT.DIS, NPho = dat3In$NPIC, NPho2 = dat3Out$NPIC, 
              NAud = dat3In$NSONG, NAud2 = dat3Out$NSONG, NAud3 = dat4In$NSONGS, NAud4 = dat4Out$NSONGS, 
              VegCover = VegCover, ArauCover = ArauCover, Alt = Altitude)#, nCell = nrow(hex.centroids), 
#cell.id = cell.id, adj = adj, num = num, sumNeigh = sumNeigh)

inits = function() {list(z = rep(1, data3$nMuni))}#, alpha=rep(0, nrow(hex.centroids)))}#, beta = c(0.3), BETA = rep(0.01,7), alpha = rep(0,data$nMuni))} 
params <- c("beta","psi","z", "BETA", "muP1", "muP2", "muP3", "muP4", "muP5", "muP6", "muP7", "muP8", 
            "Y2", "Y4", "Y6", "Y8", "BRSC2", "BRSC4", "BRSC6", "BRSC8")

mod3 = "path.../ValidationCAROut.txt"
outCAROutM1 <- bugs(data = data3, inits = inits, parameters.to.save = params, model.file = mod3,
                    n.chains = nc, n.iter = ni,n.burnin = nb, n.thin = nthin, debug = FALSE)


#### Compute fit scores ####
# Data set
Ytruth2 <- datOut$A_VINACEA 
Yhat2 <- outCAROutM1$mean$Y2

Ytruth4 <- dat2Out$A_VINACEA
Yhat4 <- outCAROutM1$mean$Y4

Ytruth6 <- dat3Out$AVINACEA
Yhat6 <- outCAROutM1$mean$Y6

Ytruth8 <- dat4Out$A_VINACEA
Yhat8 <- outCAROutM1$mean$Y8

## Deviance:
likhood2 <- (Yhat2^Ytruth2)*((1-Yhat2)^(1-Ytruth2))
DEV2 <- -(2*(sum(log(likhood2))))

likhood4 <- (Yhat4^Ytruth4)*((1-Yhat4)^(1-Ytruth4))
DEV4 <- -(2*(sum(log(likhood4))))

likhood6 <- (Yhat6^Ytruth6)*((1-Yhat6)^(1-Ytruth6))
DEV6 <- -(2*(sum(log(likhood6))))

likhood8 <- (Yhat8^Ytruth8)*((1-Yhat8)^(1-Ytruth8))
DEV8 <- -(2*(sum(log(likhood8))))

DEVtotal <- DEV2 + DEV4 + DEV6 + DEV8


###############################################
### 4. Taking out the COVS:

data6 <- list(muni1 = datIn$IDENT, muni2 = datOut$IDENT, muni3 = dat2In$IDENT, muni4 = dat2Out$IDENT, 
              muni5 = dat3In$IDENT, muni6 = dat3Out$IDENT, muni7 = dat4In$IDENT, muni8 = dat4Out$IDENT,
              Y1 = datIn$A_VINACEA, Y3 = dat2In$A_VINACEA, Y5 = dat3In$AVINACEA, Y7 = dat4In$A_VINACEA, 
              nMuni = length(muni), nObs1 = nrow(datIn), nObs2 = nrow(datOut), nObs3 = nrow(dat2In), 
              nObs4 = nrow(dat2Out), nObs5 = nrow(dat3In), nObs6 = nrow(dat3Out), nObs7 = nrow(dat4In), 
              nObs8 = nrow(dat4Out), Ytruth2 = datOut$A_VINACEA, Ytruth4 = dat2Out$A_VINACEA, 
              Ytruth6 = dat3Out$AVINACEA, Ytruth8 = dat4Out$A_VINACEA,
              TObs = datIn$EFFORT_MIN/60, TObs2 = datOut$EFFORT_MIN/60, SSee2 = dat2Out$NSPECIES, SSee = dat2In$NSPECIES, 
              TObs3 = dat2In$DURATION.M/60, TObs4 = dat2Out$DURATION.M/60, 
              RLen = dat2In$EFFORT.DIS, RLen2 = dat2Out$EFFORT.DIS, NPho = dat3In$NPIC, NPho2 = dat3Out$NPIC, 
              NAud = dat3In$NSONG, NAud2 = dat3Out$NSONG, NAud3 = dat4In$NSONGS, NAud4 = dat4Out$NSONGS, 
              nCell = nrow(hex.centroids), cell.id = cell.id, adj = adj, num = num, sumNeigh = sumNeigh) 
#VegCover = VegCover, ArauCover = ArauCover, Alt = Altitude, 

inits = function() {list(z = rep(1, data6$nMuni))}#, alpha=rep(0, nrow(hex.centroids)))}#, beta = c(0.3), BETA = rep(0.01,7), alpha = rep(0,data$nMuni))} 
params <- c("beta","psi","z", "BETA", "muP1", "muP2", "muP3", "muP4", "muP5", "muP6", "muP7", "muP8", 
            "Y2", "Y4", "Y6", "Y8", "BRSC2", "BRSC4", "BRSC6", "BRSC8", "spacesigma", "alpha")

mod6 = "path.../ValidationCOVSOut.txt"
outCOVSOutM1 <- bugs(data = data6, inits = inits, parameters.to.save = params, model.file = mod6,
                             n.chains = nc, n.iter = ni,n.burnin = nb, n.thin = nthin, debug = FALSE)


#### Compute fit scores ####
# Data set
Ytruth2 <- datOut$A_VINACEA 
Yhat2 <- outCOVSOutM1$mean$Y2

Ytruth4 <- dat2Out$A_VINACEA
Yhat4 <- outCOVSOutM1$mean$Y4

Ytruth6 <- dat3Out$AVINACEA
Yhat6 <- outCOVSOutM1$mean$Y6

Ytruth8 <- dat4Out$A_VINACEA
Yhat8 <- outCOVSOutM1$mean$Y8

## Deviance:
likhood2 <- (Yhat2^Ytruth2)*((1-Yhat2)^(1-Ytruth2))
DEV2 <- -(2*(sum(log(likhood2))))

likhood4 <- (Yhat4^Ytruth4)*((1-Yhat4)^(1-Ytruth4))
DEV4 <- -(2*(sum(log(likhood4))))

likhood6 <- (Yhat6^Ytruth6)*((1-Yhat6)^(1-Ytruth6))
DEV6 <- -(2*(sum(log(likhood6))))

likhood8 <- (Yhat8^Ytruth8)*((1-Yhat8)^(1-Ytruth8))
DEV8 <- -(2*(sum(log(likhood8))))

DEVtotal <- DEV2 + DEV4 + DEV6 + DEV8




## Individual datasets ##
#############################################
### Model 1 = Only RC data

data <- list(muni1 = datIn$IDENT, muni2 = datOut$IDENT, Y1 = datIn$A_VINACEA,
             nMuni = length(muni), nObs1 = nrow(datIn), nObs2 = nrow(datOut), 
             VegCover = VegCover, ArauCover = ArauCover, Alt = Altitude, nCell = nrow(hex.centroids), 
             cell.id = cell.id, adj = adj, num = num, sumNeigh = sumNeigh,
             TObs = datIn$EFFORT_MIN/60, TObs2 = datOut$EFFORT_MIN/60)

inits = function() {list(z = rep(1, data$nMuni))}#, alpha=rep(0, nrow(hex.centroids)))}#, beta = c(0.3), BETA = rep(0.01,7), alpha = rep(0,data$nMuni))} 
params <- c("beta", "BETA", "psi", "z", "muP1", "muP2", "Y2", "spacesigma", "alpha")

mod1 = "path.../COVSCAR_RC.txt"
outRC <- bugs(data = data, inits = inits, parameters.to.save = params, model.file = mod1,
                           n.chains = nc, n.iter = ni,n.burnin = nb, n.thin = nthin, debug = TRUE)


#### Compute fit scores ####
# Data set
Ytruth2 <- datOut$A_VINACEA 
Yhat2 <- outRC$mean$Y2

## Deviance:
likhood2 <- (Yhat2^Ytruth2)*((1-Yhat2)^(1-Ytruth2))
DEV2 <- -(2*(sum(log(likhood2))))



#############################################
### Model 2 = Only EB data

data2 <- list(muni3 = dat2In$IDENT, muni4 = dat2Out$IDENT, Y3 = dat2In$A_VINACEA, 
             nMuni = length(muni), nObs3 = nrow(dat2In), 
             nObs4 = nrow(dat2Out), Ytruth4 = dat2Out$A_VINACEA, SSee = dat2In$NSPECIES, 
             SSee2 = dat2Out$NSPECIES, TObs3 = dat2In$DURATION.M/60, TObs4 = dat2Out$DURATION.M/60, 
             RLen = dat2In$EFFORT.DIS, RLen2 = dat2Out$EFFORT.DIS,
             VegCover = VegCover, ArauCover = ArauCover, Alt = Altitude, nCell = nrow(hex.centroids), 
             cell.id = cell.id, adj = adj, num = num, sumNeigh = sumNeigh)

inits = function() {list(z = rep(1, data2$nMuni))}#, alpha=rep(0, nrow(hex.centroids)))}#, beta = c(0.3), BETA = rep(0.01,7), alpha = rep(0,data$nMuni))} 
params <- c("beta", "BETA", "psi", "z", "muP3", "muP4", "Y4", "spacesigma", "alpha")

mod2 = "path.../COVSCAR_EB.txt"
outEB <- bugs(data = data2, inits = inits, parameters.to.save = params, model.file = mod2,
              n.chains = nc, n.iter = ni,n.burnin = nb, n.thin = nthin, debug = TRUE)



#### Compute fit scores ####
# Data set
Ytruth4 <- dat2Out$A_VINACEA
Yhat4 <- outEB$mean$Y4

## Deviance:
likhood4 <- (Yhat4^Ytruth4)*((1-Yhat4)^(1-Ytruth4))
DEV4 <- -(2*(sum(log(likhood4))))



#############################################
### Model 3 = Only WA data

data3 <- list(muni5 = dat3In$IDENT, muni6 = dat3Out$IDENT, 
             Y5 = dat3In$AVINACEA, nMuni = length(muni), nObs5 = nrow(dat3In), nObs6 = nrow(dat3Out), 
             Ytruth6 = dat3Out$AVINACEA, NPho = dat3In$NPIC, NPho2 = dat3Out$NPIC, 
             NAud = dat3In$NSONG, NAud2 = dat3Out$NSONG,
             VegCover = VegCover, ArauCover = ArauCover, Alt = Altitude, nCell = nrow(hex.centroids), 
             cell.id = cell.id, adj = adj, num = num, sumNeigh = sumNeigh)

inits = function() {list(z = rep(1, data3$nMuni))}#, alpha=rep(0, nrow(hex.centroids)))}#, beta = c(0.3), BETA = rep(0.01,7), alpha = rep(0,data$nMuni))} 
params <- c("beta", "BETA", "psi", "z", "muP5", "muP6", "Y6", "spacesigma", "alpha")

mod3 = "path.../COVSCAR_WA.txt"
outWA <- bugs(data = data3, inits = inits, parameters.to.save = params, model.file = mod3,
                           n.chains = nc, n.iter = ni,n.burnin = nb, n.thin = nthin, debug = FALSE)

#### Compute fit scores ####
# Data set
Ytruth6 <- dat3Out$AVINACEA
Yhat6 <- outWA$mean$Y6

## Deviance:
likhood6 <- (Yhat6^Ytruth6)*((1-Yhat6)^(1-Ytruth6))
DEV6 <- -(2*(sum(log(likhood6))))


#############################################
### Model 4 = Only XC data

data4 <- list(muni7 = dat4In$IDENT, muni8 = dat4Out$IDENT, Y7 = dat4In$A_VINACEA, 
             nMuni = length(muni), nObs7 = nrow(dat4In), 
             nObs8 = nrow(dat4Out), Ytruth8 = dat4Out$A_VINACEA,
             NAud3 = dat4In$NSONGS, NAud4 = dat4Out$NSONGS, 
             VegCover = VegCover, ArauCover = ArauCover, Alt = Altitude, nCell = nrow(hex.centroids), 
             cell.id = cell.id, adj = adj, num = num, sumNeigh = sumNeigh)

inits = function() {list(z = rep(1, data4$nMuni))}#, alpha=rep(0, nrow(hex.centroids)))}#, beta = c(0.3), BETA = rep(0.01,7), alpha = rep(0,data$nMuni))} 
params <- c("beta", "BETA", "psi", "z", "muP7", "muP8", "Y8", "spacesigma", "alpha")

mod4 = "path.../COVSCAR_XC.txt"
outXC <- bugs(data = data4, inits = inits, parameters.to.save = params, model.file = mod4,
                           n.chains = nc, n.iter = ni,n.burnin = nb, n.thin = nthin, debug = FALSE)

#### Compute fit scores ####
# Data set
Ytruth8 <- dat4Out$A_VINACEA
Yhat8 <- outXC$mean$Y8

## Deviance:
likhood8 <- (Yhat8^Ytruth8)*((1-Yhat8)^(1-Ytruth8))
DEV8 <- -(2*(sum(log(likhood8))))


