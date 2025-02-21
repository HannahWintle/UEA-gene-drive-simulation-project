## testing mgdrive====
####################
# Load libraries
####################
library(MGDrivE)
#> Loading MGDrivE: Mosquito Gene Drive Explorer

####################
# Output Folder
####################
dir.create("mgdrive")
outFolder <- "mgdrive/run005"

####################
# Simulation Parameters
####################
# days to run the simulation
tMax <- 25

# entomological parameters
bioParameters <- list(betaK=20, tEgg=5, tLarva=6, tPupa=4, popGrowth=1.175, muAd=0.09)

# a 1-node network where mosquitoes do not leave
moveMat <- matrix(data = 1, nrow = 1, ncol = 1)

# parameters of the population equilibrium
adultPopEquilibrium <- 500
sitesNumber <- nrow(moveMat)

####################
# Basic Inheritance pattern
####################
# Mendelian cube with standard (default) parameters
cube <- cubeMEDEA()

## Safety checks
if(any(c(Teff) < 0) || any(c(Teff) > 1)){
  stop("Teff must be a probability between 0 and 1")
}

####################
# Setup releases and batch migration
####################

patchReleases <- replicate(n=sitesNumber,
                           expr={list(maleReleases=NULL,femaleReleases=NULL,
                                      eggReleases=NULL,matedFemaleReleases=NULL)},
                           simplify=FALSE)

releasesParameters <- list(releasesStart=0,  # Release starts at generation 0
                           releasesNumber=1,  # Single release event
                           releasesInterval=0,
                           releaseProportion=0.1)

releasesVector <- generateReleaseVector(driveCube=cube,
                                        releasesParameters=releasesParameters)

patchReleases[[1]]$maleReleases <- releasesVector
patchReleases[[1]]$femaleReleases <- releasesVector

# batch migration is disabled by setting the probability to 0
# This is required because of the stochastic simulations, but doesn't make sense
#  in a deterministic simulation.
batchMigration <- basicBatchMigration(batchProbs=0,
                                      sexProbs=c(.5,.5), numPatches=sitesNumber)

                                      
####################
## Define matrices
####################
## Matrix Dimensions Key: [femaleGenotype, maleGenotype, offspringGenotype]
gtype <- c('TT', 'Tt', 'tt')
size <- length(gtype)
tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) # transition matrix

## Fill tMatrix with probabilities based on Mendelian inheritance
# TT x TT -> 100% TT
tMatrix['TT','TT', 'TT'] <- 1

# TT x Tt -> 50% TT, 50% Tt
tMatrix['TT','Tt', c('TT', 'Tt')] <- c(0.5, 0.5)

# TT x tt -> 100% Tt
tMatrix['TT','tt', 'Tt'] <- 1

# Tt x TT -> 50% TT, 50% Tt
tMatrix['Tt','TT', c('TT', 'Tt')] <- c(0.5, 0.5)

# Tt x Tt -> 25% TT, 50% Tt, 25% tt
tMatrix['Tt','Tt', c('TT', 'Tt', 'tt')] <- c(0.25, 0.5, 0.25)

# Tt x tt -> 50% Tt, 50% tt
tMatrix['Tt','tt', c('Tt', 'tt')] <- c(0.5, 0.5)

# tt x TT -> 100% Tt
tMatrix['tt','TT', 'Tt'] <- 1

# tt x Tt -> 50% Tt, 50% tt
tMatrix['tt','Tt', c('Tt', 'tt')] <- c(0.5, 0.5)

# tt x tt -> 100% tt
tMatrix['tt','tt', 'tt'] <- 1

## Set the other half of the matrix (ensures symmetry)
boolMat <- upper.tri(x = tMatrix[ , ,1], diag = FALSE)
for(z in 1:size) {
  tMatrix[ , ,z][boolMat] <- t(tMatrix[ , ,z])[boolMat]
}

tMatrix[tMatrix < .Machine$double.eps] <- 0 

viabilityMask <- array(data = 1, dim = c(size, size, size), dimnames = list(gtype, gtype, gtype))

## Adjust viability based on genotype interactions
for(slice in 1:size) {
  viabilityMask['Tt',slice, ] <- c(1, 1, 1)
  viabilityMask['tt',slice, ] <- c(1, 1, 1)
}

####################
# Combine parameters and run!
####################
# set MGDrivE to run deterministic
setupMGDrivE(stochasticityON = FALSE, verbose = FALSE)

# setup parameters for the network. This builds a list of parameters required for
#  every population in the network. In ths case, we havee a network of 1 population.
netPar <- parameterizeMGDrivE(runID=1, simTime=tMax, sampTime = 1, nPatch=sitesNumber,
                              beta=bioParameters$betaK, muAd=bioParameters$muAd,
                              popGrowth=bioParameters$popGrowth, tEgg=bioParameters$tEgg,
                              tLarva=bioParameters$tLarva, tPupa=bioParameters$tPupa,
                              AdPopEQ=adultPopEquilibrium, inheritanceCube = cube)

# build network prior to run
MGDrivESim <- Network$new(params=netPar,
                          driveCube=cube,
                          patchReleases=patchReleases,
                          migrationMale=moveMat,
                          migrationFemale=moveMat,
                          migrationBatch=batchMigration,
                          directory=outFolder,
                          verbose=FALSE)

# run simulation
MGDrivESim$oneRun(verbose = T)

####################
# Post Analysis
####################
# split output by patch
#  Required for plotting later
splitOutput(readDir = outFolder, remFile = TRUE, verbose = FALSE)

# aggregate females by their mate choice
#  This reduces the female file to have the same columns as the male file
aggregateFemales(readDir = outFolder, genotypes = cube$genotypesID,
                 remFile = TRUE, verbose = FALSE)

# plot output to see effect
plotMGDrivESingle(readDir = outFolder, totalPop = TRUE, lwd = 3.5, alpha = 1)
