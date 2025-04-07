## Testing MEREA in MGDrivE
####################
# Load libraries
####################
library(MGDrivE)

####################
# Output Folder
####################
dir.create("mgdrive")
outFolder <- "mgdrive/run008"

####################
# Simulation Parameters
####################
tMax <- 365  # Simulation run time

# Entomological parameters (unchanged)
bioParameters <- list(betaK=20, tEgg=5, tLarva=6, tPupa=4, popGrowth=1.175, muAd=0.09)

# Single-node population model (no migration)
moveMat <- matrix(data = 1, nrow = 1, ncol = 1)

# Population equilibrium size
adultPopEquilibrium <- 500
sitesNumber <- nrow(moveMat)

####################
# Basic Inheritance pattern
####################
# Using the MEREA cube instead of MEDEA
cube <- cubeMEREA()

####################
# Setup releases and batch migration
####################

# Initialize release settings
patchReleases <- replicate(n=sitesNumber,
                           expr={list(maleReleases=NULL, femaleReleases=NULL,
                                      eggReleases=NULL, matedFemaleReleases=NULL)},
                           simplify=FALSE)

# Define separate release proportions for male and female
releasesParametersMale <- list(releasesStart=0, # Release starts at generation 0
                               releasesNumber=1, # Single release event 
                               releasesInterval=0,  
                               releaseProportion=250) #Can tweak male and female releases separately #Adjust threshold testing here

releasesParametersFemale <- list(releasesStart=0,  
                                 releasesNumber=1,  
                                 releasesInterval=0,  
                                 releaseProportion=250)  # female release

# Generate release vectors for both sexes
maleReleasesVector <- generateReleaseVector(driveCube=cube, releasesParameters=releasesParametersMale)
femaleReleasesVector <- generateReleaseVector(driveCube=cube, releasesParameters=releasesParametersFemale)

# Assign release vectors to appropriate sex
patchReleases[[1]]$maleReleases <- maleReleasesVector
patchReleases[[1]]$femaleReleases <- femaleReleasesVector

# batch migration is disabled by setting the probability to 0
# This is required because of the stochastic simulations, but doesn't make sense
#  in a deterministic simulation.
batchMigration <- basicBatchMigration(batchProbs=0,
                                      sexProbs=c(.5,.5), numPatches=sitesNumber)

####################
# Combine parameters and run!
####################
# Set MGDrivE to deterministic mode
setupMGDrivE(stochasticityON = FALSE, verbose = FALSE)

# Define network parameters
netPar <- parameterizeMGDrivE(runID=1, simTime=tMax, sampTime=1, nPatch=sitesNumber,
                              beta=bioParameters$betaK, muAd=bioParameters$muAd,
                              popGrowth=bioParameters$popGrowth, tEgg=bioParameters$tEgg,
                              tLarva=bioParameters$tLarva, tPupa=bioParameters$tPupa,
                              AdPopEQ=adultPopEquilibrium, inheritanceCube=cube)

# Initialize simulation network
MGDrivESim <- Network$new(params=netPar,
                          driveCube=cube,
                          patchReleases=patchReleases,
                          migrationMale=moveMat,
                          migrationFemale=moveMat,
                          migrationBatch=batchMigration,
                          directory=outFolder,
                          verbose=FALSE)

# Run simulation
MGDrivESim$oneRun(verbose = TRUE)

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
