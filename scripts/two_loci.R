## Simulate different threshold and resistance outcomes
# Release threshold reflects how intense your release is
# Resistance allele formation rate (rM) controls how often the intended drive
  # allele (M) mutates into a resistance allele (R)

####################
# Load libraries
####################
library(MGDrivE)
source("cubes/cube_MEREA_two_loci.R")
source("cubes/cube_auxiliary.R")

current_run <- "mgdrive/two_loci_test"
dir.create(current_run)

####################
# Define parameter ranges
####################
introduction_thresholds <- seq(0.1, 1.0, by = 0.1)  # From 10% to 100%
resistance_values <- seq(0, 0.5, by = 0.1)  # Resistance allele probability

data <- expand.grid(introduction_thresholds = introduction_thresholds, 
            resistance_values = resistance_values)

run <- nrow(data)

####################
# Simulation Parameters
####################
tMax <- 1000  # Total simulation time
bioParameters <- list(betaK=20, tEgg=5, tLarva=6, tPupa=4, popGrowth=1.175, muAd=0.09)
moveMat <- matrix(data = 1, nrow = 1, ncol = 1)
adultPopEquilibrium <- 500
sitesNumber <- nrow(moveMat)

####################
# Nested loop over resistance and release frequency
####################
for (i in 1:nrow(data)) {
    
    print(paste("Running: Threshold =", data$introduction_thresholds[[i]], "rM =", data$resistance_values[[i]]))
    
    # Create inheritance cube with current rM
    cube <- cubeMEREA_2L(rM = data$resistance_values[[i]]
    )
    
    # Release parameters
    releasesParameters <- list(releasesStart=0, releasesNumber=0,
                               releasesInterval=30, releaseProportion=0)
    
    # Generate release vectors
    maleReleasesVector <- NULL
    
    femaleReleasesVector <- NULL
    # Patch setup
    patchReleases <- replicate(n=sitesNumber,
                               expr={list(maleReleases=NULL, femaleReleases=NULL,
                                          eggReleases=NULL, matedFemaleReleases=NULL)},
                               simplify=FALSE)
    
   # patchReleases[[1]]$maleReleases <- maleReleasesVector
  #  patchReleases[[1]]$femaleReleases <- femaleReleasesVector
    
    # Unique output folder
    outFolder <- paste0(current_run, "/threshold_", data$introduction_thresholds[[i]], "_rM_", data$resistance_values[[i]])
    dir.create(outFolder, recursive = TRUE, showWarnings = FALSE)
    
    # MGDrivE parameters
    netPar <- parameterizeMGDrivE(
      runID=paste0("threshold_", data$introduction_thresholds[[i]], "_rM_", data$resistance_values[[i]]), simTime=tMax,
      sampTime=30, nPatch=sitesNumber, beta=bioParameters$betaK, muAd=bioParameters$muAd,
      popGrowth=bioParameters$popGrowth, tEgg=bioParameters$tEgg,
      tLarva=bioParameters$tLarva, tPupa=bioParameters$tPupa,
      AdPopEQ=adultPopEquilibrium, inheritanceCube=cube
    )
    
    netPar$AdPopRatio_F <- matrix(c(1), nrow = 1, dimnames = list(NULL, c("ZW")))
    
    netPar$AdPopRatio_M <- matrix(c(1-data$introduction_thresholds[[i]], data$introduction_thresholds[[i]]), nrow = 1, dimnames = list(NULL, c("ZZ", "MaMb")))
    
    netPar$LarPopRatio <- NULL
    
    # Initialize and run simulation
    MGDrivESim <- Network$new(
      params=netPar,
      driveCube=cube,
      patchReleases=patchReleases,
      migrationMale=moveMat,
      migrationFemale=moveMat,
      migrationBatch=basicBatchMigration(batchProbs=0, sexProbs=c(.5,.5), numPatches=sitesNumber),
      directory=outFolder,
      verbose=FALSE
    )
    
    MGDrivESim$oneRun(verbose = TRUE)
    ####################
    # Post-processing and plotting
    ####################
    splitOutput(readDir = outFolder, remFile = TRUE, verbose = FALSE)
    aggregateFemales(readDir = outFolder, genotypes = cube$genotypesID,
                     remFile = TRUE, verbose = FALSE)
    
    # Plot results
    plot_file <- file.path(outFolder, paste0("plot_threshold_", data$introduction_thresholds[[i]], "_rM_", data$resistance_values[[i]], ".png"))
    png(filename = plot_file, width = 1200, height = 800)
    plotMGDrivESingle(readDir = outFolder, totalPop = TRUE, lwd = 3.5, alpha = 1)
    dev.off()
    
    # Display the plot in RStudio Viewer
    plotMGDrivESingle(readDir = outFolder, totalPop = TRUE, lwd = 3.5, alpha = 1)
    
    print(paste("Completed: Threshold =", data$introduction_thresholds[[i]], "rM =", data$resistance_values[[i]]))
  }


## Simulation end message
print("All simulations complete.")