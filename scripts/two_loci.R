## Simulate different threshold and resistance outcomes
# Release threshold reflects how intense your release is
# Resistance allele formation rate (rM) controls how often the intended drive
  # allele (M) mutates into a resistance allele (R)

####################
# Load libraries
####################
library(MGDrivE)
source("cubes/cube_MEREA_two_loci.R")  # Use the two-locus cube with parameters
source("cubes/cube_auxiliary.R")

current_run <- "mgdrive/two_loci"
dir.create(current_run)

####################
# Define parameter ranges
####################
introduction_thresholds <- seq(0.1, 1.0, by = 0.1)  # From 10% to 100%
resistance_values <- seq(0, 0.5, by = 0.1)  # Resistance allele probability

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
for (threshold in introduction_thresholds) {
  for (rM_value in resistance_values) {
    
    print(paste("Running: Threshold =", threshold, "rM =", rM_value))
    
    # Create inheritance cube with current rM
    cube <- cubeMEREA_2L(rM = rM_value)
    
    # Calculate number of individuals to release
    release_size <- adultPopEquilibrium * threshold
    
    # Release parameters
    releasesParameters <- list(releasesStart=100, releasesNumber=3,
                               releasesInterval=30, releaseProportion=release_size)
    
    # Generate release vectors
    maleReleasesVector <- generateReleaseVector(
      driveCube = cube,
      nameGenotypes = list(
        c("MaMb", release_size / 3),
        c("MaMa", release_size / 3),
        c("MbMb", release_size / 3)
      ),
      releasesParameters = releasesParameters
    )
    
    femaleReleasesVector <- generateReleaseVector(
      driveCube = cube,
      nameGenotypes = list(
        c("MaW", release_size / 2),
        c("MbW", release_size / 2)
      ),
      releasesParameters = releasesParameters
    )
    
    # Patch setup
    patchReleases <- replicate(n=sitesNumber,
                               expr={list(maleReleases=NULL, femaleReleases=NULL,
                                          eggReleases=NULL, matedFemaleReleases=NULL)},
                               simplify=FALSE)
    
    patchReleases[[1]]$maleReleases <- maleReleasesVector
    patchReleases[[1]]$femaleReleases <- femaleReleasesVector
    
    # Unique output folder
    outFolder <- paste0(current_run, "/threshold_", threshold, "_rM_", rM_value)
    dir.create(outFolder, recursive = TRUE, showWarnings = FALSE)
    
    # MGDrivE parameters
    netPar <- parameterizeMGDrivE(
      runID=paste0("threshold_", threshold, "_rM_", rM_value), simTime=tMax,
      sampTime=30, nPatch=sitesNumber, beta=bioParameters$betaK, muAd=bioParameters$muAd,
      popGrowth=bioParameters$popGrowth, tEgg=bioParameters$tEgg,
      tLarva=bioParameters$tLarva, tPupa=bioParameters$tPupa,
      AdPopEQ=adultPopEquilibrium, inheritanceCube=cube
    )
    
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
    plot_file <- file.path(outFolder, paste0("plot_threshold_", threshold, "_rM_", rM_value, ".png"))
    png(filename = plot_file, width = 1200, height = 800)
    plotMGDrivESingle(readDir = outFolder, totalPop = TRUE, lwd = 3.5, alpha = 1)
    dev.off()
    
    # Display the plot in RStudio Viewer
    plotMGDrivESingle(readDir = outFolder, totalPop = TRUE, lwd = 3.5, alpha = 1)
    
    print(paste("Completed: Threshold =", threshold, "rM =", rM_value))
  }
}

## Simulation end message
print("All simulations complete.")