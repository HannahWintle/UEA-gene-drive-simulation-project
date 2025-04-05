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

current_run <- "mgdrive/two_loci_test003"
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
# for loop
####################
data <- expand.grid(
  introduction_thresholds = introduction_thresholds,
  resistance_values = resistance_values
)

# Loop over each row using index i
for (i in 1:nrow(data)) {
  threshold <- data$introduction_thresholds[i]
  rM_value <- data$resistance_values[i]
  
  print(paste("Running: Threshold =", threshold, "rM =", rM_value))
  
    # Create inheritance cube
    cube <- cubeMEREA_2L(rM = rM_value)
    
    # Calculate release size
    release_size <- adultPopEquilibrium * threshold
    
    # Release setup
    releasesParameters <- list(
      releasesStart=100, 
      releasesNumber=3,
      releasesInterval=30, 
      releaseProportion=release_size
      )
    
    maleReleasesVector <- generateReleaseVector(
      driveCube = cube,
      nameGenotypes = list(c("MaMb", release_size)),
      releasesParameters = releasesParameters
    )
    
    femaleReleasesVector <- generateReleaseVector(
      driveCube = cube,
      nameGenotypes = list(c("MaW", release_size)),
      releasesParameters = releasesParameters
    )
    
    patchReleases <- replicate(n = sitesNumber,
                               expr = list(maleReleases = NULL, femaleReleases = NULL,
                                           eggReleases = NULL, matedFemaleReleases = NULL),
                               simplify = FALSE)
    patchReleases[[1]]$maleReleases <- maleReleasesVector
    patchReleases[[1]]$femaleReleases <- femaleReleasesVector
    
    outFolder <- paste0(current_run, "/threshold_", threshold, "_rM_", rM_value)
    dir.create(outFolder, recursive = TRUE, showWarnings = FALSE)

    netPar <- parameterizeMGDrivE(
      runID = paste0("threshold_", threshold, "_rM_", rM_value),
      simTime = tMax,
      sampTime = 30,
      nPatch = sitesNumber,
      beta = bioParameters$betaK,
      muAd = bioParameters$muAd,
      popGrowth = bioParameters$popGrowth,
      tEgg = bioParameters$tEgg,
      tLarva = bioParameters$tLarva,
      tPupa = bioParameters$tPupa,
      AdPopEQ = adultPopEquilibrium,
      inheritanceCube = cube
    )
    
#    netPar$AdPopRatio_F <- matrix(c(1), nrow = 1, dimnames = list(NULL, c("ZW")))
    
#    netPar$AdPopRatio_M <- matrix(c(1-data$introduction_thresholds[[i]], data$introduction_thresholds[[i]]), nrow = 1, dimnames = list(NULL, c("ZZ", "MaMb")))
    
#    netPar$LarPopRatio <- NULL
    
    MGDrivESim <- Network$new(
      params = netPar,
      driveCube = cube,
      patchReleases = patchReleases,
      migrationMale = moveMat,
      migrationFemale = moveMat,
      migrationBatch = basicBatchMigration(batchProbs = 0, sexProbs = c(0.5, 0.5), numPatches = sitesNumber),
      directory = outFolder,
      verbose = FALSE
    )
    
    MGDrivESim$oneRun(verbose = TRUE)
    
    ####################
    # Post-processing
    ####################
    splitOutput(readDir = outFolder, remFile = TRUE, verbose = FALSE)
    aggregateFemales(readDir = outFolder, genotypes = cube$genotypesID, remFile = TRUE, verbose = FALSE)
    
    plot_file <- file.path(outFolder, paste0("plot_threshold_", threshold, "_rM_", rM_value, ".png"))
    png(filename = plot_file, width = 1200, height = 800)
    plotMGDrivESingle(readDir = outFolder, totalPop = TRUE, lwd = 3.5, alpha = 1)
    dev.off()
    
    print(paste("Completed: Threshold =", threshold, "rM =", rM_value))
}

## Simulation end message
print("All simulations complete.")