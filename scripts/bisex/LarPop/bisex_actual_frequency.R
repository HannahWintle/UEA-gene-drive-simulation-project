## Testing MEREA in MGDrivE
####################
# Load libraries
####################
library(MGDrivE)

source("cubes/cube_MEREA_with_resistance_allele.R")
source("cubes/cube_auxiliary.R")

current_run <- "mgdrive/bisex_runs/actual_bisex"
dir.create(current_run)

####################
# Define vectors of parameters to vary
####################
introduction_thresholds <- seq(0.1, 1.0, by = 0.1)  # from 10% to 100% introduction frequency
resistance_values <- seq(0, 0.001, by = 0.0001)  # Resistance allele probability #changed

#1.0 = 10%
#0.1 = 1%
#0.01 = 0.1%
#0.001 = 0.01%
#0.0001 = 0.001%

####################
# Simulation Parameters
####################
tMax <- 2000  # Total simulation time
bioParameters <- list(betaK=20, tEgg=5, tLarva=6, tPupa=4, popGrowth=1.175, muAd=0.09)
moveMat <- matrix(data = 1, nrow = 1, ncol = 1)
adultPopEquilibrium <- 500
sitesNumber <- nrow(moveMat)

# Set MGDrivE to deterministic mode
setupMGDrivE(stochasticityON = FALSE, verbose = FALSE)

####################
# for loop
####################
# Create data.frame with all combinations
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
  cube <- cubeMEREA(rM = rM_value)
  
  # Calculate release size
  release_size <- adultPopEquilibrium * threshold
  
  # Release setup
  releasesParameters <- list(
    releasesStart = 0,
    releasesNumber = 3,
    releasesInterval = 30,
    releaseProportion = release_size
  )
  
  maleReleasesVector <- generateReleaseVector(
    driveCube = cube,
    nameGenotypes = list(c("MM", release_size)),
    releasesParameters = releasesParameters
  )
  
  femaleReleasesVector <- generateReleaseVector(
    driveCube = cube,
    nameGenotypes = list(c("MW", release_size)),
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
  
 #Set wildtype larval population
  netPar$AdPopRatio_F <- matrix(c(1), nrow = 1, dimnames = list(NULL, c("ZW")))
  netPar$AdPopRatio_M <- matrix(c(1), nrow = 1, dimnames = list(NULL, c("ZZ")))
  netPar$LarPopRatio <- matrix(c(0, 0), nrow = 1, dimnames = list(NULL, c("ZW", "ZZ")))
  
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
  
  # Display the plot in RStudio Viewer
  plotMGDrivESingle(readDir = outFolder, totalPop = TRUE, lwd = 3.5, alpha = 1)
  
  print(paste("Completed: Threshold =", threshold, "rM =", rM_value))
}

## Simulation end message
print("All simulations complete.")