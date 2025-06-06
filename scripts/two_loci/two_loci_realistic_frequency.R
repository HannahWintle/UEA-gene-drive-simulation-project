library(MGDrivE)

source("cubes/cube_MEREA_two_loci.R")
source("cubes/cube_auxiliary.R")

current_run <- "mgdrive/two_loci/realistic_two_loci_resistance"
dir.create(current_run, recursive = TRUE)

# Parameters to vary
introduction_thresholds <- seq(0.5, 1.5, by = 0.1)
resistance_values <- c(0, 0.0001)
releasesNumber_values <- seq(6, 16, by = 1)

# Simulation constants
tMax <- 2000
bioParameters <- list(betaK = 20, tEgg = 5, tLarva = 6, tPupa = 4, popGrowth = 1.175, muAd = 0.09)
moveMat <- matrix(data = 1, nrow = 1, ncol = 1)
adultPopEquilibrium <- 500
sitesNumber <- nrow(moveMat)

setupMGDrivE(stochasticityON = FALSE, verbose = FALSE)

# Build parameter grid
param_grid <- expand.grid(
  introduction_thresholds = introduction_thresholds,
  resistance_values = resistance_values,
  releasesNumber = releasesNumber_values
)

for (i in 1:nrow(param_grid)) {
  threshold <- param_grid$introduction_thresholds[i]
  rM_value <- param_grid$resistance_values[i]
  n_releases <- param_grid$releasesNumber[i]
  
  outFolder <- paste0(current_run, "/threshold_", threshold, "_rM_", rM_value, "_nRel_", n_releases)
  if (dir.exists(outFolder)) {
    message("✅ Already completed: ", outFolder)
    next
  }
  
  print(paste("Running: Threshold =", threshold, "rM =", rM_value, "Releases =", n_releases))
  
  cube <- cubeMEREA_2L(rM = rM_value)
  release_size <- adultPopEquilibrium * threshold
  
  releasesParameters <- list(
    releasesStart = 100,
    releasesNumber = n_releases,
    releasesInterval = 30,
    releaseProportion = release_size
  )
  
  maleReleasesVector <- generateReleaseVector(
    driveCube = cube,
    releasesParameters = releasesParameters,
    nameGenotypes = list(c("MaMb", release_size))
  )
  
  femaleReleasesVector <- generateReleaseVector(
    driveCube = cube,
    releasesParameters = releasesParameters,
    nameGenotypes = list(
      c("MaW", release_size / 2),
      c("MbW", release_size / 2)                   
    )
  )
  
  patchReleases <- replicate(n = sitesNumber,
                             expr = list(maleReleases = NULL, femaleReleases = NULL,
                                         eggReleases = NULL, matedFemaleReleases = NULL),
                             simplify = FALSE)
  patchReleases[[1]]$maleReleases <- maleReleasesVector
  patchReleases[[1]]$femaleReleases <- femaleReleasesVector
  
  dir.create(outFolder, recursive = TRUE, showWarnings = FALSE)
  
  netPar <- parameterizeMGDrivE(
    runID = paste0("threshold_", threshold, "_rM_", rM_value, "_nRel_", n_releases),
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
  
  netPar$AdPopRatio_F <- matrix(c(1), nrow = 1, dimnames = list(NULL, c("ZW")))
  netPar$AdPopRatio_M <- matrix(c(1), nrow = 1, dimnames = list(NULL, c("ZZ")))
  netPar$LarPopRatio <- matrix(c(0.5, 0.5), nrow = 1, dimnames = list(NULL, c("ZW", "ZZ")))
  
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
  
  splitOutput(readDir = outFolder, remFile = TRUE, verbose = FALSE)
  aggregateFemales(readDir = outFolder, genotypes = cube$genotypesID, remFile = TRUE, verbose = FALSE)
  
  png(filename = file.path(outFolder, paste0("plot_threshold_", threshold, "_rM_", rM_value, "_nRel_", n_releases, ".png")),
      width = 400, height = 300)
  par(mar = c(4, 4, 2, 1))
  plotMGDrivESingle(readDir = outFolder, totalPop = TRUE, lwd = 2, alpha = 1)
  dev.off()
  
  print(paste("Completed: Threshold =", threshold, "rM =", rM_value, "Releases =", n_releases))
}

print("All simulations complete.")