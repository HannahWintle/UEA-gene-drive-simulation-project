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


current_run <- "mgdrive/change_threshold_resistance"
dir.create(current_run)


####################
# Define Releases
####################
releases <- seq(0.1, 1, by = 0.1) # release thresholds
resistance_rates <- seq(0, 0.5, by = 0.1) # resistance allele formation (rM)

####################
# Simulation Parameters
####################
tMax <- 1000  # Simulation run time

# Entomological parameters
bioParameters <- list(betaK=20, tEgg=5, tLarva=6, tPupa=4, popGrowth=1.175, muAd=0.09)

# Single-node population model (no migration)
moveMat <- matrix(data = 1, nrow = 1, ncol = 1)

# Population equilibrium size
adultPopEquilibrium <- 500
sitesNumber <- nrow(moveMat)

####################
# Nested loop over resistance and release frequency
####################

for (res in resistance_rates) {
  for (rel in releases) {
    
    # Generate cube with specific resistance rate (rM)
    cube <- cubeMEREA_2L(rM = res)
    
    # Set MGDrivE to deterministic mode
    setupMGDrivE(stochasticityON = FALSE, verbose = FALSE)
    
    # Output folder
    outFolder <- file.path(current_run, paste0("rM_", res, "_release_", rel))
    dir.create(outFolder, recursive = TRUE, showWarnings = FALSE)
    
    # Patch releases setup (empty)
    patchReleases <- replicate(n = sitesNumber,
                               expr = { list(maleReleases = NULL, femaleReleases = NULL,
                                             eggReleases = NULL, matedFemaleReleases = NULL) },
                               simplify = FALSE)
    
    # Define network parameters
    netPar <- parameterizeMGDrivE(
      runID = paste0("rM_", res, "_rel_", rel),
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
    
    # Assign introduction thresholds to starting adult population
    netPar$AdPopRatio_F <- matrix(c(1 - rel, rel), nrow = 1, dimnames = list(NULL, c("ZW", "MaW")))
    netPar$AdPopRatio_M <- matrix(c(1 - rel, rel), nrow = 1, dimnames = list(NULL, c("ZZ", "MaMb")))
    
    # Run simulation
    MGDrivESim <- Network$new(params = netPar,
                              driveCube = cube,
                              patchReleases = patchReleases,
                              migrationMale = moveMat,
                              migrationFemale = moveMat,
                              migrationBatch = basicBatchMigration(batchProbs = 0, sexProbs = c(.5, .5), numPatches = sitesNumber),
                              directory = outFolder,
                              verbose = FALSE)
    
    MGDrivESim$oneRun(verbose = TRUE)
    
    ####################
    # Post-processing and plotting
    ####################
    splitOutput(readDir = outFolder, remFile = TRUE, verbose = FALSE)
    aggregateFemales(readDir = outFolder, genotypes = cube$genotypesID, remFile = TRUE, verbose = FALSE)
    
    # Plot result in RStudio viewer
    plotMGDrivESingle(readDir = outFolder, totalPop = TRUE, lwd = 3.5, alpha = 1)
    
    # Optional: Save plot as PNG
    png(file.path(outFolder, paste0("plot_rM_", res, "_release_", rel, ".png")), width = 1200, height = 800)
    plotMGDrivESingle(readDir = outFolder, totalPop = TRUE, lwd = 3.5, alpha = 1)
    dev.off()
    
    print(paste("Finished rM =", res, "release =", rel))
  }
}