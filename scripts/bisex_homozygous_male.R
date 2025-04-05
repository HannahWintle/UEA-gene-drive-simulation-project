## Testing MEREA in MGDrivE
####################
# Load libraries
####################
library(MGDrivE)

source("cubes/cube_MEREA_with_resistance_allele.R")
source("cubes/cube_auxiliary.R")


current_run <- "mgdrive/bisex_homozygous_male_test"
dir.create(current_run)


####################
# Define Releases
####################
releases <- seq(50, 500, 50)  # Generates: 50, 100, 150, ..., 500

####################
# Output Folder Setup
####################
dir.create("mgdrive", showWarnings = FALSE)

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
# Basic Inheritance pattern
####################
# Using the MEREA cube instead of MEDEA
cube <- cubeMEREA()

####################
# Loop through each release scenario
####################
for (i in 1:length(releases)) {
  
  # Define release start at day 100
  releases_list <- list(releasesStart=100, releasesNumber=3, releasesInterval=30, releaseProportion=releases[[i]])
  
  releasesParametersMale <- releases_list
  releasesParametersFemale <- releases_list
  
  # Initialize release settings
  patchReleases <- replicate(n=sitesNumber,
                             expr={list(maleReleases=NULL, femaleReleases=NULL,
                                        eggReleases=NULL, matedFemaleReleases=NULL)},
                             simplify=FALSE)
  
  # Generate release vectors (original version)
  maleReleasesVector <- generateReleaseVector(driveCube=cube, 
                                              nameGenotypes = list(c("MM", releases[[i]])),
                                              releasesParameters=releasesParametersMale)
  femaleReleasesVector <- generateReleaseVector(driveCube=cube, 
                                                nameGenotypes = list(c("MW", releases[[i]])),
                                                releasesParameters=releasesParametersFemale)
  
  # Assign to patches
  patchReleases[[1]]$maleReleases <- maleReleasesVector
  patchReleases[[1]]$femaleReleases <- femaleReleasesVector
  
  # Define unique output folder for each run
  outFolder <- paste0(current_run,"/", releases[[i]])
  dir.create(outFolder, recursive = TRUE, showWarnings = FALSE)
  
  ####################
  # Combine parameters and run simulation
  ####################
  
  # Set MGDrivE to deterministic mode
  setupMGDrivE(stochasticityON = FALSE, verbose = FALSE)
  
  # Define network parameters
  netPar <- parameterizeMGDrivE(runID=i, simTime=tMax, sampTime=30, nPatch=sitesNumber,
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
                            migrationBatch=basicBatchMigration(batchProbs=0, sexProbs=c(.5,.5), numPatches=sitesNumber),
                            directory=outFolder,
                            verbose=FALSE)
  
  # Run simulation for this release scenario
  MGDrivESim$oneRun(verbose = TRUE)
  
  ####################
  # Post Analysis
  ####################
  # Split output by patch (for plotting later)
  splitOutput(readDir = outFolder, remFile = TRUE, verbose = FALSE)
  
  # Aggregate female genotypes to match male columns
  aggregateFemales(readDir = outFolder, genotypes = cube$genotypesID,
                   remFile = TRUE, verbose = FALSE)
  
  # Define file path for saving plot
  plot_file <- file.path(outFolder, paste0("plot_release_", releases[[i]], ".png"))
  
  # Save the plot
  png(filename = plot_file, width = 1200, height = 800)  # Set resolution and size
  plotMGDrivESingle(readDir = outFolder, totalPop = TRUE, lwd = 3.5, alpha = 1)
  dev.off()  # Close the PNG device
  
  print(paste("Simulation and plot saved for release size:", releases[[i]], "at", plot_file))
}
