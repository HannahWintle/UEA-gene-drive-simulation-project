# Load necessary libraries
library(tidyverse)

# Set directories for simulations
sim_dirs <- c("mgdrive/two_loci/actual_two_loci",
              "mgdrive/two_loci/test_runs/realistic_two_loci_006")

# Function to read and process individual simulation outputs
read_simulation <- function(sim_path) {
  
  # List all subdirectories corresponding to parameter combinations
  sub_dirs <- list.dirs(path = sim_path, recursive = FALSE, full.names = TRUE)
  
  # Process each subdirectory
  sim_data <- map_df(sub_dirs, function(sub_dir) {
    
    # Extract parameter values from folder names using regex
    params <- str_match(basename(sub_dir), "threshold_(\\d\\.\\d+)_rM_(\\d\\.\\d+)") %>%
      as.data.frame() %>%
      rename(threshold = V2, resistance = V3) %>%
      mutate(threshold = as.numeric(threshold),
             resistance = as.numeric(resistance))
    
    # Read aggregated female data CSV produced by aggregateFemales()
    female_file <- file.path(sub_dir, "AggFemale.csv")
    
    if (file.exists(female_file)) {
      
      df <- read_csv(female_file, col_types = cols()) %>%
        mutate(threshold = params$threshold,
               resistance = params$resistance,
               simulation = basename(sim_path),
               run_folder = basename(sub_dir))
      
      # Compute total females across genotypes (assuming genotype columns start from the 2nd column onward)
      df <- df %>%
        mutate(total_females = rowSums(across(-Generation, -threshold, -resistance, -simulation, -run_folder)))
      
      return(df)
    } else {
      # Handle cases where CSV does not exist or is missing
      warning(paste("File does not exist:", female_file))
      return(NULL)
    }
  })
  
  return(sim_data)
}

# Apply the function across all simulation directories
all_simulations_df <- map_df(sim_dirs, read_simulation)

# Inspect combined dataframe
glimpse(all_simulations_df)

# Save the combined dataframe for easy access later
write_csv(all_simulations_df, "combined_simulation_results.csv")
