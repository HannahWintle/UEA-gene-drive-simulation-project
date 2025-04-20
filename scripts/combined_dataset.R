# Load necessary libraries
library(tidyverse)

# Set directories for simulations
sim_dirs <- c("mgdrive/two_loci/actual_two_loci",
              "mgdrive/two_loci/test_runs/realistic_two_loci",
              "mgdrive/bisex_runs/actual_bisex",
              "mgdrive/bisex_runs/realistic_bisex")

# Create a directory for combined results
current_run <- "mgdrive/combined"
dir.create(current_run, recursive = TRUE, showWarnings = FALSE)

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
      
      # Compute total females across genotypes
      df <- df %>%
        mutate(total_females = rowSums(across(-Generation, -threshold, -resistance, -simulation, -run_folder)))
      
      return(df)
    } else {
      warning(paste("File does not exist:", female_file))
      return(NULL)
    }
  })
  
  return(sim_data)
}

# Combine datasets from all simulations
all_simulations_df <- map_df(sim_dirs, read_simulation)

# Inspect combined dataframe
print(glimpse(all_simulations_df))

# Save the combined dataframe explicitly to your created folder
write_csv(all_simulations_df, file.path(current_run, "combined_simulation_results.csv"))