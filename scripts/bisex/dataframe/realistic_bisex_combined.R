# 01_process_simulation_outputs.R

library(tidyverse)
library(stringr)

# Set your directory where the files are stored
base_dir <- "mgdrive/bisex_runs/realistic_bisex_runs"

# Grab all female aggregate files
csv_files <- list.files(path = base_dir, pattern = "^F_Aggregate_Runthreshold_.*\\.csv$", recursive = TRUE, full.names = TRUE)

read_female_file <- function(file_path) {
  df <- read_csv(file_path, show_col_types = FALSE)
  
  # Extract full folder path
  folder_path <- dirname(file_path)
  folder_name <- basename(folder_path)
  
  # Extract threshold, rM, and number of releases from the folder name
  threshold_str <- str_extract(folder_name, "(?<=threshold_)[0-9\\.eE\\-]+")
  rM_str <- str_extract(folder_name, "(?<=rM_)[0-9\\.eE\\-]+")
  nRel_str <- str_extract(folder_name, "(?<=nRel_)[0-9]+")
  
  # Safely convert to numeric
  threshold <- suppressWarnings(as.numeric(threshold_str))
  rM <- suppressWarnings(as.numeric(rM_str))
  releases <- suppressWarnings(as.numeric(nRel_str))
  
  df %>%
    mutate(
      generation = row_number(),
      threshold = threshold,
      rM = rM,
      releases = releases,
      total_females = rowSums(across(any_of(c("ZW", "MW", "RW")), ~ replace_na(., 0)))
    )
}

# Combine all into one dataframe
combined_females <- map_dfr(csv_files, read_female_file)

# Define output directory
current_run <- "mgdrive/bisex_runs"
dir.create(current_run, recursive = TRUE, showWarnings = FALSE)

# Save combined output
write_csv(combined_females, file.path(current_run, "realistic_bisex_combined.csv"))

# Optional: Check structure
glimpse(combined_females)