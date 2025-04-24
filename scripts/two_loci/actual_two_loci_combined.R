# 01_process_simulation_outputs.R

library(tidyverse)
library(stringr)

# Set your directory where the files are stored
base_dir <- "mgdrive/two_loci/actual_two_loci"

# Grab all female aggregate files
csv_files <- list.files(path = base_dir, pattern = "^F_Aggregate_Runthreshold_.*\\.csv$", recursive = TRUE, full.names = TRUE)

read_female_file <- function(file_path) {
  df <- read_csv(file_path, show_col_types = FALSE)
  
  # Extract full folder path
  folder_path <- dirname(file_path)
  folder_name <- basename(folder_path)
  
  # Extract threshold and rM from the folder name, not the filename
  threshold_str <- str_extract(folder_name, "(?<=threshold_)[0-9\\.eE\\-]+")
  rM_str <- str_extract(folder_name, "(?<=rM_)[0-9\\.eE\\-]+")
  
  # Safely convert to numeric
  threshold <- suppressWarnings(as.numeric(threshold_str))
  rM <- suppressWarnings(as.numeric(rM_str))
  
  df %>%
    mutate(
      generation = row_number(),
      threshold = threshold,
      rM = rM,
      total_females = rowSums(across(any_of(c("ZW", "MaW", "MbW", "RW")), ~ replace_na(., 0)))
    )
}

# Combine all into one dataframe
combined_females <- map_dfr(csv_files, read_female_file)

# Define output directory
current_run <- "mgdrive/two_loci"
dir.create(current_run, recursive = TRUE, showWarnings = FALSE)

# Save combined output
write_csv(combined_females, file.path(current_run, "actual_two_loci_combined.csv"))

# Optional: Check structure
glimpse(combined_females)