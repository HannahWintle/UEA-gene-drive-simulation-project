# Load required libraries
library(tidyverse)

# Load datasets
actual <- read_csv("mgdrive/bisex_runs/actual_bisex_combined.csv") %>%
  mutate(simulation_type = "actual")

realistic <- read_csv("mgdrive/bisex_runs/realistic_bisex_combined.csv") %>%
  mutate(simulation_type = "realistic")

# Combine and filter for rM = 0 and generation 10 to 30
combined_data <- bind_rows(actual, realistic) %>%
  filter(rM == 0, generation >= 10, generation <= 30)

# Create output directory if it doesn't exist
output_dir <- "plots/bisex"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define thresholds and custom colour palette
threshold_levels <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
threshold_palette <- c(
  "#FF6666", "#FF3333", "#CC0000", "#990000",  # Reversed reds for 0.1–0.4
  "#3366CC", "#3399FF", "#66CCFF", "#99CCFF", "#CCE5FF", "#E6F2FF"  # Blues for 0.5–1.0
)

# Convert threshold to factor with correct ordering
combined_data$threshold <- factor(combined_data$threshold, levels = threshold_levels)

# Create plot
plot_rM0_coloured <- ggplot(combined_data, aes(x = generation, y = total_females, colour = threshold)) +
  geom_line(alpha = 0.9, size = 1) +
  scale_colour_manual(values = threshold_palette) +
  facet_wrap(~ simulation_type, ncol = 1) +
  labs(
    title = "Female Population Dynamics by Introduction Threshold (rM = 0)",
    subtitle = "Generations 10–30",
    x = "Generation",
    y = "Total Females",
    colour = "Threshold"
  ) +
  theme_minimal(base_size = 14)

# Save plot
ggsave(filename = file.path(output_dir, "bisex_linegraphs_rM_0_gens10_30_coloured.png"), plot = plot_rM0_coloured, width = 12, height = 8, dpi = 300)