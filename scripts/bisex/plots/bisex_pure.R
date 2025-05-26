# Load libraries
library(tidyverse)

# Output directory
output_dir <- "plots/bisex"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load dataset
data <- read_csv("mgdrive/bisex_runs/dataframe/actual_bisex_combined.csv")

# Filter for 1 release, rM = 0 (no resistance), all thresholds, and a generation range
filtered_data <- data %>%
  filter(releases == 1, rM == 0, generation >= 1, generation <= 30)

# Map numeric thresholds to labels
threshold_levels <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
threshold_labels <- c("10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%")
filtered_data$threshold <- factor(filtered_data$threshold, levels = threshold_levels, labels = threshold_labels)

# Colour palette
threshold_palette <- c(
  "10%" = "#FAA61A",
  "20%" = "#F68B1F",
  "30%" = "#F47C20",
  "40%" = "#E85C29",
  "50%" = "#00AEEF",
  "60%" = "#199AD9",
  "70%" = "#3F88C5",
  "80%" = "#5271B3",
  "90%" = "#5E5DAA",
  "100%" = "#6B4EA0"
)

# Plot
plot_single_release <- ggplot(filtered_data, aes(x = generation, y = total_females, colour = threshold)) +
  geom_line(size = 1.2) +
  scale_colour_manual(values = threshold_palette) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Generation",
    y = "Total Wildtype Females",
    colour = "Threshold"
  ) +
  theme_classic(base_size = 22) +
  theme(
    axis.title = element_text(size = 26, face = "bold"),
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 20)
  )

# Save plot
ggsave(filename = file.path(output_dir, "single_release_all_thresholds.png"),
       plot = plot_single_release, width = 12, height = 8, dpi = 300)
