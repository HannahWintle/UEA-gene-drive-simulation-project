# Load required libraries
library(tidyverse)

# Load datasets
realistic <- read_csv("mgdrive/two_loci/realistic_two_loci_combined.csv") %>%
  mutate(simulation_type = "realistic")

# Filter for rM = 0, generation range, and 3 releases
realistic_data <- realistic %>%
  filter(rM == 0, generation >= 5, generation <= 30, releases == 8)

# Create output directory if it doesn't exist
output_dir <- "plots/two_loci"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define thresholds and custom colour palette
threshold_levels <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
threshold_labels <- c("10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%")
threshold_palette <- c(
  "10%" = "#FAA61A",  # Anchor orange
  "20%" = "#F68B1F",
  "30%" = "#F47C20",
  "40%" = "#E85C29",
  "50%" = "#00AEEF",  # Anchor blue
  "60%" = "#199AD9",
  "70%" = "#3F88C5",
  "80%" = "#5271B3",
  "90%" = "#5E5DAA",
  "100%" = "#6B4EA0"
)

# Convert threshold to factor with percentage labels
realistic_data$threshold <- factor(realistic_data$threshold, levels = threshold_levels, labels = threshold_labels)

# Create plot
plot_rM0_final <- ggplot(realistic_data, aes(x = generation, y = total_females, colour = threshold)) +
  geom_line(alpha = 0.9, size = 1.4) +
  scale_colour_manual(values = threshold_palette) +
  scale_x_continuous(limits = c(5, 30), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, max(realistic_data$total_females), by = 50)) +
  labs(
    x = "Generation",
    y = "Total Wildtype (ZW) Females",
    colour = "Threshold"
  ) +
  theme_classic(base_size = 24) +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 24),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 24),
    plot.title = element_blank()
  )

# Save plot
ggsave(filename = file.path(output_dir, "two_loci_thresholds_all_release_8.png"),
       plot = plot_rM0_final, width = 12, height = 7, dpi = 300)
