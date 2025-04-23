# Load required libraries
library(tidyverse)

# Load datasets
actual <- read_csv("mgdrive/bisex_runs/actual_bisex_combined.csv") %>%
  mutate(simulation_type = "actual")

realistic <- read_csv("mgdrive/bisex_runs/realistic_bisex_combined.csv") %>%
  mutate(simulation_type = "realistic")

# Combine and filter for rM = 0 and generation 10 to 30
combined_data <- bind_rows(actual, realistic) %>%
  filter(rM == 0.001, generation >= 10, generation <= 40)

combined_data$simulation_type <- factor(combined_data$simulation_type,
                                        levels = c("actual", "realistic"),
                                        labels = c("Controlled Scenario", "Realistic Scenario"))

# Create output directory if it doesn't exist
output_dir <- "plots/bisex"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define thresholds and custom colour palette
threshold_levels <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
threshold_labels <- c("10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%")
threshold_palette <- c(
  "#FF6666", "#FF3333", "#CC0000", "#990000",  # Light to deep red (10%–40%)
  "#3366CC", "#3399FF", "#66CCFF", "#99CCFF", "#CCE5FF", "#E6F2FF"  # Blues (50%–100%)
)

# Convert threshold to factor with percentage labels
combined_data$threshold <- factor(combined_data$threshold, levels = threshold_levels, labels = threshold_labels)

# Create plot
plot_rM0_final <- ggplot(combined_data, aes(x = generation, y = total_females, colour = threshold)) +
  geom_line(alpha = 0.9, size = 1) +
  geom_hline(data = combined_data %>% filter(simulation_type == "Controlled Scenario" & generation == 10),
             aes(yintercept = 0), colour = "black") +
  scale_colour_manual(values = threshold_palette) +
  scale_x_continuous(limits = c(10, 40), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 50, 100, 150, 200, 250)) +
  facet_wrap(~ simulation_type, ncol = 1) +
  labs(
    title = "Female Population Dynamics by Introduction Threshold (High Resistance)",
    x = "Generation",
    y = "Total Females",
    colour = "Threshold"
  ) +
  theme_classic(base_size = 14) +
  theme(strip.background = element_blank())

# Save plot
ggsave(filename = file.path(output_dir, "bisex_linegraphs_rM_high.png"),
       plot = plot_rM0_final, width = 12, height = 8, dpi = 300)