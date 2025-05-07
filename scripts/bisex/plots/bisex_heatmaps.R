library(tidyverse)

# Load both datasets
actual <- read_csv("mgdrive/bisex_runs/actual_bisex_combined.csv") %>%
  mutate(simulation_type = "actual")

realistic <- read_csv("mgdrive/bisex_runs/realistic_bisex_combined.csv") %>%
  mutate(simulation_type = "realistic")

# Combine datasets
data_combined <- bind_rows(actual, realistic)

# Filter to final generation only
final_data <- data_combined %>%
  group_by(simulation_type, threshold, rM) %>%
  filter(generation == max(generation)) %>%
  ungroup()

# Convert simulation_type to nicer labels
final_data$simulation_type <- factor(final_data$simulation_type,
                                     levels = c("actual", "realistic"),
                                     labels = c("Controlled Scenario", "Realistic Scenario"))

# Create heatmap plot
heatmap_plot <- ggplot(final_data, aes(x = threshold, y = rM, fill = total_females)) +
  geom_tile() +
  facet_wrap(~ simulation_type) +
  scale_x_continuous(breaks = seq(0.1, 1.0, by = 0.1)) +
  scale_fill_viridis_c(option = "plasma", trans = "log10", name = "Total Females\n(Final Generation)") +
  labs(
    title = "Female Population Size at Final Generation",
    x = "Introduction Threshold",
    y = "Resistance Allele Formation Rate (rM)"
  ) +
  theme_classic(base_size = 14) +
  theme(strip.background = element_blank())

# Display the plot
print(heatmap_plot)

# Save the heatmap
output_dir <- "plots/bisex"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(output_dir, "bisex_heatmap.png"),
  plot = heatmap_plot,
  width = 12,
  height = 8,
  dpi = 300
)
