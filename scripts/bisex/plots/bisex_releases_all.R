# Load required libraries
library(tidyverse)

# Load only the realistic dataset
realistic <- read_csv("mgdrive/bisex_runs/dataframe/realistic_bisex_combined.csv") %>%
  mutate(simulation_type = "realistic")

# Filter for relevant parameters
filtered_data <- realistic %>%
  filter(rM == 0, generation >= 1, generation <= 25, threshold == 0.5, releases %in% 1:6)

# Output directory
output_dir <- "plots/bisex"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define release levels and colours
releases_levels <- c(1, 2, 3, 4, 5, 6)
releases_labels <- c("1", "2", "3", "4", "5", "6")
releases_palette <- c(
  "1" = "#F47C20",
  "2" = "#E85C29",
  "3" = "#00AEEF",
  "4" = "#199AD9",
  "5" = "#3F88C5",
  "6" = "#5271B3"
)

filtered_data$releases <- factor(filtered_data$releases, levels = releases_levels, labels = releases_labels)

# Plot without facet
plot_realistic_only <- ggplot(filtered_data, aes(x = generation, y = ZW, colour = releases)) +
  geom_line(size = 1.4) +
  scale_colour_manual(values = releases_palette) +
  scale_x_continuous(limits = c(1, 25), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, max(filtered_data$ZW), by = 50)) +
  labs(
    x = "Generation",
    y = "Wildtype (ZW) Females",
    colour = "Number of Releases"
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
ggsave(filename = file.path(output_dir, "realistic_all_releases.png"),
       plot = plot_realistic_only, width = 12, height = 7, dpi = 300)