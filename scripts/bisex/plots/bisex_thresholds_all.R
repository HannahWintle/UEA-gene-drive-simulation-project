# Load required libraries
library(tidyverse)

# Load datasets
# Load data
data <- read_csv("mgdrive/bisex_runs/dataframe/realistic_bisex_combined.csv")

# Combine and filter for rM = 0 and generation 10 to 30
subset_data <- data %>%
  filter(rM == 0, generation >= 5, generation <= 30, releases == 1)

# Pivot to long format: ZW = wildtype, MW = gene drive, RW = resistant females
allele_data <- subset_data %>%
  select(generation, ZW, MW) %>%
  pivot_longer(cols = c("ZW", "MW"), names_to = "genotype", values_to = "count") %>%
  mutate(genotype = recode(genotype,
                           "ZW" = "Wildtype females",
                           "MW" = "MEREA females"))

# Get y-axis upper limit
y_max <- max(allele_data$count) * 1.05

# Create output directory if it doesn't exist
output_dir <- "plots/bisex"
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
combined_data$threshold <- factor(combined_data$threshold, levels = threshold_levels, labels = threshold_labels)


# Plot
plot_rM0_final <- ggplot(allele_data, aes(x = generation, y = count, colour = threshold)) +
  geom_line(alpha = 0.9, size = 1.4) +
  scale_colour_manual(values = threshold_palette) +
  scale_x_continuous(limits = c(5, 30), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 50, 100, 150, 200, 250)) +
  facet_wrap(~ genotype, ncol = 1) +
  labs(
    x = "Generation",
    y = "Total Females",
    colour = "Threshold"
  ) +
  theme_classic(base_size = 24) +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 24),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 24),
    strip.background = element_blank(),
    strip.text = element_text(size = 24),
    plot.title = element_blank()
  )

# Save plot
ggsave(filename = file.path(output_dir, "bisex_threshold_all.png"),
       plot = plot_rM0_final, width = 12, height = 7, dpi = 300)