library(tidyverse)

# Load data
actual <- read_csv("mgdrive/two_loci/dataframe/actual_two_loci_combined.csv") %>%
  mutate(simulation_type = "Mathematically Ideal Scenario")

realistic <- read_csv("mgdrive/two_loci/dataframe/realistic_two_loci_combined.csv") %>%
  mutate(simulation_type = "More Realistic Scenario")

# Filter for realistic, rM == 0, and final generation only
combined_data <- realistic %>%
  filter(rM == 0, generation == max(generation)) %>%
  mutate(
    ZW = replace_na(ZW, 0),
    percent_ZW = pmin((ZW / 250) * 100, 100),
    threshold = factor(threshold),
    releases = factor(releases)
  )

# Plot
heatmap_plot <- ggplot(combined_data, aes(x = threshold, y = releases, fill = percent_ZW)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(
    option = "plasma",
    name = "% ZW Females",
    limits = c(0, 100),
    na.value = "grey90"
  ) +
  labs(
    title = "ZW Female Population at Final Generation (Realistic Scenario, rM = 0)",
    x = "Introduction Threshold",
    y = "Number of Releases",
    fill = "% ZW Females"
  ) +
  theme_light(base_size = 24) +
  theme(
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 24),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 24),
    plot.title = element_text(size = 26, face = "bold")
  )

# Save the plot
ggsave("plots/two_loci/two_loci_heatmap.png",
       plot = heatmap_plot,
       width = 12, height = 8, dpi = 300)
