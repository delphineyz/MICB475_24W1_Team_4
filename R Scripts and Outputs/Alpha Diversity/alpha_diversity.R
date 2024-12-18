# Includes alpha diversity for CRH vs VRH and Humidity Levels
# Run nasa_data_processing.R first
library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(broom)
library(vegan)
library(phyloseq)

# CRH vs. VRH

# Extract sample data from nasa_rare
sample_data_df <- data.frame(sample_data(nasa_rare))

# Compute richness metrics (Observed, Shannon)
richness_metrics <- estimate_richness(nasa_rare, measures = c("Observed", "Shannon"))

# Add treatment and humidity bin information
richness_metrics$treatment <- sample_data_df$treatment
richness_metrics$humidity_bin <- sample_data_df$humidity_bin

# Reshape the data
long_data <- richness_metrics %>%
  pivot_longer(
    cols = c(Observed, Shannon),  # Metrics to include
    names_to = "Metric",                           # Column for metric names
    values_to = "Value"                            # Column for metric values
  )

long_data$treatment <- factor(long_data$treatment,
                              levels = c("Constant Relative Humidity", "Varying Relative Humidity"),
                              labels = c("CRH", "VRH"))

# Create the combined plot
alpha_s1a <- ggplot(long_data, aes(x = treatment, y = Value, fill = treatment)) +
  geom_boxplot(alpha = 0.3, outlier.shape = 21, outlier.color = "black", outlier.size = 2) +
  scale_fill_manual(values = c(
    "CRH" = "#FF5733",
    "VRH" = "#33FF57"
  )) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = TRUE) +
  facet_wrap(~ Metric, scales = "free_y") +
  theme_classic(base_size = 15) +
  theme(
    plot.background = element_rect(fill = "white", color = "black", linewidth = 1),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  labs(
    x = "Treatment",
    y = "Value",
    title = "Alpha Diversity Metrics"
  )

# Print the combined plot
print(alpha_s1a)



# Low, Medium, High

# Define comparisons
comparisons <- list(
  c("Low", "High"),
  c("Medium", "High"),
  c("Low", "Medium")
)

# Reorder the levels of humidity_bin
long_data$humidity_bin <- factor(long_data$humidity_bin, levels = c("Low", "Medium", "High"))

# Create the combined plot
alpha_1a <- ggplot(long_data, aes(x = humidity_bin, y = Value, fill = humidity_bin)) +
  geom_boxplot(alpha = 0.3, outlier.shape = 21, outlier.color = "black", outlier.size = 2) +
  scale_fill_manual(name = "Humidity Bin", values = c("High" = "#FF5733", "Medium" = "#3333FF", "Low" = "#33FF57")) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = FALSE, comparisons = comparisons, bracket.size= 0.3, size = 3) +
  facet_wrap(~ Metric, scales = "free_y") +
  theme_classic(base_size = 15) +
  theme(
    plot.background = element_rect(fill = "white", color = "black", linewidth = 1),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  labs(
    x = "Humidity Bin",
    y = "Value",
    title = "Alpha Diversity Metrics"
  )

# Print the combined plot
print(alpha_1a)
