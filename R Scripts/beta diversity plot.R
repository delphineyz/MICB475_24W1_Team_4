bc_dm <- phyloseq::distance(nasa_rare, method = "bray")
pcoa_bc <- ape::pcoa(bc_dm)
pcoa_scores <- as.data.frame(pcoa_bc$vectors)
colnames(pcoa_scores) <- paste0("Axis.", seq_len(ncol(pcoa_scores)))

# Add metadata
pcoa_scores$humidity_bin <- sample_data(nasa_rare)$humidity_bin

# Define comparisons
comparisons <- list(c("High", "Medium"), c("High", "Low"), c("Medium", "Low"))

# Plot for Axis 1
p1 <- ggplot(pcoa_scores, aes(x = humidity_bin, y = Axis.1, fill = humidity_bin)) +
  geom_boxplot(alpha = 0.3, outlier.shape = 21, outlier.color = "black", outlier.size = 2) +
  scale_fill_manual(values = c("High" = "#FF5733", "Medium" = '#3333FF', "Low" = "#33FF57")) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif") +
  theme_classic(base_size = 15) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "none"
  ) +
  ggtitle("PCoA Axis 1 by Humidity Bin") +
  ylab("PCoA Axis 1")

p1

# Plot for Axis 2
p2 <- ggplot(pcoa_scores, aes(x = humidity_bin, y = Axis.2, fill = humidity_bin)) +
  geom_boxplot(alpha = 0.3, outlier.shape = 21, outlier.color = "black", outlier.size = 2) +
  scale_fill_manual(values = c("High" = "#FF5733", "Medium" = "#3333FF", "Low" = "#33FF57")) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif") +
  theme_classic(base_size = 15) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "none"
  ) +
  ggtitle("PCoA Axis 2 by Humidity Bin") +
  ylab("PCoA Axis 2")

p2

# Combine plots using patchwork or gridExtra
library(patchwork)
p1 + p2 + plot_layout(ncol = 1)  # Combine the two plots vertically

# Combined
pcoa_long <- pivot_longer(pcoa_scores, cols = starts_with("Axis"), names_to = "Axis", values_to = "Score")

# Compute the mean score of Axis 1 and Axis 2
pcoa_scores$mean_axis <- rowMeans(pcoa_scores[, c("Axis.1", "Axis.2")])

# Create the combined box plot
pcoa_box <- ggplot(pcoa_scores, aes(x = humidity_bin, y = mean_axis, fill = humidity_bin)) +
  geom_boxplot(alpha = 0.3, outlier.shape = 21, outlier.color = "black", outlier.size = 2) +
  scale_fill_manual(values = c("High" = "#FF5733", "Medium" = "#3333FF", "Low" = "#33FF57")) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif") +
  theme_classic(base_size = 15) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "none"
  ) +
  ggtitle("Combined PCoA Axes by Humidity Bin") +
  ylab("Mean PCoA Score (Axis 1 & 2)")
pcoa_box



# Extract PCoA scores and eigenvalues (variance explained)
pcoa_scores <- as.data.frame(pcoa_bc$vectors)  # Replace with appropriate function for your PCoA object
colnames(pcoa_scores) <- paste0("Axis.", seq_len(ncol(pcoa_scores)))
eigenvalues <- pcoa_bc$values$Relative_eig  # Relative eigenvalues (proportion of variance explained)

# Ensure eigenvalues match the number of axes
eigenvalues <- eigenvalues[1:ncol(pcoa_scores)]

# Scale the scores by variance explained
for (i in seq_along(eigenvalues)) {
  pcoa_scores[, i] <- pcoa_scores[, i] * eigenvalues[i]
}

# Combine axes into a single weighted score
pcoa_scores$Weighted_Score <- rowSums(pcoa_scores[, 1:2])  # Use the first two axes for simplicity

# Add metadata
pcoa_scores$humidity_bin <- sample_data(nasa_rare)$humidity_bin  # Replace with your metadata object

# Define PERMANOVA p-values (replace with your real results)
permanova_pvalues <- c(
  `High vs Low` = 0.021,
  `High vs Medium` = 0.097,
  `Low vs Medium` = 0.929
)

# Convert p-values to significance labels
permanova_labels <- sapply(permanova_pvalues, function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
})

# Adjusted comparisons data for proper brackets
comparisons <- data.frame(
  Group1 = c("High", "High", "Low"),
  Group2 = c("Low", "Medium", "Medium"),
  Label = c("*", "ns", "ns"),
  y_position = max(pcoa_scores$Weighted_Score) + c(0.2, 0.15, 0.1)  # Adjust Y positions for labels
)

# Plot with full brackets
ggplot(pcoa_scores, aes(x = humidity_bin, y = Weighted_Score, fill = humidity_bin)) +
  geom_boxplot(alpha = 0.3, outlier.shape = 21, outlier.color = "black", outlier.size = 2) +
  scale_fill_manual(values = c("High" = "#FF5733", "Medium" = "#3333FF", "Low" = "#33FF57")) +
  
  # Horizontal lines for brackets
  geom_segment(
    data = comparisons,
    aes(
      x = as.numeric(factor(Group1)), 
      xend = as.numeric(factor(Group2)),
      y = y_position,
      yend = y_position
    ),
    inherit.aes = FALSE, size = 0.5
  ) +
  
  # Vertical lines for brackets
  geom_segment(
    data = comparisons,
    aes(
      x = as.numeric(factor(Group1)), 
      xend = as.numeric(factor(Group1)),
      y = y_position - 0.02,
      yend = y_position
    ),
    inherit.aes = FALSE, size = 0.5
  ) +
  geom_segment(
    data = comparisons,
    aes(
      x = as.numeric(factor(Group2)), 
      xend = as.numeric(factor(Group2)),
      y = y_position - 0.02,
      yend = y_position
    ),
    inherit.aes = FALSE, size = 0.5
  ) +
  
  # Significance labels above brackets
  geom_text(
    data = comparisons,
    aes(
      x = (as.numeric(factor(Group1)) + as.numeric(factor(Group2))) / 2,  # Midpoint of groups
      y = y_position + 0.01,
      label = Label
    ),
    inherit.aes = FALSE, size = 5, vjust = 0
  ) +
  
  theme_classic(base_size = 15) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "none"
  ) +
  ggtitle("Weighted PCoA Axes by Humidity Bin with PERMANOVA") +
  ylab("Weighted PCoA Score")
