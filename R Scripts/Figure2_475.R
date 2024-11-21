# This is the code to do AIM 2, it includes code to do alpha and beta diversity metrics
#Please run the nasa_data_processing.R script before running this script

### Figure 1.a: Alpha Diversity comparison of metadata treatment column ----

# Step 1: Calculate Faith's Phylogenetic Diversity (PD)
phylo_dist <- pd(t(otu_table(nasa_rare)), phy_tree(nasa_rare), include.root = FALSE)

# Step 2: Add Faith's PD to the sample data in the phyloseq object
sample_data(nasa_rare)$Faith_PD <- phylo_dist$PD

# Step 3: Generate the plot for standard alpha diversity indices
gg_richness <- plot_richness(
  nasa_rare, 
  x = "humidity_bin", 
  measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
) +
  geom_boxplot(aes(fill = humidity_bin), alpha = 0.3, outlier.shape = 21, outlier.color = "black", outlier.size = 2) +
  scale_fill_manual(values = c("High" = "#FF5733", "Medium" = '#3333FF',"Low" = "#33FF57")) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = TRUE) +
  theme_classic(base_size = 15) +
  theme(
    plot.background = element_rect(fill = "white", color = "black", linewidth = 1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  ggtitle("Alpha Diversity of Constant vs. Varying Relative Humidity") +
  ylab("Alpha Diversity Index")

# Step 4: Generate the separate plot for Faith's PD with adjusted dimensions
plot_faith <- ggplot(sample_data(nasa_rare), aes(x = humidity_bin, y = Faith_PD, fill = humidity_bin)) +
  geom_boxplot(alpha = 0.3, outlier.shape = 21, outlier.color = "black", outlier.size = 2) +
  scale_fill_manual(values = c("High" = "#FF5733", "Medium" = '#3333FF',"Low" = "#33FF57")) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = TRUE) +
  theme_classic(base_size = 15) +
  theme(
    plot.background = element_rect(fill = "white", color = "black", linewidth = 1),
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.ticks.x = element_blank(), # Remove x-axis ticks
    axis.title.x = element_blank(), # Remove x-axis title
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.title.y = element_blank(), # Remove y-axis title
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  ggtitle("") # Remove the plot title

# Add a custom label for Faith's PD
faith_label <- ggplot() + 
  annotate("text", x = 0.5, y = 0.5, label = "Faith's PD", size = 5, fontface = "bold") + 
  theme_void()

# Plot 1: Alpha Diversity measures
print(gg_richness)
# Plot 2: Faith's Phylogenetic Diversity (PD)
print(plot_faith)



### Figure 1.b ----
# Calculate Bray-Curtis distance and perform PCoA ordination
bc_dm <- phyloseq::distance(nasa_rare, method = "bray")
pcoa_bc <- phyloseq::ordinate(nasa_rare, method = "PCoA", distance = bc_dm)

# Plot with enhanced aesthetics and ellipses
gg_pcoa <- plot_ordination(nasa_rare, pcoa_bc, color = "humidity_bin") +
  geom_point(size = 4, alpha = 0.7, shape = 21, stroke = 0.5) +  # Custom point size, transparency, and border
  
  # Add transparent ellipses for each treatment group
  stat_ellipse(aes(fill = humidity_bin), geom = "polygon", alpha = 0.2, level = 0.95) +  # 95% confidence level ellipse
  
  # Bright color palette
  scale_color_manual(values = c("High" = "#FF5733", "Medium" = '#3333FF',"Low" = "#33FF57")) +
  scale_fill_manual(values = c("High" = "#FF5733", "Medium" = '#3333FF',"Low" = "#33FF57")) +
  
  # Black border and clean background
  theme_classic(base_size = 15) +
  theme(
    plot.background = element_rect(fill = "white", color = "black", size = 1),  # Black border
    panel.grid.major = element_line(color = "grey90", size = 0.2),  # Light grid lines for readability
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  ) +
  
  # Add titles and axis labels
  ggtitle("Bray-Curtis PCoA of Beta Diversity by Treatment") +
  xlab("PCoA Axis 1") + 
  ylab("PCoA Axis 2")

# Display the enhanced PCoA plot with ellipses
gg_pcoa

### Figure 1.c ---- 
# Convert to relative abundance and group by phylum
nasa_RA <- transform_sample_counts(nasa_rare, function(x) x/sum(x))
nasa_phylum <- tax_glom(nasa_RA, taxrank = "Phylum", NArm=FALSE)
# plot the new plot
gg_taxa <- plot_bar(nasa_phylum, fill="Phylum") + 
  facet_wrap(.~humidity_bin, scales = "free_x")
gg_taxa

# Melt data and calculate the mean relative abundance of each Phylum by treatment
nasa_phylum_data <- psmelt(nasa_phylum) %>%
  group_by(humidity_bin, Phylum) %>%
  summarize(mean_abundance = mean(Abundance), .groups = "drop")

# Plot mean relative abundance by Phylum within each treatment
gg_taxa <- ggplot(nasa_phylum_data, aes(x = humidity_bin, y = mean_abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent_format()) +  # Display as percentages
  theme_classic(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.background = element_rect(fill = "white", color = "black", size = 1)
  ) +
  ylab("Mean Relative Abundance (%)") +
  ggtitle("Mean Relative Abundance of Phylum by Treatment")

# Display the plot
gg_taxa
