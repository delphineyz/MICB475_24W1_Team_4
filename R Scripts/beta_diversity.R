#Please run the nasa_data_processing.R script before running this script
library(ggpubr)
library(broom)
library(vegan)

### Figure S1B ----
    # Calculate Bray-Curtis distance and perform PCoA ordination
    bc_dm <- phyloseq::distance(nasa_rare, method = "bray")
    pcoa_bc <- phyloseq::ordinate(nasa_rare, method = "PCoA", distance = bc_dm)
    
    # Plot with enhanced aesthetics and ellipses
    pcoa_s1b <- plot_ordination(nasa_rare, pcoa_bc, color = "treatment") +
      geom_point(size = 4, alpha = 0.7, shape = 21, stroke = 0.5) +  # Custom point size, transparency, and border
      
      # Add transparent ellipses for each treatment group
      stat_ellipse(aes(fill = treatment), geom = "polygon", alpha = 0.2, level = 0.95) +  # 95% confidence level ellipse
      
      # Bright color palette
      scale_color_manual(name = "Treatment", values = c("Constant Relative Humidity" = "#FF5733", "Varying Relative Humidity" = "#33FF57"), labels = c("Constant Relative Humidity" = "CRH", "Varying Relative Humidity" = "VRH")) +
      scale_fill_manual(name = "Treatment",values = c("Constant Relative Humidity" = "#FF5733", "Varying Relative Humidity" = "#33FF57"), labels = c("Constant Relative Humidity" = "CRH", "Varying Relative Humidity" = "VRH")) +
      
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
    pcoa_s1b

    
### Taxa Bar Plot - CRH vs. VRH
    # Convert to relative abundance and group by phylum
    nasa_RA <- transform_sample_counts(nasa_rare, function(x) x/sum(x))
    nasa_phylum <- tax_glom(nasa_RA, taxrank = "Phylum", NArm=FALSE)
    # plot the new plot
    gg_taxa <- plot_bar(nasa_phylum, fill="Phylum") + 
      facet_wrap(.~treatment, scales = "free_x")
    
    # Melt data and calculate the mean relative abundance of each Phylum by treatment
    nasa_phylum_data <- psmelt(nasa_phylum) %>%
      group_by(treatment, Phylum) %>%
      summarize(mean_abundance = mean(Abundance), .groups = "drop")
    
    # Plot mean relative abundance by Phylum within each treatment
    taxa_aim1 <- ggplot(nasa_phylum_data, aes(x = treatment, y = mean_abundance, fill = Phylum)) +
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
    taxa_aim1
    
### Figure 1B ----
# Reorder the levels of humidity_bin
sample_data(nasa_rare)$humidity_bin <- factor(
  sample_data(nasa_rare)$humidity_bin,
  levels = c("Low", "Medium", "High")  # Specify the desired order
)

# Calculate Bray-Curtis distance and perform PCoA ordination
bc_dm <- phyloseq::distance(nasa_rare, method = "bray")
pcoa_bc <- phyloseq::ordinate(nasa_rare, method = "PCoA", distance = bc_dm)

# Plot with enhanced aesthetics and ellipses
pcoa_1b <- plot_ordination(nasa_rare, pcoa_bc, color = "humidity_bin") +
  geom_point(size = 4, alpha = 0.7, shape = 21, stroke = 0.5) +  # Custom point size, transparency, and border
  
  # Add transparent ellipses for each treatment group
  stat_ellipse(aes(fill = humidity_bin), geom = "polygon", alpha = 0.2, level = 0.95) +  # 95% confidence level ellipse
  
  # Bright color palette
  scale_color_manual(name = "Humidity Bin", values = c("High" = "#FF5733", "Medium" = '#3333FF',"Low" = "#33FF57")) +
  scale_fill_manual(name = "Humidity Bin", values = c("High" = "#FF5733", "Medium" = '#3333FF',"Low" = "#33FF57")) +
  
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
pcoa_1b

### Taxa Bar Plot - Low, Medium, High 
# Convert to relative abundance and group by phylum
nasa_RA <- transform_sample_counts(nasa_rare, function(x) x/sum(x))
nasa_phylum <- tax_glom(nasa_RA, taxrank = "Phylum", NArm=FALSE)
# plot the new plot
taxa_humidity <- plot_bar(nasa_phylum, fill="Phylum") + 
  facet_wrap(.~humidity_bin, scales = "free_x")
taxa_humidity

# Melt data and calculate the mean relative abundance of each Phylum by treatment
nasa_phylum_data <- psmelt(nasa_phylum) %>%
  group_by(humidity_bin, Phylum) %>%
  summarize(mean_abundance = mean(Abundance), .groups = "drop")

# Plot mean relative abundance by Phylum within each treatment
taxa_humidity <- ggplot(nasa_phylum_data, aes(x = humidity_bin, y = mean_abundance, fill = Phylum)) +
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
taxa_humidity


    