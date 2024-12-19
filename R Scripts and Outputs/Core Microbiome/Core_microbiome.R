



#### Core Microbiome----

# Step 2: Convert to relative abundance
nasa_RA <- transform_sample_counts(nasa_rare, function(x) x / sum(x))

# Step 3: Subset the phyloseq object by humidity_bin categories
nasa_low <- subset_samples(nasa_RA, humidity_bin == "Low")
nasa_medium <- subset_samples(nasa_RA, humidity_bin == "Medium")
nasa_high <- subset_samples(nasa_RA, humidity_bin == "High")

# Prune samples with zero counts to avoid errors
nasa_low <- prune_samples(sample_sums(nasa_low) > 0, nasa_low)
nasa_medium <- prune_samples(sample_sums(nasa_medium) > 0, nasa_medium)
nasa_high <- prune_samples(sample_sums(nasa_high) > 0, nasa_high)

# Step 4: Define detection and prevalence thresholds
detection_threshold <- 0.001  # Minimum relative abundance
prevalence_threshold <- 0.5   # Core taxa must be present in at least 50% of samples

# Step 5: Identify core ASVs for each subset using core_members()
core_low_ASVs <- core_members(nasa_low, detection = detection_threshold, prevalence = prevalence_threshold)
core_medium_ASVs <- core_members(nasa_medium, detection = detection_threshold, prevalence = prevalence_threshold)
core_high_ASVs <- core_members(nasa_high, detection = detection_threshold, prevalence = prevalence_threshold)

# Step 6: Create lists of core ASVs for Venn diagram
core_list <- list(
  "Low Humidity" = core_low_ASVs,
  "Medium Humidity" = core_medium_ASVs,
  "High Humidity" = core_high_ASVs
)

# Step 7: Generate the Venn diagram
venn_plot <- ggVennDiagram(core_list, label = "count") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "Venn Diagram of Core Microbiome by Humidity Bin")

# Save the Venn diagram as a PNG file
ggsave("core_microbiome_venn.png", venn_plot, width = 8, height = 6)

# Step 8: Display the Venn diagram
print(venn_plot)























# Step 9: Prune phyloseq objects to only keep core ASVs for each category
nasa_core_low <- prune_taxa(core_low_ASVs, nasa_low)
nasa_core_medium <- prune_taxa(core_medium_ASVs, nasa_medium)
nasa_core_high <- prune_taxa(core_high_ASVs, nasa_high)


# Define a consistent color palette for the genera
genus_colors <- c(
  "g__Corynebacterium" = "#F8766D",
  "g__Pseudomonas" = "#00BA38",
  "g__Staphylococcus" = "#619CFF",
  "g__Streptococcus" = "#C77CFF"
)



# Extract taxonomy table and convert to data frame
taxonomy_df <- as.data.frame(tax_table(nasa_RA)) %>%
  rownames_to_column(var = "ASV")

# Create a combined data frame for the core microbiome across humidity bins
core_microbiome_data <- data.frame(
  ASV = unique(c(
    rownames(tax_table(nasa_core_low)),
    rownames(tax_table(nasa_core_medium)),
    rownames(tax_table(nasa_core_high))
  )),
  Low = 0,
  Medium = 0,
  High = 0
)

# Fill in relative abundances for each humidity bin
core_microbiome_data$Low <- ifelse(
  core_microbiome_data$ASV %in% rownames(tax_table(nasa_core_low)),
  taxa_sums(nasa_core_low)[core_microbiome_data$ASV] / sum(taxa_sums(nasa_core_low)),
  0
)

core_microbiome_data$Medium <- ifelse(
  core_microbiome_data$ASV %in% rownames(tax_table(nasa_core_medium)),
  taxa_sums(nasa_core_medium)[core_microbiome_data$ASV] / sum(taxa_sums(nasa_core_medium)),
  0
)

core_microbiome_data$High <- ifelse(
  core_microbiome_data$ASV %in% rownames(tax_table(nasa_core_high)),
  taxa_sums(nasa_core_high)[core_microbiome_data$ASV] / sum(taxa_sums(nasa_core_high)),
  0
)

# Normalize relative abundances
core_microbiome_data <- core_microbiome_data %>%
  mutate(
    Total = Low + Medium + High,
    Low = Low / Total,
    Medium = Medium / Total,
    High = High / Total
  ) %>%
  select(-Total)

# Map ASV to genera
core_microbiome_data <- core_microbiome_data %>%
  left_join(taxonomy_df %>% select(ASV, Genus), by = "ASV") %>%
  mutate(Genus = ifelse(is.na(Genus), "Unclassified", Genus))

# Generate the ternary plot
ternary_plot <- ggtern(data = core_microbiome_data, aes(x = Low, y = Medium, z = High, color = Genus)) +
  geom_point(size = 4, alpha = 0.8) +             # Add points for each genus
  geom_text(aes(label = Genus), hjust = -0.3, vjust = -0.3, size = 3, check_overlap = TRUE) +  # Add genus labels
  scale_color_manual(values = genus_colors, na.value = "gray") +  # Use consistent colors for genera
  theme_classic(base_size = 15) +                 # Classic theme for clarity
  labs(
    title = "Core Microbiome Composition (Ternary Plot)",
    x = "Low Humidity",
    y = "Medium Humidity",
    z = "High Humidity"
  ) +
  theme(
    legend.position = "right",                     # Display legend for clarity
    plot.title = element_text(size = 16, face = "bold")
  )

# Display the plot
print(ternary_plot)

# Generate the ternary plot
ternary_plot <- ggtern(data = core_microbiome_data, aes(x = Low, y = Medium, z = High, color = Genus)) +
  geom_point(size = 4, alpha = 0.8) +             # Add points for each genus
  geom_text(aes(label = Genus), hjust = -0.2, vjust = -0.2, size = 3, check_overlap = TRUE) +  # Add genus labels
  scale_color_manual(values = genus_colors, na.value = "gray") +  # Use consistent colors for genera
  theme_bw(base_size = 15) +                      # Use a clean black and white theme
  labs(
    title = "Core Microbiome Composition (Ternary Plot)",
    x = "Low Humidity",
    y = "Medium Humidity",
    z = "High Humidity"
  ) +
  theme(
    legend.position = "right",                     # Display legend for clarity
    plot.title = element_text(size = 16, face = "bold"),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)  # Add space around the plot
  )

# Display the corrected ternary plot
print(ternary_plot)

                                   










                                   
