# This is the code to do DESeq and Core Microbiome Analysis
# Please run the nasa_data_processing.R script prior to working on this file
library(microbiome)
library(ggVennDiagram)
library(sf)
library(ggtern)
library(patchwork)

#### DESeq ####

## NOTE: If you get a zeros error, then you need to add '1' count to all reads
nasa_plus1 <- transform_sample_counts(nasa_rare, function(x) x+1)
nasa_deseq <- phyloseq_to_deseq2(nasa_plus1, ~`humidity_bin`)
DESEQ_nasa <- DESeq(nasa_deseq)
res1 <- results(DESEQ_nasa, tidy=TRUE, 
                #this will ensure that No is your reference group
                contrast = c("humidity_bin","Low","Medium"))
View(res1)

res2 <- results(DESEQ_nasa, tidy=TRUE, 
                #this will ensure that No is your reference group
                contrast = c("humidity_bin","Low","High"))
View(res2)



# Add a contrast label to each result set
res1 <- res1 %>% mutate(contrast = "Low vs Medium")
res2 <- res2 %>% mutate(contrast = "Low vs High")

# Combine results and filter for significant ASVs
combined_res <- bind_rows(res1, res2) 

combined_res_sig <- combined_res %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  dplyr::rename(ASV = row)
combined_res

# Extract the significant ASV names from combined_res_sig
sig_ASVs <- combined_res_sig$ASV

# Now filter the original combined_res to include only the significant ASVs
filtered_combined_res <- combined_res %>%
  dplyr::rename(ASV = row) %>%
  filter(ASV %in% sig_ASVs)

# Prune phyloseq object to keep only ASVs of interest and add taxonomy information
nasa_filtered <- prune_taxa(sig_ASVs, nasa_rare)

# Add taxonomy data to the filtered results
sigASVs_final <- tax_table(nasa_filtered) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(filtered_combined_res, by = "ASV")

# Create the bar plot with flipped axes
deseq_bar <- ggplot(sigASVs_final, aes(x = log2FoldChange, y = Genus, fill = contrast)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  scale_y_discrete(labels = function(y) ifelse(y == "g__Chloroplast", "p__Cyanobacteria", y)) +
  geom_errorbar(aes(xmin = log2FoldChange - lfcSE, xmax = log2FoldChange + lfcSE), 
                position = position_dodge(width = 0.7), width = 0.3) +
  scale_fill_manual(values = c("Low vs Medium" = "skyblue", "Low vs High" = "salmon")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1)) +
  ylab("Genus") +
  xlab("log2 Fold Change") +
  theme_minimal() +
  labs(title = "Fold Change Comparison by Contrast (Filtered by Significant Genuses)", fill = "Contrast")
deseq_bar



# Filter combined_res to remove any NA values in padj
volcano_data <- combined_res %>%
  dplyr::rename(ASV = row) %>%
  filter(!is.na(padj))

# Create a new column to flag significant ASVs based on thresholds
volcano_data <- volcano_data %>%
  mutate(
    significant = padj < 0.05 & abs(log2FoldChange) > 1,
    shape = ifelse(contrast == "Low vs Medium", 16, 17)  # Different shapes for each contrast
  )

# Create the volcano plot
volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = significant, shape = contrast)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
  scale_shape_manual(values = c("Low vs Medium" = 16, "Low vs High" = 17)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgreen") +
  theme_minimal() +
  labs(
    title = "Volcano Plot with Contrasts",
    x = "log2 Fold Change",
    y = "-log10(Adjusted p-value)",
    color = "Significant",
    shape = "Contrast"
  ) +
  theme(legend.position = "right")
volcano_plot




















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

# Step 10: Plot bar charts of core ASVs' relative abundance
plot_bar(nasa_core_low, fill = "Genus") +
  ggtitle("Core Microbiome: Low Humidity") +
  theme_minimal()

plot_bar(nasa_core_medium, fill = "Genus") +
  ggtitle("Core Microbiome: Medium Humidity") +
  theme_minimal()

plot_bar(nasa_core_high, fill = "Genus") +
  ggtitle("Core Microbiome: High Humidity") +
  theme_minimal()


# Define a consistent color palette for the genera
genus_colors <- c(
  "g__Corynebacterium" = "#F8766D",
  "g__Pseudomonas" = "#00BA38",
  "g__Staphylococcus" = "#619CFF",
  "g__Streptococcus" = "#C77CFF"
)

# Step 1: Create individual bar plots for each humidity category with consistent colors
low_plot <- plot_bar(nasa_core_low, fill = "Genus") +
  ggtitle("Low Humidity") +
  scale_fill_manual(values = genus_colors) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

medium_plot <- plot_bar(nasa_core_medium, fill = "Genus") +
  ggtitle("Medium Humidity") +
  scale_fill_manual(values = genus_colors) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

high_plot <- plot_bar(nasa_core_high, fill = "Genus") +
  ggtitle("High Humidity") +
  scale_fill_manual(values = genus_colors) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Step 2: Combine plots using patchwork
combined_plot <- low_plot / medium_plot / high_plot +
  plot_layout(ncol = 1) + 
  plot_annotation(
    title = "Core Microbiome Composition by Humidity Bin",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

# Display the combined plot
print(combined_plot)
