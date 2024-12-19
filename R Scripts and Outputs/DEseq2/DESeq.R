# Please run the nasa_data_processing.R script prior to working on this file
library(DESeq2)
library(phyloseq)
library(tidyverse)
library(ape)
library(sf)
library(ggtern)
library(patchwork)
library(dplyr)
library(tibble)

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
