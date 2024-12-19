# Run nasa_data_processing.R and diversity scripts first
library(dplyr)
library(vegan)
library(gt)
library(purrr)
library(tibble)

# Function to format and create table
format_p_values <- function(df, p_col) {
  df %>%
    mutate(
      Significance = case_when(
        !!sym(p_col) < 0.001 ~ "***",
        !!sym(p_col) < 0.01 ~ "**",
        !!sym(p_col) < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
}
create_gt_table <- function(df, title, p_columns, label_list) {
  gt_table <- df %>%
    gt() %>%
    fmt_number(
      columns = all_of(p_columns),
      decimals = 3
    ) %>%
    tab_header(title = title) %>%
    tab_options(
      table.font.names = "Times New Roman",
      table.font.size = 12
    )
    for (col_name in names(label_list)) {
    gt_table <- gt_table %>%
      cols_label(!!col_name := label_list[[col_name]])
  }
  
  return(gt_table)
}

# Fig S1A
# Perform Wilcoxon Rank-Sum Test
wilcox_s1a <- richness_metrics %>%
  summarise(
    Observed_p = wilcox.test(Observed ~ treatment, data = ., exact = FALSE)$p.value,
    Shannon_p = wilcox.test(Shannon ~ treatment, data = ., exact = FALSE)$p.value,
  )

wilcox_s1a <- wilcox_s1a %>%
  format_p_values("Observed_p") %>%
  format_p_values("Shannon_p")

s1a_table <- create_gt_table(
  wilcox_s1a,
  "CRH vs. VRH Alpha Diversity",
  c("Observed_p", "Shannon_p"),
  label_list = list(
    Observed_p = "Observed Diversity P-Value",
    Shannon_p = "Shannon Diversity P-Value",
    Significance = "Significance"
  )
)
s1a_table

# Fig S1B
# Perform PERMANOVA
permanova_s1b <- adonis2(bc_dm ~ treatment, data = sample_data_df, permutations = 999)

# View results
permanova_s1b <- as.data.frame(permanova_s1b) %>%
  format_p_values("Pr(>F)")

s1b_table <- create_gt_table(
  permanova_s1b,
  "CRH vs. VRH Beta Diversity",
  c("Pr(>F)"),
  label_list = list(
    `Pr(>F)` = "P-Value",
    Significance = "Significance"
  )
)
s1b_table

# Fig S1C
# Melt phyloseq object to get raw abundance data
raw_abundance_data <- psmelt(nasa_phylum)

# Perform Wilcoxon rank-sum test for each Phylum
wilcox_s1c <- raw_abundance_data %>%
  group_by(Phylum) %>%
  summarise(
    p_value = wilcox.test(
      Abundance ~ treatment,
      data = cur_data(),
      exact = FALSE
    )$p.value,
    .groups = "drop"
  )

# View results
s1c_table <- wilcox_s1c %>%
  format_p_values("p_value") %>%
  create_gt_table(
    .,
    "Taxa P-Values (CRH vs. VRH)",
    c("p_value"),
    label_list = list(
      Phylum = "Phylum",
      p_value = "P-Value",
      Significance = "Significance"
    )
  )
s1c_table

# Fig 1A
# Define comparisons
comparisons <- list(
  c("Low", "High"),
  c("Medium", "High"),
  c("Low", "Medium")
)

# Function to compute Wilcoxon p-values for each metric and comparison
compute_pairwise_wilcox <- function(data, comparisons) {
  results <- list()
  
  for (metric in unique(data$Metric)) {
    for (comp in comparisons) {
      group1 <- comp[1]
      group2 <- comp[2]
      
      # Subset data for the pairwise comparison
      subset_data <- data %>%
        filter(Metric == metric, humidity_bin %in% c(group1, group2))
      
      # Perform Wilcoxon test
      p_value <- wilcox.test(Value ~ humidity_bin, data = subset_data, exact = FALSE)$p.value
      
      # Save results
      results <- append(results, list(
        data.frame(
          Metric = metric,
          Comparison = paste(group1, "vs", group2),
          p_value = p_value
        )
      ))
    }
  }
  
  # Combine results into a single data frame
  results_df <- do.call(rbind, results)
  
  # Add significance stars
  results_df <- results_df %>%
    mutate(
      Significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  
  return(results_df)
}

# Compute p-values
wilcox_1a <- compute_pairwise_wilcox(long_data, comparisons)

# View results
column_labels <- list(
  Metric = "Diversity Metric",
  Comparison = "Group Comparison",
  p_value = "P-Value",
  Significance = "Significance"
)

# Generate the table
label_list <- list(
  Metric = "Diversity Metric",
  Comparison = "Group Comparison",
  p_value = "P-Value",
  Significance = "Significance"
)

# Generate the table
wilcox_1a_table <- create_gt_table(
  wilcox_1a,
  "Humidity Bin Alpha Diversity",
  p_columns = c("p_value"),
  label_list = label_list
)
wilcox_1a_table

# Fig 2B
# Pairwise PERMANOVA function
pairwise.adonis2 <- function(dist, grouping, permutations = 999) {
  
# Ensure grouping is a factor
  grouping <- factor(grouping)
  groups <- levels(grouping)
  results <- list()
  
# Loop through all unique pairwise comparisons
  for (i in 1:(length(groups) - 1)) {
    for (j in (i + 1):length(groups)) {
      group_i <- groups[i]
      group_j <- groups[j]
      
      # Subset the data
      subset_indices <- grouping %in% c(group_i, group_j)
      subset_dist <- as.dist(as.matrix(dist)[subset_indices, subset_indices])
      subset_grouping <- grouping[subset_indices]
      
      # Ensure the subset has exactly two levels
      if (length(levels(factor(subset_grouping))) < 2) {
        warning(paste("Skipping comparison:", group_i, "vs", group_j, "- insufficient levels"))
        next
      }
      
      # Perform PERMANOVA for the subset
      result <- adonis2(subset_dist ~ subset_grouping, permutations = permutations)
      results[[paste(group_i, "vs", group_j)]] <- result
    }
  }
  
  return(results)
}

# Perform PERMANOVA
pairwise_permanova <- function(distance_matrix, group_vector, permutations = 999) {
  pairwise_results <- pairwise.adonis2(distance_matrix, group_vector, permutations = permutations)
  return(pairwise_results)
}
permanova_1b <- pairwise_permanova(bc_dm, sample_data_df$humidity_bin)

# Extract results
extract_p_values <- function(permanova_list) {
  data.frame(
    Comparison = names(permanova_list),
    P_Value = sapply(permanova_list, function(result) {
      if (!is.null(result)) {
        # Attempt to extract p-value (Pr(>F)) from the Model row
        model_row <- as.data.frame(result)$`Pr(>F)`[1]
        return(model_row)
      } else {
        return(NA) # Handle null or missing data
      }
    }),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      Significance = case_when(
        P_Value < 0.001 ~ "***",
        P_Value < 0.01 ~ "**",
        P_Value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
}
permanova_p_values <- extract_p_values(permanova_1b)

# View Results
labels <- list(
  Comparison = "Group Comparison",
  P_Value = "P-Value",
  Significance = "Significance"
)

permanova_1b_table <- permanova_p_values %>%
  gt() %>%
  fmt_number(columns = "P_Value", decimals = 3) %>%
  cols_label(!!!labels) %>%
  tab_header(title = "Humidity Bin Beta Diversity") %>%
  tab_options(
    table.font.names = "Times New Roman",
    table.font.size = 12
  )
permanova_1b_table

# Fig 1c
# Perform Kruskal-Wallis test for each Phylum
kruskal_results <- raw_abundance_data %>%
  group_by(Phylum) %>%
  summarise(
    p_value = kruskal.test(Abundance ~ humidity_bin, data = cur_data())$p.value
  ) %>%
  mutate(
    Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# View results
kruskal_results <- kruskal_results %>%
  format_p_values("p_value") %>%
  create_gt_table(
    .,
    "Taxa P-Values (Humidity Bin)",
    c("p_value"),
    label_list = list(
      Phylum = "Phylum",
      p_value = "P-Value",
      Significance = "Significance"
    )
  )
kruskal_results
