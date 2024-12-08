library(dplyr)
library(vegan)

# Fig 1A
# Perform Wilcoxon Rank-Sum Test
wilcox_1a <- combined_metrics %>%
  summarise(
    Observed_p = wilcox.test(Observed ~ treatment, data = ., exact = FALSE)$p.value,
    Chao1_p = wilcox.test(Chao1 ~ treatment, data = ., exact = FALSE)$p.value,
    Shannon_p = wilcox.test(Shannon ~ treatment, data = ., exact = FALSE)$p.value,
    Faith_PD_p = wilcox.test(`Faith's PD` ~ treatment, data = ., exact = FALSE)$p.value
  )

# View results
print(wilcox_1a)

# Fig 1B
# Perform PERMANOVA
permanova_1b <- adonis2(bc_dm ~ treatment, data = sample_data_df, permutations = 999)

# View results
print(permanova_1b)

# Fig 1C
# Melt phyloseq object to get raw abundance data
raw_abundance_data <- psmelt(nasa_phylum)

# Perform Wilcoxon rank-sum test for each Phylum
wilcox_1c <- raw_abundance_data %>%
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
print(wilcox_1c)

# Fig 2a
# Define comparisons
comparisons <- list(
  c("Low", "High"),
  c("Medium", "High"),
  c("Low", "Medium")
)

# Function to compute Wilcoxon p-values for each metric and comparison
compute_p_values <- function(data, comparisons) {
  results <- data %>%
    group_by(Metric) %>%
    summarize(
      p_values = list(
        lapply(comparisons, function(comp) {
          wilcox.test(
            Value ~ humidity_bin,
            data = filter(data, humidity_bin %in% comp),
            exact = FALSE
          )$p.value
        })
      ),
      comparisons = list(comparisons)
    )
  return(results)
}

# Compute p-values
wilcox_2a_results <- compute_p_values(long_data, comparisons)

# Unnest results
wilcox_2a <- wilcox_2a_results %>%
  unnest(cols = c(p_values, comparisons)) %>%
  rowwise() %>%
  mutate(
    Comparison = paste(comparisons[[1]], "vs", comparisons[[2]])
  ) %>%
  ungroup() %>%
  select(Metric, Comparison, p_values)

# Unlist the p_values column
wilcox_2a <- wilcox_2a %>%
  mutate(p_values = unlist(p_values))

# View results
print(wilcox_2a)

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
permanova_2b <- pairwise_permanova(bc_dm, sample_data_df$humidity_bin)

# View results
print(permanova_2b)

# Fig 2c
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
print(kruskal_results)


