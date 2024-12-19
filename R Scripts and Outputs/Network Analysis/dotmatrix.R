# Step 1: Generate the SpiecEasi network
set.seed(42)

otu_data <- as.matrix(otu_table(mouse_filtered))
se.mb <- spiec.easi(otu_data, method = 'mb', pulsar.params = list(rep.num = 20))

# Step 2: Extract the adjacency matrix from the SpiecEasi object
adj_matrix <- as.matrix(getRefit(se.mb))

# Step 3: Convert the adjacency matrix to a long format
chord_data <- as.data.frame(as.table(adj_matrix))
chord_data <- chord_data[chord_data$Freq > 0 & chord_data$Var1 != chord_data$Var2, ]

# Step 4: Add taxonomy information for better grouping
tax <- tax_table(mouse_filtered)


# Ensure that the Phylum column is properly extracted
chord_data$Genus_from <- sapply(chord_data$Var1, function(x) tax[x, "Genus"])
chord_data$Genus_to <- sapply(chord_data$Var2, function(x) tax[x, "Genus"])

# Check if there are any NA values in the Phylum columns and remove them
chord_data <- chord_data[!is.na(chord_data$Genus_from) & !is.na(chord_data$Genus_to), ]















# Prepare the data for the triangular dot matrix plot
dot_data <- chord_data %>%
  filter(Genus_from != Genus_to) %>%
  mutate(Pair = paste(pmin(Genus_from, Genus_to), pmax(Genus_from, Genus_to), sep = "_")) %>%
  distinct(Pair, .keep_all = TRUE)



# Define the genera of interest
highlighted_genera <- c(
  "g_Veillonella", "g_UCG-002", "g_Streptococcus", "g_Prevotella", 
  "g_Peptoniphilus", "g_Murdochiella", "g_Fusobacterium", "g_Fastidiosipila", 
  "g_Corynebacterium", "g_Chloroplast", "g_Brevibacterium", "g_Anaerococcus", "g_Finegoldia"
)

# Create the plot with bolded text for specific genera on the axes
ggplot(dot_data, aes(x = Genus_from, y = Genus_to, size = Freq)) +
  geom_point(color = "red") +  # Set all dots to the same color
  scale_size_continuous(range = c(4, 10), guide = "none") +  # Remove frequency legend
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1, size = 12,
      face = ifelse(dot_data$Genus_from %in% highlighted_genera, "bold", "plain")
    ),
    axis.text.y = element_text(
      size = 12,
      face = ifelse(dot_data$Genus_to %in% highlighted_genera, "bold", "plain")
    ),
    axis.ticks = element_line(),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  coord_fixed(ratio = 0.5, clip = "off") +  # Adjust coordinate ratio
  ggtitle("Triangular Dot Matrix Plot of Microbial Co-occurrences by Phylum")
