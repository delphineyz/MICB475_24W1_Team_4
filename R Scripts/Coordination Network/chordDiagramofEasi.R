# Load necessary libraries
library(phyloseq)
library(SpiecEasi)
library(igraph)
library(circlize)
library(dplyr)
library(gridExtra)


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
chord_data$Phylum_from <- sapply(chord_data$Var1, function(x) tax[x, "Phylum"])
chord_data$Phylum_to <- sapply(chord_data$Var2, function(x) tax[x, "Phylum"])

# Check if there are any NA values in the Phylum columns and remove them
chord_data <- chord_data[!is.na(chord_data$Phylum_from) & !is.na(chord_data$Phylum_to), ]

# Step 5: Prepare colors for each Phylum
unique_phyla <- unique(c(chord_data$Phylum_from, chord_data$Phylum_to))
num_phyla <- length(unique_phyla)

# Generate colors using a distinct color palette
colors <- rainbow(num_phyla)
names(colors) <- unique_phyla

# Step 6: Dynamically adjust gaps to match the number of unique sectors
gaps <- rep(2, num_phyla)
gaps[num_phyla] <- 15  # Add a larger gap at the end for better separation

# Step 7: Clear previous circos plot
circos.clear()

# Step 8: Set circos parameters with the correct gap length
circos.par(gap.after = gaps)

# Step 9: Create the chord diagram grouped by Phylum
chordDiagram(
  x = chord_data[, c("Phylum_from", "Phylum_to", "Freq")],
  transparency = 0.6,
  annotationTrack = "grid",
  preAllocateTracks = 1,
  grid.col = colors
)

# Step 10: Add labels for Phylum around the circle with adjusted text
circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    circos.text(
      x = CELL_META$xcenter, y = CELL_META$ylim[1] + 1,
      labels = sector.name, facing = "clockwise", niceFacing = TRUE,
      adj = c(0, 0.5), cex = 0.7, col = "black"
    )
  },
  bg.border = NA
)

# Add a title to the chord diagram
title("Chord Diagram of Microbial Co-occurrence by Phylum")

# Step 11: Add a legend to explain colors for each phylum with a white background
legend(
  "topright", 
  legend = names(colors), 
  fill = colors,
  border = "white", 
  cex = 0.8, 
  bty = "o",    # Box type, set to 'o' to display the box
  box.col = "white",  # Box border color
  bg = "white",       # Background color
  title = "Phylum"
)









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

# Step 5: Prepare colors for each Phylum
unique_genus <- unique(c(chord_data$Genus_from, chord_data$Genus_to))
num_phyla <- length(unique_genus)

# Generate colors using a distinct color palette
colors <- rainbow(num_phyla)
names(colors) <- unique_phyla

# Step 6: Dynamically adjust gaps to match the number of unique sectors
gaps <- rep(2, num_phyla)
gaps[num_phyla] <- 15  # Add a larger gap at the end for better separation

# Step 7: Clear previous circos plot
circos.clear()

# Step 8: Set circos parameters with the correct gap length
circos.par(gap.after = gaps)

# Step 9: Create the chord diagram grouped by Phylum
chordDiagram(
  x = chord_data[, c("Genus_from", "Genus_to", "Freq")],
  transparency = 0.6,
  annotationTrack = "grid",
  preAllocateTracks = 1,
  grid.col = colors
)

# Step 10: Add labels for Phylum around the circle with adjusted text
circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    circos.text(
      x = CELL_META$xcenter, y = CELL_META$ylim[1] + 1,
      labels = sector.name, facing = "clockwise", niceFacing = TRUE,
      adj = c(0, 0.5), cex = 0.7, col = "black"
    )
  },
  bg.border = NA
)

# Add a title to the chord diagram
title("Chord Diagram of Microbial Co-occurrence by Phylum")

# Step 11: Add a legend to explain colors for each phylum with a white background
legend(
  "topright", 
  legend = names(colors), 
  fill = colors,
  border = "white", 
  cex = 0.8, 
  bty = "o",    # Box type, set to 'o' to display the box
  box.col = "white",  # Box border color
  bg = "white",       # Background color
  title = "Phylum"
)