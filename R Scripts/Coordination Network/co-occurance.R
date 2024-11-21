### Lets try building some co-occurance networks
library(SpiecEasi)
library(igraph)
library(phyloseq)
library(ggplot2)
library(tidygraph)
library(ggraph)
library(viridis)
library(RColorBrewer)
library(dplyr)




# Step 1: Filter to include only samples with 'humidity_category' equal to 'Low'
mouse_low_humidity <- subset_samples(mouse_rare, humidity_category == "Low")
mouse_med_humidity <- subset_samples(mouse_rare, humidity_category == "Medium")
mouse_high_humidity <- subset_samples(mouse_rare, humidity_category == "High")








# Step 1: Filter by genus with a minimum abundance threshold
filter_by_genus <- function(physeq_obj, min_count = 1, min_samples = 1) {
  # Aggregate to genus level
  physeq_genus <- tax_glom(physeq_obj, taxrank = "Genus")
  
  # Remove low-abundance genera
  physeq_genus <- prune_taxa(taxa_sums(physeq_genus) > min_count, physeq_genus)
  
  # Keep only genera that appear in at least 'min_samples' samples
  physeq_genus <- filter_taxa(physeq_genus, function(x) sum(x > 0) >= min_samples, prune = TRUE)
  
  return(physeq_genus)
}

# Step 2: Filter the dataset at the genus level using the function above
mouse_filtered <- filter_by_genus(mouse_high_humidity)



genera_to_label <- c(
  "g_Veillonella", "g_UCG-002", "g_Streptococcus", "g_Prevotella", 
  "g_Peptoniphilus", "g_Murdochiella", "g_Fusobacterium", "g_Fastidiosipila", 
  "g_Corynebacterium", "g_Chloroplast", "g_Brevibacterium", "g_Anaerococcus", "g_Finegoldia"
)


generate_network_igraph <- function(physeq_obj, seed = 42, genera_to_label) {
  # Step 1: Apply data transformations and generate SpiecEasi network
  print("Applying data transformations...")
  otu_data <- as.matrix(otu_table(physeq_obj))
  se.mb <- spiec.easi(otu_data, method = 'mb', pulsar.params = list(rep.num = 40))
  print("Selecting model with pulsar using stars...")
  print("Fitting final estimate with mb... done")
  print("SpiecEasi analysis completed")
  
  # Step 2: Convert to igraph object
  adj_matrix <- getRefit(se.mb)
  ig_network <- adj2igraph(adj_matrix, vertex.attr = list(name = taxa_names(physeq_obj)))
  print("igraph Network Created")
  
  # Step 3: Add taxonomic information to vertices
  taxa_info <- as.data.frame(tax_table(physeq_obj))
  V(ig_network)$Genus <- taxa_info$Genus[match(V(ig_network)$name, rownames(taxa_info))]
  V(ig_network)$Phylum <- taxa_info$Phylum[match(V(ig_network)$name, rownames(taxa_info))]
  
  # Step 4: Check abundance data
  otu_abundance <- taxa_sums(physeq_obj)
  valid_names <- intersect(V(ig_network)$name, names(otu_abundance))
  otu_abundance <- otu_abundance[valid_names]
  
  # Step 5: Assign sizes to vertices
  V(ig_network)$size <- sapply(V(ig_network)$name, function(x) {
    value <- otu_abundance[x]
    if (is.list(value)) value <- unlist(value)[1]
    if (!is.null(value) && !is.na(value) && is.numeric(value)) {
      return(log10(value + 1) * 4)
    } else {
      return(1)
    }
  })
  
  # Ensure sizes are numeric
  if (any(sapply(V(ig_network)$size, is.list))) {
    V(ig_network)$size <- sapply(V(ig_network)$size, function(x) if (is.list(x)) unlist(x)[1] else x)
  }
  V(ig_network)$size <- as.numeric(V(ig_network)$size)
  
  # Step 6: Generate colors for the phyla
  unique_phyla <- unique(V(ig_network)$Phylum)
  phylum_colors <- rainbow(length(unique_phyla))
  names(phylum_colors) <- unique_phyla
  V(ig_network)$color <- phylum_colors[V(ig_network)$Phylum]
  
  # Step 7: Assign label and border colors for genera of interest and their neighbors
  V(ig_network)$label_color <- "black"
  V(ig_network)$border_color <- "black"
  
  # Step 8: Label nodes of interest in red and their neighbors in black
  V(ig_network)$label <- NA
  
  # Label the specified genera in red
  V(ig_network)$label[V(ig_network)$Genus %in% genera_to_label] <- V(ig_network)$Genus[V(ig_network)$Genus %in% genera_to_label]
  V(ig_network)$label_color[V(ig_network)$Genus %in% genera_to_label] <- "red"
  
  # Find nodes directly connected to the specified genera
  connected_to_labeled <- unique(unlist(sapply(genera_to_label, function(g) {
    if (g %in% V(ig_network)$Genus) {
      neighbor_vertices <- neighbors(ig_network, V(ig_network)$Genus == g)
      return(V(ig_network)[neighbor_vertices]$name)
    }
  })))
  
  # Label the neighbors in black
  V(ig_network)$label[V(ig_network)$name %in% connected_to_labeled] <- V(ig_network)$Genus[V(ig_network)$name %in% connected_to_labeled]
  V(ig_network)$label_color[V(ig_network)$name %in% connected_to_labeled] <- "black"
  
  # Step 9: Ensure degrees are numeric
  node_degrees <- degree(ig_network)
  if (any(sapply(node_degrees, is.list))) {
    node_degrees <- sapply(node_degrees, function(x) if (is.list(x)) unlist(x)[1] else x)
  }
  node_degrees <- as.numeric(node_degrees)
  
  # Remove nodes with zero degree (unconnected)
  ig_network <- delete_vertices(ig_network, V(ig_network)[node_degrees == 0])
  
  # Step 10: Plot the network
  set.seed(seed)
  layout <- layout_with_fr(ig_network)
  
  plot(
    ig_network,
    layout = layout,
    vertex.label = V(ig_network)$label,
    vertex.label.color = V(ig_network)$label_color,
    vertex.label.cex = 0.6,
    vertex.size = V(ig_network)$size,
    vertex.color = V(ig_network)$color,
    vertex.frame.color = V(ig_network)$border_color,
    edge.width = 0.5,
    edge.color = "gray50",
    main = "Microbial Network by Phylum"
  )
  
  # Add legend
  legend(
    "topleft",
    legend = names(phylum_colors),
    col = phylum_colors,
    pch = 19,
    cex = 0.8,
    title = "Phylum"
  )
  
  print("Plotting completed")
}
# Run the function
generate_network_igraph(mouse_filtered, genera_to_label = genera_to_label)






























# Define genera of interest
genera_to_label <- c(
  "g_Veillonella", "g_UCG-002", "g_Streptococcus", "g_Prevotella", 
  "g_Peptoniphilus", "g_Murdochiella", "g_Fusobacterium", "g_Fastidiosipila", 
  "g_Corynebacterium", "g_Chloroplast", "g_Brevibacterium", "g_Anaerococcus", "g_Finegoldia"
)

# Generate network with igraph
generate_network_igraph <- function(physeq_obj, seed = 42, genera_to_label) {
  print("Applying data transformations...")
  
  # Extract OTU table
  otu_data <- as.matrix(otu_table(physeq_obj))
  
  # Apply SpiecEasi with adjusted parameters
  se.mb <- spiec.easi(
    otu_data, 
    method = 'mb', 
    pulsar.params = list(rep.num = 20),
    nlambda = 50,
    lambda.min.ratio = 1e-3
  )
  print("SpiecEasi analysis completed")
  
  # Extract adjacency matrix
  adj_matrix <- getRefit(se.mb)
  
  # Check if the adjacency matrix is empty
  if (sum(adj_matrix) == 0) {
    stop("Adjacency matrix is empty. Adjust SpiecEasi parameters.")
  }
  
  # Convert to igraph
  ig_network <- adj2igraph(adj_matrix, vertex.attr = list(name = taxa_names(physeq_obj)))
  
  # Add taxonomic information
  taxa_info <- as.data.frame(tax_table(physeq_obj))
  V(ig_network)$Genus <- taxa_info$Genus[match(V(ig_network)$name, rownames(taxa_info))]
  V(ig_network)$Phylum <- taxa_info$Phylum[match(V(ig_network)$name, rownames(taxa_info))]
  
  # Handle missing values in taxonomic information
  V(ig_network)$Genus[is.na(V(ig_network)$Genus)] <- "Unknown"
  V(ig_network)$Phylum[is.na(V(ig_network)$Phylum)] <- "Unknown"
  
  # Assign sizes to vertices based on abundance
  otu_abundance <- taxa_sums(physeq_obj)
  valid_names <- intersect(V(ig_network)$name, names(otu_abundance))
  otu_abundance <- otu_abundance[valid_names]
  V(ig_network)$size <- log10(otu_abundance + 1) * 4
  
  # Generate colors for phyla
  unique_phyla <- unique(V(ig_network)$Phylum)
  phylum_colors <- brewer.pal(min(length(unique_phyla), 8), "Set1")
  names(phylum_colors) <- unique_phyla
  V(ig_network)$color <- phylum_colors[V(ig_network)$Phylum]
  
  # Assign labels to specified genera
  V(ig_network)$label <- NA
  V(ig_network)$label_color <- "black"
  V(ig_network)$label[V(ig_network)$Genus %in% genera_to_label] <- V(ig_network)$Genus[V(ig_network)$Genus %in% genera_to_label]
  V(ig_network)$label_color[V(ig_network)$Genus %in% genera_to_label] <- "red"
  
  # Remove unconnected nodes
  ig_network <- delete_vertices(ig_network, V(ig_network)[degree(ig_network) == 0])
  
  # Plot the network
  set.seed(seed)
  layout <- layout_with_fr(ig_network)
  
  plot(
    ig_network,
    layout = layout,
    vertex.label = V(ig_network)$label,
    vertex.label.color = V(ig_network)$label_color,
    vertex.label.cex = 0.6,
    vertex.size = V(ig_network)$size,
    vertex.color = V(ig_network)$color,
    vertex.frame.color = "gray30",
    edge.width = 0.5,
    edge.color = "gray50",
    main = "Microbial Network by Phylum"
  )
  
  # Add legend for phyla colors
  legend(
    "topleft",
    legend = names(phylum_colors),
    col = phylum_colors,
    pch = 19,
    cex = 0.8,
    bty = "n",
    title = "Phylum",
    bg = "white"
  )
  
  print("Plotting completed")
}

# Run the function
generate_network_igraph(mouse_filtered, genera_to_label = genera_to_label)






































































# Assign sizes to vertices with extra checks for list issues
V(ig_network)$size <- sapply(V(ig_network)$name, function(x) {
  value <- otu_abundance[x]
  
  # Handle cases where the value is a list
  if (is.list(value)) {
    print(paste("List detected for:", x, value))
    value <- unlist(value)[1]  # Convert to numeric if it's a list
  }
  
  # Convert to numeric if still not numeric after unlisting
  value <- as.numeric(value)
  
  # Assign size based on abundance, ensuring it's a numeric value
  if (!is.null(value) && !is.na(value) && is.finite(value)) {
    return(log10(value + 1) * 4)
  } else {
    return(1)  # Default size if abundance is missing or invalid
  }
})

# Ensure sizes are numeric
if (!is.numeric(V(ig_network)$size)) {
  stop("Node sizes are not numeric")
}

# Check if any sizes are lists (final safeguard)
if (any(sapply(V(ig_network)$size, is.list))) {
  V(ig_network)$size <- sapply(V(ig_network)$size, function(x) if (is.list(x)) unlist(x)[1] else x)
}

# Debug: Print the first few node sizes to confirm
print("Node sizes after assignment:")
print(head(V(ig_network)$size))











































































# Function to filter by genus with a minimum abundance threshold
filter_by_genus <- function(physeq_obj, min_count = 10, min_samples = 3) {
  # Aggregate to genus level
  physeq_genus <- tax_glom(physeq_obj, taxrank = "Genus")
  
  # Remove low-abundance genera
  physeq_genus <- prune_taxa(taxa_sums(physeq_genus) > min_count, physeq_genus)
  
  # Keep only genera that appear in at least 'min_samples' samples
  physeq_genus <- filter_taxa(physeq_genus, function(x) sum(x > 0) >= min_samples, prune = TRUE)
  
  return(physeq_genus)
}

# Subset by humidity category and filter for meaningful genera
mouse_low <- subset_samples(mouse_rare, humidity_category == "Low")
mouse_medium <- subset_samples(mouse_rare, humidity_category == "Medium")
mouse_high <- subset_samples(mouse_rare, humidity_category == "High")

# Filter each subset for genus level analysis
mouse_low_filtered <- filter_by_genus(mouse_low)
mouse_medium_filtered <- filter_by_genus(mouse_medium)
mouse_high_filtered <- filter_by_genus(mouse_high)




generate_network_ggraph <- function(physeq_obj) {
  # Run SpiecEasi analysis
  se.mb <- spiec.easi(physeq_obj, method = 'mb', lambda.min.ratio = 1e-2, nlambda = 20, pulsar.params = list(rep.num = 2))
  
  # Convert to igraph object
  ig_network <- adj2igraph(getRefit(se.mb), vertex.attr = list(name = taxa_names(physeq_obj)))
  
  # Extract the taxonomy table
  tax <- tax_table(physeq_obj)
  
  # Ensure taxonomy table is valid and has the correct dimensions
  if (is.null(tax) || ncol(tax) < 1) {
    stop("Taxonomy table is not properly set.")
  }
  
  # Extract Genus information while handling missing data
  OTU_ids <- V(ig_network)$name
  genus_info <- tax[OTU_ids, "Genus"]
  
  # Check for NAs and handle missing genera
  if (any(is.na(genus_info))) {
    genus_info[is.na(genus_info)] <- "Unknown"
  }
  
  # Assign genus names to the nodes
  V(ig_network)$Genus <- as.character(genus_info)
  
  # Calculate node sizes based on OTU abundance
  otu_abundance <- taxa_sums(physeq_obj)
  V(ig_network)$size <- log10(otu_abundance[OTU_ids] + 1) * 3
  
  # Convert to tidy graph object
  network_tidy <- as_tbl_graph(ig_network)
  
  # Generate the network plot using ggraph
  ggraph(network_tidy, layout = "fr") +
    geom_edge_link(aes(edge_alpha = 0.5), color = "gray70") +
    geom_node_point(aes(color = Genus, size = size), alpha = 0.8) +
    scale_color_viridis_d(na.translate = FALSE) +  # Automatically assign colors
    theme_minimal() +
    labs(title = "Microbial Network", color = "Genus") +
    theme(legend.position = "right")
}

# Run the analysis for the filtered low humidity dataset
generate_network_ggraph(mouse_low_filtered)



# Function to filter genera by abundance
filter_by_abundance <- function(physeq_obj, min_abundance = 50) {
  # Aggregate at the genus level
  physeq_genus <- tax_glom(physeq_obj, taxrank = "Genus")
  
  # Filter genera based on the minimum abundance
  physeq_genus <- prune_taxa(taxa_sums(physeq_genus) > min_abundance, physeq_genus)
  return(physeq_genus)
}

# Apply the function to your filtered low humidity dataset
filtered_mouse_low <- filter_by_abundance(mouse_low_filtered, min_abundance = 100)

# Generate a simplified network plot using your filtered data
generate_network_ggraph(filtered_mouse_low)















# Subset phyloseq object for Low humidity category
mouse_rare_low <- subset_samples(mouse_rare, humidity_category == "Low")
mouse_rare_low <- prune_taxa(taxa_sums(mouse_rare_low) > 0, mouse_rare_low)

# Subset phyloseq object for Medium humidity category
mouse_rare_medium <- subset_samples(mouse_rare, humidity_category == "Medium")
mouse_rare_medium <- prune_taxa(taxa_sums(mouse_rare_medium) > 0, mouse_rare_medium)

# Subset phyloseq object for High humidity category
mouse_rare_high <- subset_samples(mouse_rare, humidity_category == "High")
mouse_rare_high <- prune_taxa(taxa_sums(mouse_rare_high) > 0, mouse_rare_high)



# Function to generate and plot the network
generate_network <- function(physeq_obj, color_palette) {
  # Run SpiecEasi analysis
  se.mb <- spiec.easi(physeq_obj, method = 'mb', lambda.min.ratio = 1e-2, nlambda = 20, pulsar.params = list(rep.num = 5))
  
  # Convert to igraph object
  ig_network <- adj2igraph(getRefit(se.mb), vertex.attr = list(name = taxa_names(physeq_obj)))
  
  # Add taxonomy information to the vertices
  V(ig_network)$Phylum <- tax_table(physeq_obj)[V(ig_network)$name, "Phylum"]
  
  # Calculate node sizes based on OTU abundance
  otu_abundance <- taxa_sums(physeq_obj)
  V(ig_network)$size <- log10(otu_abundance[V(ig_network)$name] + 1) * 3
  
  # Assign colors to nodes based on their Phylum
  V(ig_network)$color <- color_palette[V(ig_network)$Phylum]
  
  # Plot the network
  plot(ig_network, 
       vertex.label = NA,
       vertex.size = V(ig_network)$size,
       vertex.color = V(ig_network)$color,
       edge.width = 0.5,
       edge.color = "gray50",
       main = "Microbial Network")
  
  # Add legend
  legend("topright", legend = names(color_palette), col = color_palette, pch = 19, bty = "n", cex = 0.8)
}

# Define a color palette for different phyla
phyla_colors <- c(
  "Chloroflexi" = "#A3A500", "Acidobacteria" = "#E76BF3", 
  "Bacteroidetes" = "#E7B800", "Verrucomicrobia" = "#FC4E07",
  "Actinobacteria" = "#00BFC4", "Fusobacteria" = "#F8766D",
  "Planctomycetes" = "#C77CFF", "Firmicutes" = "#00BA38",
  "Proteobacteria" = "#619CFF"
)


# Generate network for Low humidity category
generate_network(mouse_rare_low, phyla_colors)














# Function to generate network graphs for a specific humidity category
generate_network <- function(physeq_obj, humidity_level, color_palette) {
  # Ensure the input object is a phyloseq object
  if (!inherits(physeq_obj, "phyloseq")) {
    stop("Input object is not a phyloseq object")
  }
  
  # Ensure the sample data contains the 'humidity_category' column
  if (!"humidity_category" %in% colnames(sample_data(physeq_obj))) {
    stop("The sample data does not contain 'humidity_category'")
  }
  
  # Convert the humidity_category column to a factor
  sample_data(physeq_obj)$humidity_category <- factor(sample_data(physeq_obj)$humidity_category)
  
  # Subset the phyloseq object based on the humidity level
  subset_physeq <- subset_samples(physeq_obj, humidity_category == humidity_level)
  
  # Check if there are samples left after subsetting
  if (nsamples(subset_physeq) == 0) {
    stop(paste("No samples found for humidity category:", humidity_level))
  }
  
  # Filter out taxa with zero counts after subsetting
  subset_physeq <- prune_taxa(taxa_sums(subset_physeq) > 0, subset_physeq)
  
  # Run SpiecEasi analysis
  se.mb <- spiec.easi(subset_physeq, method = 'mb', lambda.min.ratio = 1e-2, nlambda = 20, pulsar.params = list(rep.num = 50))
  
  # Convert to igraph object
  ig_network <- adj2igraph(getRefit(se.mb), vertex.attr = list(name = taxa_names(subset_physeq)))
  
  # Add taxonomy information to the vertices
  V(ig_network)$Phylum <- tax_table(subset_physeq)[V(ig_network)$name, "Phylum"]
  
  # Calculate node sizes based on OTU abundance
  otu_abundance <- taxa_sums(subset_physeq)
  V(ig_network)$size <- log10(otu_abundance[V(ig_network)$name] + 1) * 3
  
  # Assign colors to nodes based on their Phylum
  V(ig_network)$color <- color_palette[V(ig_network)$Phylum]
  
  # Plot the network
  plot(ig_network, 
       vertex.label = NA,
       vertex.size = V(ig_network)$size,
       vertex.color = V(ig_network)$color,
       edge.width = 0.5,
       edge.color = "gray50",
       main = paste("Microbial Network - Humidity Category:", humidity_level))
  
  # Add legend
  legend("topright", legend = names(color_palette), col = color_palette, pch = 19, bty = "n", cex = 0.8)
}

# Define a color palette for different phyla
phyla_colors <- c(
  "Chloroflexi" = "#A3A500", "Acidobacteria" = "#E76BF3", 
  "Bacteroidetes" = "#E7B800", "Verrucomicrobia" = "#FC4E07",
  "Actinobacteria" = "#00BFC4", "Fusobacteria" = "#F8766D",
  "Planctomycetes" = "#C77CFF", "Firmicutes" = "#00BA38",
  "Proteobacteria" = "#619CFF"
)

# Test the function with the "Low" humidity category
generate_network(mouse_rare, "Low", color_palette = phyla_colors)






















# Step 1: Subset the phyloseq object for the "Low" humidity category
subset_low <- subset_samples(mouse_rare, humidity_category == "Low")
# Run SpiecEasi on the "Low" humidity subset
set.seed(123)
spiec_out_low <- spiec.easi(subset_low, method = 'mb', lambda.min.ratio = 1e-2,
                            nlambda = 20, pulsar.params = list(rep.num = 20, ncores = 4))

# Extract adjacency matrix and convert it to an igraph object
adj_mat_low <- as.matrix(getRefit(spiec_out_low))
net_low <- graph_from_adjacency_matrix(adj_mat_low, mode = "undirected", diag = FALSE)

# Add metadata to the nodes
V(net_low)$name <- taxa_names(subset_low)
V(net_low)$Phylum <- as.character(tax_table(subset_low)[, "Phylum"])

















# Function to create and plot a network for a given condition
create_network <- function(physeq, condition) {
  # Subset data based on humidity category
  subset_physeq <- subset_samples(physeq, humidity_category == condition)
  
  # Run SpiecEasi
  set.seed(123)
  spiec_out <- spiec.easi(subset_physeq, method = 'mb', lambda.min.ratio = 1e-2,
                          nlambda = 20, pulsar.params = list(rep.num = 20, ncores = 4))
  
  # Extract adjacency matrix and convert to igraph object
  adj_mat <- as.matrix(getRefit(spiec_out))
  net <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
  V(net)$name <- taxa_names(subset_physeq)
  V(net)$Phylum <- tax_table(subset_physeq)[, "Phylum"]
  
  # Plot network using ggraph
  ggraph(as_tbl_graph(net), layout = "fr") + 
    geom_edge_link(color = "grey70") + 
    geom_node_point(aes(color = Phylum), size = 5) + 
    theme_minimal() + 
    labs(title = paste("Microbial Network for", condition, "Humidity"), color = "Phylum") +
    theme(legend.position = "right")
}


# Create and plot networks for each condition
plot_low <- create_network(mouse_rare, "Low")
plot_medium <- create_network(mouse_rare, "Medium")
plot_high <- create_network(mouse_rare, "High")

# Display the plots
plot_low
plot_medium
plot_high












