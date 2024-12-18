### Lets try building some co-occurrence networks
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
nasa_low_humidity <- subset_samples(nasa_rare, humidity_bin == "Low")
nasa_med_humidity <- subset_samples(nasa_rare, humidity_bin == "Medium")
nasa_high_humidity <- subset_samples(nasa_rare, humidity_bin == "High")








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
mouse_filtered <- filter_by_genus(nasa_med_humidity)



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
  #legend(
   # "topleft",
   # legend = names(phylum_colors),
   # col = phylum_colors,
   # pch = 19,
   # cex = 0.8,
    #title = "Phylum"
  #)
  
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





