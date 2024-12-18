### load the necessary packages ----
library(phyloseq)
library(tidyverse)
library(ape)
library(vegan)
library(picante)
library(dplyr)
library(DESeq2)

### load the data ----
meta <- read_delim("475_files/updated_nasa_dust_metadata.tsv", delim="\t")
otu <- read_delim("475_files/feature-table.txt", delim="\t", skip=1)
tax <- read_delim("475_files/taxonomy.tsv", delim="\t")
phylotree <- read.tree("475_files/tree.nwk")


# Format OTU table ----
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 

# Format sample metadata ----
samp_df <- as.data.frame(meta[,-1])
rownames(samp_df)<- meta$'sample-id'
SAMP <- sample_data(samp_df)

#### Formatting taxonomy ----
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix()
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$`Feature ID`
TAX <- tax_table(tax_mat)

#### Make the phyloseq object ----
phylo_nasa <- phyloseq(OTU, SAMP, TAX, phylotree)


#### Processing and Rarefying ----
phylo_nasa_filt <- subset_taxa(phylo_nasa,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
phylo_nasa_filt_nolow <- filter_taxa(phylo_nasa_filt, function(x) sum(x)>5, prune = TRUE)
phylo_nasa_filt_nolow_samps <- prune_samples(sample_sums(phylo_nasa_filt_nolow)>100, phylo_nasa_filt_nolow)
phylo_nasa_filt_nolow_samps_sub <- subset_samples(phylo_nasa_filt_nolow_samps, !is.na(humidity_bin) )


rarecurve(t(as.data.frame(otu_table(phylo_nasa_filt_nolow_samps_sub))), cex=0.1)
nasa_rare <- rarefy_even_depth(phylo_nasa_filt_nolow_samps_sub, rngseed = 1, sample.size = 709)



