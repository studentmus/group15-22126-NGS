library(phyloseq)
library(tidyverse)


# ==== Loading data ====

# Load relative abundance data
motus_rel_abundance <- read.table(
  "data/motus_relab_merged.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
)

# Fix sample names
parse_name <- function(s) {
  sub(".*(SRR[0-9]+).*", "\\1", s)
}
motus_rel_abundance <- rename_with(motus_rel_abundance, parse_name)
glimpse(motus_rel_abundance)

# Load sample metadata
samples <- read.csv("data/samples.csv", header = FALSE)
samples <- data.frame(
  sample_id = samples$V1,
  diet = factor(samples$V2, levels = c("vegan", "omnivore"))
)
knitr::kable(samples, "simple")


# ==== Converting data into phyloseq format ====
# See: https://joey711.github.io/phyloseq/import-data

# Function to parse mOTU taxonomy
parse_motus_taxonomy <- function(motus_abundance) {

  # Taxonomy string structure:
  #   d__Domain;p__Phylum;c__Class;o__Order;f__Family;g__Group;s__Species
  # Unassigned taxonomy:
  #   d__;p___A;c__;o__;f__;g__;s__
  motus_names <- motus_abundance$Taxonomy

  # Initialize taxonomy matrix
  tax_mat <- matrix(NA, nrow = length(motus_names), ncol = 7)
  colnames(tax_mat) <- c(
    "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"
  )
  # Use the mOTU4.0_XXXXXX as rownames
  rownames(tax_mat) <- rownames(motus_abundance)

  for (i in seq_along(motus_names)) {
    name <- motus_names[i]

    parts <- strsplit(name, ";")[[1]]
    if (length(parts) != 7) {
      stop("Error: Taxonomy strings should have 7 parts")
    }
    parts_clean <- gsub("^[dpcofgs]__", "", parts)
    # Handle unknown category at each level
    parts_clean[str_starts(parts_clean, "Unknown")] <- NA

    if (parts_clean[1] == "") {
      tax_mat[i, "Domain"] <- "Unclassified"
    } else {
      # Assign to taxonomy levels
      tax_mat[i, ] <- parts_clean
    }
  }

  # Replace empty strings with NA
  tax_mat[tax_mat == ""] <- NA

  return(tax_mat)
}

# Create taxonomy table
TAX <- tax_table(parse_motus_taxonomy(motus_rel_abundance))
head(TAX)

# Create phyloseq OTU table
# Note: taxa_are_rows = TRUE because mOTUs are in rows
otu_mat <- as.matrix(select(motus_rel_abundance, -Taxonomy))
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
head(OTU)

# Create sample metadata
# Rows must match column names in OTU table
metadata <- data.frame(
  sample_id = colnames(otu_mat),
  row.names = colnames(otu_mat)
)
metadata$diet <- samples$diet[match(metadata$sample_id, samples$sample_id)]

# Convert to phyloseq sample_data
SAMP <- sample_data(metadata)
head(SAMP)

# Create phyloseq object
physeq <- phyloseq(OTU, TAX, SAMP)

# Check the phyloseq object
print(physeq)
# Should show:
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ X taxa and 12 samples ]
# sample_data() Sample Data:       [ 12 samples by 2 sample variables ]
# tax_table()   Taxonomy Table:    [ X taxa by 7 taxonomic ranks ]
