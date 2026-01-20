library(phyloseq)
library(ggplot2)
library(tidyverse)


figures_dir <- "data/figures/"

# ==== Loading data ====

# Load relative abundance data
motus_rel_abundance <- read.table(
  "data/motus_relab_merged.tsv",
  header=TRUE,
  row.names=1,
  sep="\t",
)
# Remove taxonomy info
# motus_rel_abundance <- select(motus_rel_abundance, -Taxonomy)

# Fix sample names
parse_name <- function(s) {
  sub(".*(SRR[0-9]+).*", "\\1", s)
}
motus_rel_abundance <- rename_with(motus_rel_abundance, parse_name)
glimpse(motus_rel_abundance)
dim(motus_rel_abundance)

# Load sample metadata
samples <- read.csv("data/samples.csv", header = FALSE)
samples <- data.frame(
  sample_id = samples$V1,
  diet = factor(samples$V2, levels = c("vegan", "omnivore"))
)
knitr::kable(samples, "simple")


# ==== Converting data into phyloseq format ====

# Function to parse mOTU taxonomy
parse_motus_taxonomy <- function(motus_abundance) {
  
  motus_names = motus_abundance$Taxonomy

  # Initialize taxonomy matrix
  tax_mat <- matrix(NA, nrow = length(motus_names), ncol = 7)
  colnames(tax_mat) <- c(
    "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
  )
  rownames(tax_mat) <- rownames(motus_abundance)
  
  for(i in seq_along(motus_names)) {
    name <- motus_names[i]
    
    parts <- strsplit(name, ";")[[1]]
    if (length(parts) != 7) {
      stop("Error: Taxonomy strings should have 7 parts")
    }
    parts_clean <- gsub("^[dpcofgs]__", "", parts)
    
    if (parts_clean[1] == "") {
      tax_mat[i, "Kingdom"] <- "Unclassified"
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
rownames(OTU)

# Create sample metadata
# Rows must match column names in OTU table
metadata <- data.frame(
  sample_id = colnames(otu_mat),
  row.names = colnames(otu_mat)
)
# merge(metadata, samples, by = "sample_id")
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


# ==== Barplots of sample composition at phylum-level ====

# Aggregate to phylum level
physeq_phylum <- tax_glom(physeq, taxrank = "Phylum")

# Transform to relative abundance - skip (already relative)
# physeq_phylum_rel <- transform_sample_counts(physeq_phylum, function(x) x / sum(x))

# Basic bar plot
plot_bar(physeq_phylum, fill = "Phylum") +
  theme_bw() +
  labs(
    title = "Phylum-level composition",
    y = "Relative Abundance"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig_path <- paste0(figures_dir, "/phylum_composition.png")
ggsave(fig_path, last_plot(), width = 10, height = 6, dpi = 300)
cat("Plot saved to:", fig_path, "\n\n")

# Bar plot with samples ordered by diet
plot_bar(physeq_phylum, fill = "Phylum") +
  facet_wrap(~diet, scales="free_x", nrow=1) +
  theme_bw() +
  labs(
    title = "Phylum-level composition by diet group",
    y = "Relative Abundance"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig_path <- paste0(figures_dir, "/phylum_composition_by_diet.png")
ggsave(fig_path, last_plot(), width = 10, height = 6, dpi = 300)
cat("Plot saved to:", fig_path, "\n\n")
