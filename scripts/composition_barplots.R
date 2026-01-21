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


# ==== Barplots of sample composition at phylum-level ====

save_plot <- function(fig_filename) {
  fig_path <- paste0(figures_dir, fig_filename)
  ggsave(fig_path, last_plot(), width = 10, height = 6, dpi = 300)
  cat("Plot saved to:", fig_path, "\n\n")
}

# Aggregate to phylum level, keep the unclassified mOTUs
physeq_phylum <- tax_glom(physeq, taxrank = "Phylum", NArm = FALSE)

# Basic bar plot
plot_bar(physeq_phylum, fill = "Phylum") +
  theme_bw() +
  labs(
    title = "Phylum-level composition",
    y = "Relative Abundance"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot("phylum_composition.png")

# Bar plot with samples grouped by diet
plot_bar(physeq_phylum, fill = "Phylum") +
  facet_wrap(~diet, scales = "free_x", nrow = 1) +
  theme_bw() +
  labs(
    title = "Phylum-level composition by diet group",
    y = "Relative Abundance"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot("phylum_composition_by_diet.png")


# Aggregate to class level
physeq_class <- tax_glom(physeq, taxrank = "Class", NArm = FALSE)
plot_bar(physeq_class, fill = "Class", ) +
  facet_wrap(~diet, scales = "free_x", nrow = 1) +
  theme_bw() +
  labs(
    title = "Class-level composition by diet group",
    y = "Relative Abundance"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot("class_composition_by_diet.png")


# Aggregate to order level
physeq_order <- tax_glom(physeq, taxrank = "Order", NArm = FALSE)
plot_bar(physeq_order, fill = "Order", ) +
  facet_wrap(~diet, scales = "free_x", nrow = 1) +
  theme_bw() +
  labs(
    title = "Order-level composition by diet group",
    y = "Relative Abundance"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot("order_composition_by_diet.png")

# Aggregate to family level
physeq_family <- tax_glom(physeq, taxrank = "Family", NArm = FALSE)
plot_bar(physeq_family, fill = "Family", ) +
  facet_wrap(~diet, scales = "free_x", nrow = 1) +
  theme_bw() +
  labs(
    title = "Family-level composition by diet group",
    y = "Relative Abundance"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# migcer: too large legend, not worth it
# save_plot("family_composition_by_diet.png")

# migcer: too fine-grained to plot all, would need to take top 10
# # Aggregate to genus level
# physeq_genus <- tax_glom(physeq, taxrank = "Genus", NArm = FALSE)
# plot_bar(physeq_genus, fill = "Genus", ) +
#   facet_wrap(~diet, scales = "free_x", nrow = 1) +
#   theme_bw() +
#   labs(
#     title = "Genus-level composition by diet group",
#     y = "Relative Abundance"
#   ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ==== Make family-level composition plot of top 10 families + Other ====
# Keep only top 10 most abundant families
top_families <- names(sort(taxa_sums(physeq_family), decreasing = TRUE)[1:10])
top_families

# Subset to top families
physeq_top_fam <- prune_taxa(top_families, physeq_family)
physeq_top_fam

# Group remaining families as "Other"
# First get all data
ps_melt <- psmelt(physeq_family)

# Identify top phyla
top_fam_names <- ps_melt %>%
  group_by(Family) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  arrange(desc(mean_abundance)) %>%
  top_n(10) %>%
  pull(Family)
top_fam_names

# Label others
ps_melt <- ps_melt %>%
  mutate(
    Family_grouped = ifelse(
      Family %in% top_fam_names, as.character(Family), "Other"
    )
  )
names(ps_melt)

# Convert family groups into factor where Other is at the end
fam_groups <- unique(ps_melt$Family_grouped)
ps_melt$Family_grouped <- factor(
  ps_melt$Family_grouped,
  levels = c(setdiff(sort(fam_groups), "Other"), "Other")
)
glimpse(ps_melt)

# Create custom color palette
library(scales)
n_phyla <- length(fam_groups)
colors <- hue_pal(h = c(0, 360) + 15, c = 100, l = 65)(n_phyla)
colors[n_phyla] <- "#999999"  # Make "Other" gray

# Bar plot with samples ordered by diet
ggplot(ps_melt, aes(x = Sample, y = Abundance, fill = Family_grouped)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colors) +
  facet_grid(~ diet, scales = "free_x", space = "free_x") +
  theme_bw() +
  labs(
    title = "Family-level composition by diet (top 10 families)",
    x = "Sample",
    y = "Relative Abundance",
    fill = "Family"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    strip.background = element_rect(fill = "lightgrey")
  )

save_plot("family_top10_composition_by_diet.png")
