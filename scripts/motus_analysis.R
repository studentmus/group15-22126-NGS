## ======================================================
## mOTUs metagenomics analysis – FINAL STABLE VERSION
## ======================================================

rm(list = ls())

## ----------------------------
## Libraries
## ----------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(vegan)
})

## ----------------------------
## Paths
## ----------------------------
motus_file <- "data/motus_species_relab.strict.tsv"
meta_file  <- "data/samples.csv"
outdir     <- "results"

dir.create(outdir, showWarnings = FALSE)

## ----------------------------
## 1. Load mOTUs table
## ----------------------------
motus <- read.table(
  motus_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE,
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

## First column is taxon
stopifnot(colnames(motus)[1] == "taxon")

## ------------------------------------------------
## Merge duplicated taxa (mOTUs-specific requirement)
## ------------------------------------------------

motus_df <- as.data.frame(motus)

motus_df$taxon <- rownames(motus_df)

motus_df <- motus_df %>%
  group_by(taxon) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")

rownames(motus_df) <- motus_df$taxon
motus_df$taxon <- NULL

motus <- as.matrix(motus_df)
storage.mode(motus) <- "numeric"


## Convert to numeric matrix
motus <- as.matrix(motus)
storage.mode(motus) <- "numeric"
## ------------------------------------------------
## Remove non-abundance (taxonomy) rows
## ------------------------------------------------

# Keep only rows where at least one sample has a numeric value > 0
keep_rows <- apply(motus, 1, function(x) {
  any(is.finite(x) & x > 0)
})

motus <- motus[keep_rows, ]

cat("After removing taxonomy rows:\n")
cat("Number of species:", nrow(motus), "\n")
cat("Number of samples:", ncol(motus), "\n")


cat("Number of species:", nrow(motus), "\n")
cat("Number of samples:", ncol(motus), "\n")

## Safety check
if (ncol(motus) < 2) {
  stop("Not enough samples for multivariate analysis")
}

## Remove species with zero abundance
motus <- motus[rowSums(motus) > 0, ]

## ----------------------------
## 2. Load metadata
## ----------------------------
## Synchronise metadata
## ----------------------------

## ----------------------------
## Add headers
## ----------------------------

meta <- read.csv(
  "data/samples.csv",
  header = FALSE,
  stringsAsFactors = FALSE
)

# samples.csv has no header: column 1 = sample, column 2 = group
colnames(meta) <- c("sample", "group")


# Create a data.frame of samples from motus matrix
sample_df <- data.frame(
  sample = colnames(motus),
  stringsAsFactors = FALSE
)

# Left join metadata to motus samples
meta <- dplyr::left_join(sample_df, meta, by = "sample")

# Final sanity check
stopifnot(!any(is.na(meta$group)))
stopifnot(all(meta$sample == colnames(motus)))

cat("Samples retained for analysis:\n")
print(meta)

meta$group <- trimws(meta$group)


## ----------------------------
## 3. PCA analysis
## ----------------------------
motus_t <- t(motus)

pca <- prcomp(motus_t, scale. = TRUE)

pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = meta$group
)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(
    title = "PCA of microbial community composition",
    x = "PC1",
    y = "PC2"
  )

ggsave(
  file.path(outdir, "PCA_motus.png"),
  p_pca,
  width = 6,
  height = 5
)

## ----------------------------
## 4. Beta diversity (Bray–Curtis, PCoA)
## ----------------------------
bray <- vegdist(motus_t, method = "bray")
pcoa <- cmdscale(bray, k = 2, eig = TRUE)

pcoa_df <- data.frame(
  Sample = rownames(motus_t),
  PCoA1 = pcoa$points[, 1],
  PCoA2 = pcoa$points[, 2],
  Group = meta$group
)

p_pcoa <- ggplot(pcoa_df, aes(PCoA1, PCoA2, color = Group)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(
    title = "PCoA (Bray–Curtis distance)",
    x = "PCoA1",
    y = "PCoA2"
  )

ggsave(
  file.path(outdir, "PCoA_bray_motus.png"),
  p_pcoa,
  width = 6,
  height = 5
)

## =====================================================
## 5. Exploratory differential abundance (FINAL)
## =====================================================

# Explicit group names
g1 <- "vegan"
g2 <- "omnivore"

# Clean group labels just in case
meta$group <- trimws(meta$group)

# Get column indices (THIS IS THE KEY FIX)
idx_g1 <- which(meta$group == g1)
idx_g2 <- which(meta$group == g2)

# Sanity checks
stopifnot(length(idx_g1) > 0)
stopifnot(length(idx_g2) > 0)
stopifnot(max(c(idx_g1, idx_g2)) <= ncol(motus))

# Compute group means
mean_g1 <- rowMeans(motus[, idx_g1, drop = FALSE])
mean_g2 <- rowMeans(motus[, idx_g2, drop = FALSE])

# Build differential abundance table
diff_df <- data.frame(
  Species = seq_len(nrow(motus)),
  Mean_vegan = mean_g1,
  Mean_omnivore = mean_g2,
  Difference = mean_g1 - mean_g2,
  stringsAsFactors = FALSE
)


# Top species by absolute difference
top_species <- diff_df %>%
  arrange(desc(abs(Difference))) %>%
  head(10)

top_taxa <- top_species$Species

# Write results
write.csv(
  top_species,
  file.path(outdir, "top_differential_species.csv"),
  row.names = FALSE
)


## Boxplot for top species
top_mat <- motus[top_taxa, , drop = FALSE]
rownames(top_mat) <- top_taxa

top_long <- data.frame(
  Abundance = as.vector(top_mat),
  Species = rep(top_taxa, times = ncol(top_mat)),
  Group = rep(meta$group, each = length(top_taxa)),
  stringsAsFactors = FALSE
)

p_box <- ggplot(top_long, aes(Group, Abundance, fill = Group)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ Species, scales = "free_y") +
  theme_bw() +
  labs(title = "Top differentially abundant species")

ggsave(
  file.path(outdir, "top_species_boxplot.png"),
  p_box,
  width = 10,
  height = 6
)

## ----------------------------
## 6. Session info
## ----------------------------
writeLines(
  capture.output(sessionInfo()),
  file.path(outdir, "sessionInfo.txt")
)

cat("Analysis completed successfully.\n")

