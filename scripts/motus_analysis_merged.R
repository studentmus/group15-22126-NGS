## ======================================================
## mOTUs metagenomics analysis
## Input: motus_merged.tsv (counts)
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
motus_file <- "data/motus_merged.tsv"
meta_file  <- "data/samples.csv"
outdir     <- "results"

dir.create(outdir, showWarnings = FALSE)

## =====================================================
## 1. Load mOTUs merged table
## =====================================================

motus_raw <- read.table(
  motus_file,
  header = TRUE,
  sep = "\t",
  comment.char = "#",   # skip tool_version line
  check.names = FALSE,
  stringsAsFactors = FALSE
)

## First column = mOTU ID
rownames(motus_raw) <- motus_raw[, 1]

## Drop mOTU ID + Taxonomy columns
motus <- motus_raw[, -(1:2)]

## Clean sample names (path -> SRR ID)
colnames(motus) <- basename(colnames(motus))
colnames(motus) <- sub("_trimmed_nonhuman.motus.name$", "", colnames(motus))

## Convert to numeric matrix
motus <- as.matrix(motus)
storage.mode(motus) <- "numeric"

## Remove mOTUs with zero counts across all samples
motus <- motus[rowSums(motus) > 0, ]

cat("mOTUs:", nrow(motus), "\n")
cat("Samples:", ncol(motus), "\n")

## Safety check
if (ncol(motus) < 2) {
  stop("Not enough samples for multivariate analysis")
}

## =====================================================
## 2. Load and synchronise metadata
## =====================================================

meta <- read.csv(
  meta_file,
  header = FALSE,
  stringsAsFactors = FALSE
)

## samples.csv: column 1 = sample, column 2 = group
colnames(meta) <- c("sample", "group")

## Reorder metadata to match motus columns
sample_df <- data.frame(
  sample = colnames(motus),
  stringsAsFactors = FALSE
)

meta <- dplyr::left_join(sample_df, meta, by = "sample")

## Sanity checks
stopifnot(!any(is.na(meta$group)))
stopifnot(all(meta$sample == colnames(motus)))

meta$group <- trimws(meta$group)

cat("Samples retained for analysis:\n")
print(meta)

## =====================================================
## 3. Alpha diversity (richness & Shannon)
## =====================================================

motus_t <- t(motus)

alpha_df <- data.frame(
  Sample   = rownames(motus_t),
  Richness = specnumber(motus_t),
  Shannon  = diversity(motus_t, index = "shannon"),
  Group    = meta$group
)

p_rich <- ggplot(alpha_df, aes(Group, Richness, fill = Group)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  labs(title = "Species richness per sample")

ggsave(
  file.path(outdir, "alpha_richness.png"),
  p_rich,
  width = 6,
  height = 5
)

p_shan <- ggplot(alpha_df, aes(Group, Shannon, fill = Group)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  labs(title = "Shannon diversity per sample")

ggsave(
  file.path(outdir, "alpha_shannon.png"),
  p_shan,
  width = 6,
  height = 5
)

## =====================================================
## 4. PCA (community composition)
## =====================================================

pca <- prcomp(motus_t, scale. = TRUE)

pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1    = pca$x[, 1],
  PC2    = pca$x[, 2],
  Group  = meta$group
)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(
    title = "PCA of mOTUs community composition",
    x = "PC1",
    y = "PC2"
  )

ggsave(
  file.path(outdir, "PCA_motus.png"),
  p_pca,
  width = 6,
  height = 5
)

## =====================================================
## 5. Beta diversity (Bray–Curtis + PCoA)
## =====================================================

bray <- vegdist(motus_t, method = "bray")
pcoa <- cmdscale(bray, k = 2, eig = TRUE)

pcoa_df <- data.frame(
  Sample = rownames(motus_t),
  PCoA1  = pcoa$points[, 1],
  PCoA2  = pcoa$points[, 2],
  Group  = meta$group
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
## 6. Exploratory differential abundance (mean comparison)
## =====================================================

## Explicit group labels
g1 <- unique(meta$group)[1]
g2 <- unique(meta$group)[2]

idx_g1 <- which(meta$group == g1)
idx_g2 <- which(meta$group == g2)

stopifnot(length(idx_g1) > 0, length(idx_g2) > 0)

mean_g1 <- rowMeans(motus[, idx_g1, drop = FALSE])
mean_g2 <- rowMeans(motus[, idx_g2, drop = FALSE])

diff_df <- data.frame(
  mOTU        = rownames(motus),
  Mean_group1 = mean_g1,
  Mean_group2 = mean_g2,
  Difference  = mean_g1 - mean_g2,
  stringsAsFactors = FALSE
)

top_motus <- diff_df %>%
  arrange(desc(abs(Difference))) %>%
  head(10)

write.csv(
  top_motus,
  file.path(outdir, "top_differential_motus.csv"),
  row.names = FALSE
)

## Boxplot of top mOTUs
top_mat <- motus[top_motus$mOTU, , drop = FALSE]

top_long <- data.frame(
  Abundance = as.vector(top_mat),
  mOTU      = rep(rownames(top_mat), times = ncol(top_mat)),
  Group     = rep(meta$group, each = nrow(top_mat)),
  stringsAsFactors = FALSE
)

p_box <- ggplot(top_long, aes(Group, Abundance, fill = Group)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ mOTU, scales = "free_y") +
  theme_bw() +
  labs(title = "Top differentially abundant mOTUs (exploratory)")

ggsave(
  file.path(outdir, "top_mOTUs_boxplot.png"),
  p_box,
  width = 10,
  height = 6
)

## =====================================================
## 7. Session info
## =====================================================

writeLines(
  capture.output(sessionInfo()),
  file.path(outdir, "sessionInfo.txt")
)

cat("Analysis completed successfully.\n")
