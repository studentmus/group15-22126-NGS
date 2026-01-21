# Some background reading:
# - https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2017.02224/full
# - https://microbiome.github.io/OMA/docs/devel/pages/community_similarity.html
# PERMANOVA is based on group centroids, and it assumes that within-group
# variation is smaller than between-group variation. We need to perform PERMDISP
# together with PERMANOVA to investigate this.

# Goal:
# - Test if the diet groups are significant in explaining the variation between
# samples
#
# Plan:
# - Use Bray-Curtis distance - prone to compositionality bias
# - Use Aitchison distance (CLR-transformed data) to avoid compositionality bias
# - Look at dispersion within groups (PERMDISP)
# - Look at dispersion and centroids between groups (PERMANOVA)
# (e.g., test if groups have same dispersion and centroids)

# Data loading and transformations
source("scripts/load_data_phyloseq.R")
ls()
# Keep only the useful data objects
remove(list = setdiff(ls(), c("motus_rel_abundance", "samples", "physeq")))
ls()

library(phyloseq)
library(vegan)
library(compositions)
library(tidyverse)


# Helper function to save figures
save_plot <- function(fig_filename) {
  figures_dir <- "data/figures/"
  fig_path <- paste0(figures_dir, fig_filename)
  ggsave(fig_path, last_plot(), width = 10, height = 6, dpi = 300)
  cat("Plot saved to:", fig_path, "\n\n")
}

# ==== V1 - using Bray-Curtis dissimilarity, species-level =====

# Remove Taxonomy column and transpose
motus_rel_abundance <- select(motus_rel_abundance, -Taxonomy)
motus_rel_abundance_t <- t(motus_rel_abundance)

# Bray-Curtis dissimilarity
bray <- vegdist(motus_rel_abundance_t, method = "bray")

# PERMANOVA
# Prepare metadata - rows in the same order as in the distance matrix
metadata <- data.frame(
  sample_id = rownames(motus_rel_abundance_t),
  row.names = rownames(motus_rel_abundance_t)
)
metadata$diet <- samples$diet[match(metadata$sample_id, samples$sample_id)]

adonis_res <- adonis2(bray ~ diet, data = metadata, permutations = 999)
adonis_res

# PERMDIST (Beta-dispersion, e.g. within-group variation)
disp <- betadisper(bray, metadata$diet)
permutest(disp, permutations = 999)

# Interpretation:
# - Significant result -> within-group variation differs
# - Non-significant result -> within-group variation does not differ

# Visualise within-group variation
boxplot(
  disp,
  xlab = "Diet",
  ylab = "Distance to group centroid",
  main = "Within-group dispersion (Bray-Curtis)"
)

# ==== V2 - using Bray-Curtis dissimilarity, species-level, phyloseq =====

bray_ps <- phyloseq::distance(physeq, method = "bray")

# PERMANOVA
metadata_ps <- data.frame(sample_data(physeq))
adonis2(
  bray_ps ~ diet,
  data = metadata_ps,
  permutations = 999
)

# PERMDIST (Beta-dispersion, e.g. within-group variation)
disp_ps <- betadisper(bray_ps, metadata_ps$diet)
permutest(disp_ps)

# Visualise within-group variation
boxplot(
  disp_ps,
  xlab = "Diet",
  ylab = "Distance to group centroid",
  main = "Within-group dispersion (Bray-Curtis)"
)

# ==== V3 - Use Aitchison distance, species-level ====

# Transpose (need samples in rows), and ensure no zeroes
otu_table_clr_t <- t(otu_table(physeq))
otu_table_clr_t[otu_table_clr_t == 0] <- 1e-6
dim(otu_table_clr_t)

# CLR transformation (`clr` applies transformation to input matrix rows)
otu_table_clr_t <- clr(otu_table_clr_t)
# otu_table_clr <- t(apply(otu_table_clr_t, 1, function(x) clr(x)))
dim(otu_table_clr_t)

# Compute Aitchison distance (Euclidean dist on CLR)
aitchison_dist <- dist(otu_table_clr_t, method = "euclidean")
aitchison_dist

# PERMANOVA
adonis_res_ait <- adonis2(
  aitchison_dist ~ diet,
  data = metadata_ps,
  permutations = 999
)
adonis_res_ait

# PERMDIST (Beta-dispersion, e.g. within-group variation)
disp_ait <- betadisper(aitchison_dist, metadata_ps$diet)
permutest(disp_ait)

# Visualise within-group variation
boxplot(
  disp_ait,
  xlab = "Diet",
  ylab = "Distance to group centroid",
  main = "Within-group dispersion (Aitchison distance)"
)

# Less grim plot
df_disp <- data.frame(
  Diet = disp_ait$group,
  Distance = disp_ait$distances
)

ggplot(df_disp, aes(x = Diet, y = Distance, fill = Diet)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.08, size = 2, alpha = 0.8) +
  # scale_fill_manual(values = c("vegan" = "#66c2a5", "omnivore" = "#fc8d62")) +
  labs(
    title = "Within-group dispersion (Aitchison distance)",
    x = "Diet",
    y = "Distance to group centroid"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
save_plot("within_group_variation_aitchinson.png")


# ==== PCA to explore diversity between diet groups, species-level ====

pca <- prcomp(otu_table_clr_t, center = TRUE, scale. = FALSE)

# Sanity check: distances in PCA should match Aitchkinson distances
# - Correlation coef should be close to 1
d1 <- dist(otu_table_clr_t)
d2 <- dist(pca$x)
cor(as.vector(d1), as.vector(d2))

# Plot first two PCs
pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Diet = metadata_ps$diet
)
variance_explained <- summary(pca)$importance[2, 1:2] * 100

ggplot(pca_df, aes(x = PC1, y = PC2, color = Diet)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCA of CLR-transformed mOTU abundances",
    x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 1), "%)")
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    # legend.title = element_blank()
  )
save_plot("pca_clr_transformed.png")

# ==== V4 - Use Aitchison distance, family-level ====

# Aggregate at family level
physeq_family <- tax_glom(physeq, taxrank = "Family", NArm = FALSE)
metadata_ps_fam <- data.frame(sample_data(physeq_family))

# Transpose (need samples in rows), and ensure no zeroes
otu_table_fam_clr_t <- t(otu_table(physeq_family))
otu_table_fam_clr_t[otu_table_fam_clr_t == 0] <- 1e-6
dim(otu_table_fam_clr_t)

# CLR transformation (`clr` applies transformation to input matrix rows)
otu_table_fam_clr_t <- clr(otu_table_fam_clr_t)
dim(otu_table_fam_clr_t)

# Compute Aitchison distance (Euclidean dist on CLR)
aitchison_dist_fam <- dist(otu_table_fam_clr_t, method = "euclidean")
aitchison_dist_fam

# PERMANOVA
adonis_res_fam <- adonis2(
  aitchison_dist_fam ~ diet,
  data = metadata_ps_fam,
  permutations = 999
)
adonis_res_fam

# PERMDIST (Beta-dispersion, e.g. within-group variation)
disp_fam <- betadisper(aitchison_dist_fam, metadata_ps_fam$diet)
permutest(disp_fam)

# Visualise within-group variation
df_disp_fam <- data.frame(
  Diet = disp_fam$group,
  Distance = disp_fam$distances
)

ggplot(df_disp_fam, aes(x = Diet, y = Distance, fill = Diet)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.08, size = 2, alpha = 0.8) +
  # scale_fill_manual(values = c("vegan" = "#66c2a5", "omnivore" = "#fc8d62")) +
  labs(
    title = "Within-group dispersion (Aitchison distance) at Family-level",
    x = "Diet",
    y = "Distance to group centroid"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
save_plot("within_group_variation_aitchinson_family.png")
