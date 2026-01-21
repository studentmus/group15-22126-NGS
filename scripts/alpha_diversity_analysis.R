source("scripts/load_data_phyloseq.R")
ls()

library(vegan)


# Remove Taxonomy column and transpose
motus_rel_abundance <- select(motus_rel_abundance, -Taxonomy)
motus_rel_abundance_t <- t(motus_rel_abundance)
dim(motus_rel_abundance_t)

# Calculate Shannon indices
alpha_div <- data.frame(
  sample_id = rownames(motus_rel_abundance_t),
  shannon = diversity(motus_rel_abundance_t, index = "shannon")
)
# Add metadata
alpha_div <- merge(alpha_div, samples, by = "sample_id")
glimpse(alpha_div)

# Do Wilcoxon rank-sum test
wilcox.test(shannon ~ diet, data = alpha_div, exact = FALSE)

# Plot by diet
ggplot(alpha_div, aes(x = diet, y = shannon, fill = diet)) +
  geom_boxplot() +
  labs(
    title = "Alpha diversity (Shannon index) by diet group",
    x = "Diet",
    y = "Shannon index",
  ) +
  theme(legend.position = "none")
