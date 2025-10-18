#!/usr/bin/env Rscript
# Simple R script to generate upset plot showing overlap of significant DEGs
# This script can be run standalone: Rscript generate_upset_plot.R
# Run from the code/09_beautiful_graphics directory

library(tidyverse)
library(UpSetR)

# Load significant gene lists
cat("Loading gene lists...\n")

lrt_genes <- read_csv("../../output/06_deg/interactions/LRT_sig_genes.csv", show_col_types = FALSE) %>%
  pull(gene_id)

interaction_genes <- read_csv("../../output/06_deg/interactions/sig_interaction_genes.csv", show_col_types = FALSE) %>%
  pull(x)

leachate_genes <- read_csv("../../output/09_deg_multifactorial/sig_leachate.csv", show_col_types = FALSE) %>%
  pull(gene)

# Print counts
cat("\nSet sizes:\n")
cat("  LRT significant genes:", length(lrt_genes), "\n")
cat("  Interaction effect genes:", length(interaction_genes), "\n")
cat("  Main leachate effect genes:", length(leachate_genes), "\n")

# Create a list for UpSetR
gene_list <- list(
  LRT = lrt_genes,
  Interaction = interaction_genes,
  Leachate = leachate_genes
)

# Get all unique genes
all_genes <- unique(c(lrt_genes, interaction_genes, leachate_genes))
cat("  Total unique genes:", length(all_genes), "\n")

# Create and save the upset plot
cat("\nGenerating upset plot...\n")
output_file <- "../../output/09_deg_multifactorial/upset_plot_deg_overlap.png"

png(
  filename = output_file,
  width = 10,
  height = 6,
  units = "in",
  res = 300,
  bg = "white"
)

upset(
  fromList(gene_list),
  nsets = 3,
  order.by = "freq",
  sets = c("Leachate", "Interaction", "LRT"),
  keep.order = TRUE,
  point.size = 3.5,
  line.size = 1.5,
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Set Size",
  text.scale = c(1.5, 1.3, 1.2, 1.2, 1.5, 1.2),
  mb.ratio = c(0.6, 0.4)
)

dev.off()

cat("Upset plot saved to:", output_file, "\n")

# Calculate and print overlaps
upset_data_matrix <- fromList(gene_list)

lrt_only <- sum(upset_data_matrix$LRT == 1 & upset_data_matrix$Interaction == 0 & upset_data_matrix$Leachate == 0)
interaction_only <- sum(upset_data_matrix$LRT == 0 & upset_data_matrix$Interaction == 1 & upset_data_matrix$Leachate == 0)
leachate_only <- sum(upset_data_matrix$LRT == 0 & upset_data_matrix$Interaction == 0 & upset_data_matrix$Leachate == 1)
lrt_interaction <- sum(upset_data_matrix$LRT == 1 & upset_data_matrix$Interaction == 1 & upset_data_matrix$Leachate == 0)
lrt_leachate <- sum(upset_data_matrix$LRT == 1 & upset_data_matrix$Interaction == 0 & upset_data_matrix$Leachate == 1)
interaction_leachate <- sum(upset_data_matrix$LRT == 0 & upset_data_matrix$Interaction == 1 & upset_data_matrix$Leachate == 1)
all_three <- sum(upset_data_matrix$LRT == 1 & upset_data_matrix$Interaction == 1 & upset_data_matrix$Leachate == 1)

cat("\nOverlap summary:\n")
cat("  LRT only:", lrt_only, "\n")
cat("  Interaction only:", interaction_only, "\n")
cat("  Leachate only:", leachate_only, "\n")
cat("  LRT & Interaction:", lrt_interaction, "\n")
cat("  LRT & Leachate:", lrt_leachate, "\n")
cat("  Interaction & Leachate:", interaction_leachate, "\n")
cat("  All three (LRT & Interaction & Leachate):", all_three, "\n")

cat("\nDone!\n")
