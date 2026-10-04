#!/usr/bin/env Rscript

# WGCNA of the plasma proteome
# The uploaded workflow selected the 3,500 most variable proteins by MAD,
# removed sample "[105]", and used a signed bicor network with power = 10.

suppressPackageStartupMessages({
  library(readxl)
  library(WGCNA)
  library(dplyr)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()
set.seed(123)

data_dir <- file.path("data", "plasma")
out_dir <- file.path("results", "wgcna")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(data_dir, "All_sample_median_normalized_log2.xlsx")
metadata_file <- file.path(data_dir, "sample_metadata.csv")

n_top <- 3500
outlier_sample <- "[105]"
soft_power <- 10
deep_split <- 4
min_module_size <- 20
merge_cut_height <- 0.10

expr_raw <- read_excel(input_file)
protein_id <- as.character(expr_raw[[1]])

expr_mat <- as.data.frame(expr_raw[, -1, drop = FALSE])
rownames(expr_mat) <- make.unique(protein_id)
expr_mat[] <- lapply(expr_mat, function(x) as.numeric(as.character(x)))

if (nrow(expr_mat) < n_top) stop("Fewer than 3,500 proteins are available.")

if (file.exists(metadata_file)) {
  sample_info <- read.csv(metadata_file, stringsAsFactors = FALSE)
  if (!all(c("Sample", "Disease_status") %in% names(sample_info))) {
    stop("sample_metadata.csv must contain Sample and Disease_status columns.")
  }
  rownames(sample_info) <- sample_info$Sample
  if (!all(colnames(expr_mat) %in% rownames(sample_info))) {
    stop("Some expression-matrix samples are missing from sample_metadata.csv.")
  }
  sample_info <- sample_info[colnames(expr_mat), , drop = FALSE]
} else {
  if (ncol(expr_mat) != 110) {
    stop("No sample_metadata.csv found and the expression matrix does not contain 110 samples.")
  }
  sample_info <- data.frame(
    Sample = colnames(expr_mat),
    Disease_status = c(rep("non-BM", 64), rep("BM", 46)),
    stringsAsFactors = FALSE
  )
  rownames(sample_info) <- sample_info$Sample
}

protein_mad <- apply(expr_mat, 1, mad, na.rm = TRUE)
top_proteins <- names(sort(protein_mad, decreasing = TRUE))[seq_len(n_top)]
datExpr <- as.data.frame(t(expr_mat[top_proteins, , drop = FALSE]))

gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  sample_info <- sample_info[rownames(datExpr), , drop = FALSE]
}

pdf(file.path(out_dir, "sample_clustering.pdf"), width = 10, height = 6)
plot(hclust(dist(datExpr)), main = "Sample clustering", xlab = "", sub = "")
dev.off()

if (!outlier_sample %in% rownames(datExpr)) {
  stop(sprintf("Expected outlier sample '%s' was not found. Edit outlier_sample if needed.", outlier_sample))
}

datExpr <- datExpr[rownames(datExpr) != outlier_sample, , drop = FALSE]
sample_info <- sample_info[rownames(datExpr), , drop = FALSE]

write.csv(datExpr, file.path(out_dir, "wgcna_input_top3500_after_outlier_removal.csv"))

powers <- c(1:20)
sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  corFnc = "bicor",
  corOptions = list(use = "pairwise.complete.obs", maxPOutliers = 0.05),
  verbose = 0
)

pdf(file.path(out_dir, "soft_threshold_selection.pdf"), width = 9, height = 4.5)
par(mfrow = c(1, 2))
plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft-threshold power",
  ylab = "Signed R^2",
  type = "n"
)
text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers,
  cex = 0.8
)
abline(h = 0.8, lty = 2)
plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft-threshold power",
  ylab = "Mean connectivity",
  type = "n"
)
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.8)
dev.off()

net <- blockwiseModules(
  datExpr,
  power = soft_power,
  networkType = "signed",
  TOMType = "signed",
  corType = "bicor",
  maxPOutliers = 0.05,
  deepSplit = deep_split,
  minModuleSize = min_module_size,
  reassignThreshold = 0,
  mergeCutHeight = merge_cut_height,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 2
)

module_labels <- net$colors
module_colors <- labels2colors(module_labels)
names(module_colors) <- colnames(datExpr)

MEs <- orderMEs(net$MEs)

protein_module <- data.frame(
  Protein = colnames(datExpr),
  Module_label = module_labels,
  Module_color = module_colors,
  stringsAsFactors = FALSE
)

write.csv(protein_module, file.path(out_dir, "protein_module_assignment.csv"), row.names = FALSE)
write.csv(MEs, file.path(out_dir, "module_eigengenes.csv"))

pdf(file.path(out_dir, "gene_dendrogram_and_modules.pdf"), width = 12, height = 6)
plotDendroAndColors(
  net$dendrograms[[1]],
  module_colors[net$blockGenes[[1]]],
  "Module",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
dev.off()

trait_data <- data.frame(
  non_BM = as.integer(sample_info$Disease_status == "non-BM"),
  BM = as.integer(sample_info$Disease_status == "BM"),
  row.names = rownames(sample_info)
)
trait_data <- trait_data[rownames(MEs), , drop = FALSE]

module_trait_cor <- cor(MEs, trait_data, use = "pairwise.complete.obs", method = "pearson")
module_trait_p <- corPvalueStudent(module_trait_cor, nSamples = nrow(MEs))
module_trait_fdr <- matrix(
  p.adjust(as.vector(module_trait_p), method = "BH"),
  nrow = nrow(module_trait_p),
  ncol = ncol(module_trait_p),
  dimnames = dimnames(module_trait_p)
)

write.csv(module_trait_cor, file.path(out_dir, "module_trait_correlation.csv"))
write.csv(module_trait_p, file.path(out_dir, "module_trait_pvalue.csv"))
write.csv(module_trait_fdr, file.path(out_dir, "module_trait_FDR.csv"))

cor_long <- as.data.frame(as.table(module_trait_cor), stringsAsFactors = FALSE)
p_long <- as.data.frame(as.table(module_trait_p), stringsAsFactors = FALSE)
fdr_long <- as.data.frame(as.table(module_trait_fdr), stringsAsFactors = FALSE)
names(cor_long) <- c("Module", "Trait", "Correlation")
names(p_long) <- c("Module", "Trait", "Pvalue")
names(fdr_long) <- c("Module", "Trait", "FDR")

module_trait_results <- cor_long %>%
  left_join(p_long, by = c("Module", "Trait")) %>%
  left_join(fdr_long, by = c("Module", "Trait")) %>%
  arrange(Trait, FDR, desc(abs(Correlation)))

write.csv(module_trait_results, file.path(out_dir, "module_trait_results_long.csv"), row.names = FALSE)

text_matrix <- matrix(
  sprintf("%.2f\nFDR=%.3g", module_trait_cor, module_trait_fdr),
  nrow = nrow(module_trait_cor),
  dimnames = dimnames(module_trait_cor)
)

pdf(file.path(out_dir, "module_trait_heatmap.pdf"), width = 6, height = max(5, nrow(MEs) * 0.35))
labeledHeatmap(
  Matrix = module_trait_cor,
  xLabels = colnames(trait_data),
  yLabels = rownames(module_trait_cor),
  ySymbols = rownames(module_trait_cor),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = text_matrix,
  setStdMargins = FALSE,
  cex.text = 0.7,
  zlim = c(-1, 1),
  main = "Module-trait relationships"
)
dev.off()

saveRDS(net, file.path(out_dir, "wgcna_network.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
