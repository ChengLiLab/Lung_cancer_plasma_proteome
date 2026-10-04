#!/usr/bin/env Rscript

# Human tissue expression context and tissue-set enrichment of prioritized proteins.
# Tissue.xlsx format:
#   column 1 = FileName
#   column 2 = Major tissue type
#   remaining columns = "Gene ID / UniProt" protein-abundance columns
# Proteins.xlsx: first column contains target gene IDs.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
})

set.seed(123)

data_dir <- file.path("data", "tissue")
out_dir <- file.path("results", "tissue_expression")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

proteins_file <- file.path(data_dir, "Proteins.xlsx")
tissue_file <- file.path(data_dir, "Tissue.xlsx")

pseudo_count <- 1e-10
n_permutations <- 5000

target_genes <- read_excel(proteins_file)[[1]] %>%
  as.character() %>%
  str_trim() %>%
  .[!is.na(.) & . != ""] %>%
  unique()

raw <- read_excel(tissue_file)
if (ncol(raw) < 3) stop("Tissue.xlsx must contain two metadata columns and protein columns.")

metadata <- raw[, 1:2]
names(metadata) <- c("FileName", "Tissue")

expr_raw <- raw[, -(1:2), drop = FALSE]
protein_map <- tibble(
  original_col = names(expr_raw),
  GeneID = str_trim(str_extract(names(expr_raw), "^[^/]+")),
  UniProt = ifelse(
    str_detect(names(expr_raw), "/"),
    str_trim(str_extract(names(expr_raw), "(?<=/).*")),
    NA_character_
  )
)

expr_raw[] <- lapply(expr_raw, function(x) suppressWarnings(as.numeric(str_replace_all(as.character(x), ",", ""))))

# Collapse multiple columns mapping to the same gene by sample-level mean.
expr_long <- bind_cols(metadata, expr_raw) %>%
  pivot_longer(-c(FileName, Tissue), names_to = "original_col", values_to = "expression") %>%
  left_join(protein_map, by = "original_col")

expr_gene <- expr_long %>%
  group_by(FileName, Tissue, GeneID) %>%
  summarise(
    expression = ifelse(all(is.na(expression)), NA_real_, mean(expression, na.rm = TRUE)),
    .groups = "drop"
  )

expr_wide <- expr_gene %>%
  select(FileName, Tissue, GeneID, expression) %>%
  pivot_wider(names_from = GeneID, values_from = expression)

meta <- expr_wide %>% select(FileName, Tissue)
expr_mat <- expr_wide %>% select(-FileName, -Tissue) %>% as.data.frame()
rownames(expr_mat) <- meta$FileName
expr_mat <- as.matrix(expr_mat)
storage.mode(expr_mat) <- "numeric"

matched <- intersect(target_genes, colnames(expr_mat))
unmatched <- setdiff(target_genes, colnames(expr_mat))
write.csv(data.frame(Gene = matched), file.path(out_dir, "matched_target_genes.csv"), row.names = FALSE)
write.csv(data.frame(Gene = unmatched), file.path(out_dir, "unmatched_target_genes.csv"), row.names = FALSE)

if (length(matched) < 2) stop("Fewer than two target genes matched Tissue.xlsx.")

log_expr <- log2(expr_mat + pseudo_count)

z_by_gene <- function(m) {
  out <- apply(m, 2, function(x) {
    if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) return(rep(NA_real_, length(x)))
    as.numeric(scale(x))
  })
  out <- as.matrix(out)
  rownames(out) <- rownames(m)
  colnames(out) <- colnames(m)
  out
}

row_mean_na <- function(m) {
  apply(m, 1, function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE))
}
col_mean_na <- function(m) {
  apply(m, 2, function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE))
}

z_expr <- z_by_gene(log_expr)
matched <- intersect(matched, colnames(z_expr)[colSums(!is.na(z_expr)) > 0])

# Tissue-level heatmap for target genes.
tissue_mean <- expr_gene %>%
  filter(GeneID %in% matched) %>%
  mutate(log2_expression = log2(expression + pseudo_count)) %>%
  group_by(Tissue, GeneID) %>%
  summarise(mean_log2_expression = mean(log2_expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Tissue, values_from = mean_log2_expression)

heat_mat <- tissue_mean %>%
  column_to_rownames("GeneID") %>%
  as.matrix()

heat_mat_z <- t(scale(t(heat_mat)))
pdf(file.path(out_dir, "target_protein_tissue_heatmap.pdf"), width = 10, height = max(5, nrow(heat_mat_z) * 0.25))
pheatmap(heat_mat_z, cluster_rows = TRUE, cluster_cols = TRUE, border_color = NA)
dev.off()

# Sample-level target-set score.
target_z <- z_expr[, matched, drop = FALSE]
score_df <- meta %>%
  mutate(
    GeneSetScore = row_mean_na(target_z),
    DetectedTargetGenes = rowSums(
      !is.na(expr_mat[, matched, drop = FALSE]) &
        expr_mat[, matched, drop = FALSE] > 0
    )
  ) %>%
  filter(!is.na(GeneSetScore), !is.na(Tissue), Tissue != "")

write.csv(score_df, file.path(out_dir, "sample_target_set_scores.csv"), row.names = FALSE)

# One-sided tissue-vs-rest Wilcoxon tests.
tissues <- sort(unique(score_df$Tissue))
one_vs_rest <- lapply(tissues, function(tt) {
  x <- score_df$GeneSetScore[score_df$Tissue == tt]
  y <- score_df$GeneSetScore[score_df$Tissue != tt]
  p <- tryCatch(wilcox.test(x, y, alternative = "greater")$p.value, error = function(e) NA_real_)
  data.frame(
    Tissue = tt,
    N = length(x),
    Mean_score = mean(x, na.rm = TRUE),
    Median_score = median(x, na.rm = TRUE),
    P_value = p
  )
}) %>% bind_rows() %>%
  mutate(FDR = p.adjust(P_value, method = "BH")) %>%
  arrange(FDR, desc(Mean_score))

write.csv(one_vs_rest, file.path(out_dir, "tissue_one_vs_rest_wilcox.csv"), row.names = FALSE)

# Permutation enrichment using tissue-level mean z-scores.
tissue_gene_z <- lapply(tissues, function(tt) {
  idx <- which(meta$Tissue == tt)
  v <- col_mean_na(z_expr[idx, , drop = FALSE])
  data.frame(Tissue = tt, GeneID = names(v), MeanZ = as.numeric(v))
}) %>% bind_rows() %>%
  pivot_wider(names_from = GeneID, values_from = MeanZ)

tissue_names <- tissue_gene_z$Tissue
tissue_gene_mat <- tissue_gene_z %>% select(-Tissue) %>% as.data.frame()
rownames(tissue_gene_mat) <- tissue_names
tissue_gene_mat <- as.matrix(tissue_gene_mat)
storage.mode(tissue_gene_mat) <- "numeric"

background <- colnames(tissue_gene_mat)[colSums(!is.na(tissue_gene_mat)) > 0]
matched <- intersect(matched, background)

observed <- row_mean_na(tissue_gene_mat[, matched, drop = FALSE])

null_scores <- replicate(n_permutations, {
  random_genes <- sample(background, length(matched), replace = FALSE)
  row_mean_na(tissue_gene_mat[, random_genes, drop = FALSE])
})

empirical_p <- vapply(seq_along(observed), function(i) {
  if (is.na(observed[i])) return(NA_real_)
  (1 + sum(null_scores[i, ] >= observed[i], na.rm = TRUE)) / (n_permutations + 1)
}, numeric(1))

null_mean <- rowMeans(null_scores, na.rm = TRUE)
null_sd <- apply(null_scores, 1, sd, na.rm = TRUE)
enrichment_z <- (observed - null_mean) / null_sd
enrichment_z[!is.finite(enrichment_z)] <- NA_real_

perm <- data.frame(
  Tissue = rownames(tissue_gene_mat),
  Observed_set_score = observed,
  Null_mean = null_mean,
  Null_sd = null_sd,
  Enrichment_Z = enrichment_z,
  Empirical_P = empirical_p,
  Target_genes_used = rowSums(!is.na(tissue_gene_mat[, matched, drop = FALSE]))
) %>%
  mutate(
    FDR = p.adjust(Empirical_P, method = "BH"),
    Significant_enrichment = Observed_set_score > 0 & Enrichment_Z > 0 & FDR < 0.05
  ) %>%
  arrange(FDR, desc(Enrichment_Z))

write.csv(perm, file.path(out_dir, "tissue_permutation_enrichment.csv"), row.names = FALSE)
write.csv(
  perm %>% filter(Significant_enrichment),
  file.path(out_dir, "significant_tissue_enrichment.csv"),
  row.names = FALSE
)

p <- perm %>%
  filter(Significant_enrichment) %>%
  ggplot(aes(x = reorder(Tissue, Enrichment_Z), y = Enrichment_Z)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "Enrichment Z") +
  theme_classic()

ggsave(file.path(out_dir, "significant_tissue_enrichment.pdf"), p, width = 6.5, height = 4.5)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
