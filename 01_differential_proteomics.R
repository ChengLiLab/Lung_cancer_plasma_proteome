#!/usr/bin/env Rscript

# Differential plasma proteomics: BM vs non-BM
# Input format: first column = UniProt ID; remaining columns = samples.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(ggrepel)
})

set.seed(123)

data_dir <- file.path("data", "plasma")
out_dir <- file.path("results", "differential_proteomics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bm_file <- file.path(data_dir, "BM.xlsx")
non_bm_file <- file.path(data_dir, "NO_BM.xlsx")

log2fc_cutoff <- 0.25
fdr_cutoff <- 0.05
pseudo_count <- 1e-10

read_protein_matrix <- function(path) {
  x <- read_excel(path)
  ids <- as.character(x[[1]])
  mat <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- ids
  mat
}

bm_mat <- read_protein_matrix(bm_file)
non_bm_mat <- read_protein_matrix(non_bm_file)

if (!identical(rownames(bm_mat), rownames(non_bm_mat))) {
  stop("BM and non-BM files must contain proteins in the same row order.")
}

uniprot_ids <- rownames(bm_mat)

gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = uniprot_ids,
  keytype = "UNIPROT",
  column = "SYMBOL",
  multiVals = "first"
)

ensembl_id <- mapIds(
  org.Hs.eg.db,
  keys = uniprot_ids,
  keytype = "UNIPROT",
  column = "ENSEMBL",
  multiVals = "first"
)

bm_mean <- rowMeans(bm_mat, na.rm = TRUE)
non_bm_mean <- rowMeans(non_bm_mat, na.rm = TRUE)
log2fc <- log2(bm_mean / non_bm_mean)

welch_p <- vapply(seq_len(nrow(bm_mat)), function(i) {
  x <- log2(bm_mat[i, ] + pseudo_count)
  y <- log2(non_bm_mat[i, ] + pseudo_count)
  if (sum(is.finite(x)) < 2 || sum(is.finite(y)) < 2) return(NA_real_)
  tryCatch(t.test(x, y, var.equal = FALSE)$p.value, error = function(e) NA_real_)
}, numeric(1))

fdr <- p.adjust(welch_p, method = "BH")

results <- data.frame(
  UniProt = uniprot_ids,
  ENSEMBL = unname(ensembl_id),
  Gene = unname(gene_symbol),
  BM_mean = bm_mean,
  non_BM_mean = non_bm_mean,
  log2FC = log2fc,
  p_value = welch_p,
  FDR = fdr,
  stringsAsFactors = FALSE
) %>%
  mutate(
    status = case_when(
      FDR < fdr_cutoff & log2FC > log2fc_cutoff ~ "Up",
      FDR < fdr_cutoff & log2FC < -log2fc_cutoff ~ "Down",
      TRUE ~ "Not significant"
    )
  )

write.csv(results, file.path(out_dir, "all_proteins_statistics.csv"), row.names = FALSE)

dep <- results %>% filter(status != "Not significant")
write.csv(dep, file.path(out_dir, "differentially_abundant_proteins.csv"), row.names = FALSE)

top_labels <- bind_rows(
  dep %>% filter(status == "Up") %>% arrange(desc(log2FC)) %>% slice_head(n = 5),
  dep %>% filter(status == "Down") %>% arrange(log2FC) %>% slice_head(n = 5)
)

volcano <- results %>%
  mutate(minus_log10_FDR = -log10(pmax(FDR, .Machine$double.xmin)))

p <- ggplot(volcano, aes(log2FC, minus_log10_FDR, color = status)) +
  geom_point(alpha = 0.7, size = 1.6) +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed") +
  geom_text_repel(
    data = top_labels %>% mutate(minus_log10_FDR = -log10(pmax(FDR, .Machine$double.xmin))),
    aes(label = ifelse(is.na(Gene) | Gene == "", UniProt, Gene)),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(x = "log2 fold change (BM / non-BM)", y = "-log10(FDR)", color = NULL) +
  theme_classic()

ggsave(file.path(out_dir, "volcano_BM_vs_nonBM.pdf"), p, width = 6.5, height = 5.5)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
