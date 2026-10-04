#!/usr/bin/env Rscript

# Integrated single-cell analysis for lung cancer primary lesions and bone metastases.
# This script consolidates the uploaded single-cell workflow:
# QC/integration -> marker identification -> cell annotation -> cell composition ->
# candidate-gene expression -> Milo differential abundance -> CellChat.

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(Matrix)
  library(SingleCellExperiment)
  library(MiloR)
  library(ggbeeswarm)
  library(CellChat)
})

set.seed(20241106)

config_dir <- "config"
out_dir <- file.path("results", "single_cell")
obj_dir <- file.path(out_dir, "objects")
tab_dir <- file.path(out_dir, "tables")
plot_dir <- file.path(out_dir, "plots")
dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

sample_file <- file.path(config_dir, "scRNA_samples.csv")
mapping_file <- file.path(config_dir, "celltype_mapping.csv")

qc_min_features <- 300
qc_max_features <- 7000
qc_max_percent_mt <- 20
n_variable_features <- 2000
harmony_dims <- 1:20
cluster_resolution <- 0.8
milo_prop <- 0.10
cellchat_min_cells <- 10
cellchat_downsample <- 1000

samples <- read.csv(sample_file, stringsAsFactors = FALSE)
required_sample_cols <- c("library_id", "biological_sample", "dataset", "disease", "data_dir")
if (!all(required_sample_cols %in% names(samples))) {
  stop("scRNA_samples.csv must contain: ", paste(required_sample_cols, collapse = ", "))
}

read_10x_expression <- function(path) {
  x <- Read10X(data.dir = path)
  if (is.list(x)) {
    if ("Gene Expression" %in% names(x)) return(x[["Gene Expression"]])
    return(x[[1]])
  }
  x
}

# ------------------------------------------------------------------
# 1. QC, merge, normalization, Harmony integration, clustering
# ------------------------------------------------------------------
seurat_list <- lapply(seq_len(nrow(samples)), function(i) {
  s <- samples[i, ]
  counts <- read_10x_expression(s$data_dir)

  obj <- CreateSeuratObject(
    counts = counts,
    project = s$library_id,
    min.cells = 3,
    min.features = 200
  )

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj$library_id <- s$library_id
  obj$sample_id <- s$biological_sample
  obj$dataset <- s$dataset
  obj$disease <- s$disease

  subset(
    obj,
    subset = nFeature_RNA >= qc_min_features &
      nFeature_RNA <= qc_max_features &
      percent.mt <= qc_max_percent_mt
  )
})

names(seurat_list) <- samples$library_id

seu <- if (length(seurat_list) == 1) {
  seurat_list[[1]]
} else {
  merge(
    seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = names(seurat_list)
  )
}

DefaultAssay(seu) <- "RNA"
try({
  seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
}, silent = TRUE)

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = n_variable_features, verbose = FALSE)
seu <- ScaleData(seu, features = rownames(seu), verbose = FALSE)
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 50, verbose = FALSE)

seu <- RunHarmony(
  seu,
  group.by.vars = "library_id",
  reduction = "pca",
  reduction.save = "harmony",
  plot_convergence = FALSE
)

seu <- FindNeighbors(seu, reduction = "harmony", dims = harmony_dims, k.param = 20, verbose = FALSE)
seu <- FindClusters(seu, resolution = cluster_resolution, algorithm = 1, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = harmony_dims, verbose = FALSE)

saveRDS(seu, file.path(obj_dir, "seurat_harmony_clustered.rds"))

p_integration <- DimPlot(seu, reduction = "umap", group.by = "dataset") +
  ggtitle("Integrated scRNA-seq by dataset") +
  theme_classic()
ggsave(file.path(plot_dir, "UMAP_by_dataset.pdf"), p_integration, width = 6.5, height = 5)

# ------------------------------------------------------------------
# 2. Cluster markers
# ------------------------------------------------------------------
DefaultAssay(seu) <- "RNA"
try({
  seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
}, silent = TRUE)

markers <- FindAllMarkers(
  seu,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = FALSE
)
write.csv(markers, file.path(tab_dir, "cluster_markers.csv"), row.names = FALSE)

top5 <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()
write.csv(top5, file.path(tab_dir, "top5_markers_per_cluster.csv"), row.names = FALSE)

# ------------------------------------------------------------------
# 3. Cell-type annotation
# ------------------------------------------------------------------
mapping <- read.csv(mapping_file, stringsAsFactors = FALSE)
if (!all(c("cluster", "celltype") %in% names(mapping))) {
  stop("celltype_mapping.csv must contain cluster and celltype columns.")
}

seu$celltype <- mapping$celltype[
  match(as.character(seu$seurat_clusters), as.character(mapping$cluster))
]

if (anyNA(seu$celltype)) {
  warning("Some clusters are not represented in celltype_mapping.csv.")
}

saveRDS(seu, file.path(obj_dir, "seurat_annotated.rds"))

p_celltype <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.15
) + ggtitle("Cell types") + theme_classic()

ggsave(file.path(plot_dir, "UMAP_celltypes.pdf"), p_celltype, width = 7, height = 5.5)

composition <- seu@meta.data %>%
  count(disease, celltype, name = "n") %>%
  group_by(disease) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

write.csv(composition, file.path(tab_dir, "celltype_composition_by_disease.csv"), row.names = FALSE)

p_comp <- ggplot(composition, aes(disease, proportion, fill = celltype)) +
  geom_col(width = 0.75) +
  labs(x = NULL, y = "Cell proportion", fill = "Cell type") +
  theme_classic()
ggsave(file.path(plot_dir, "celltype_composition_by_disease.pdf"), p_comp, width = 7, height = 5)

# ------------------------------------------------------------------
# 4. Candidate-gene expression within matched cell types
# ------------------------------------------------------------------
candidate_genes <- c("POSTN", "ITGAM", "GLG1", "RALY", "P4HA1", "CYBB", "PTPA")

deg_by_celltype <- lapply(sort(unique(na.omit(seu$celltype))), function(ct) {
  sub <- subset(seu, subset = celltype == ct)
  if (length(unique(sub$disease)) < 2) return(NULL)

  n_by_group <- table(sub$disease)
  if (!all(c("BM", "NO_BM") %in% names(n_by_group)) || any(n_by_group[c("BM", "NO_BM")] < 3)) {
    return(NULL)
  }

  Idents(sub) <- "disease"
  res <- FindMarkers(
    sub,
    ident.1 = "BM",
    ident.2 = "NO_BM",
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    verbose = FALSE
  )
  res$gene <- rownames(res)
  res$celltype <- ct
  res
}) %>% bind_rows()

if (nrow(deg_by_celltype) > 0) {
  candidate_deg <- deg_by_celltype %>%
    filter(gene %in% candidate_genes)
  write.csv(candidate_deg, file.path(tab_dir, "candidate_genes_BM_vs_nonBM_by_celltype.csv"), row.names = FALSE)
}

# ------------------------------------------------------------------
# 5. Milo differential-abundance analysis
# Biological sample IDs, rather than library IDs, are used as replicates.
# ------------------------------------------------------------------
sce <- as.SingleCellExperiment(seu, assay = "RNA")

harmony_emb <- Embeddings(seu, "harmony")
common_cells <- intersect(colnames(sce), rownames(harmony_emb))
sce <- sce[, common_cells]
harmony_emb <- harmony_emb[common_cells, , drop = FALSE]

d_value <- min(30, ncol(harmony_emb))
reducedDim(sce, "PCA") <- harmony_emb[, seq_len(d_value), drop = FALSE]

if ("umap" %in% Reductions(seu)) {
  umap_emb <- Embeddings(seu, "umap")[common_cells, , drop = FALSE]
  reducedDim(sce, "UMAP") <- umap_emb
}

milo <- Milo(sce)
k_value <- min(30, max(10, floor(ncol(milo) * 0.01)))

milo <- buildGraph(milo, k = k_value, d = d_value, reduced.dim = "PCA")
milo <- makeNhoods(
  milo,
  prop = milo_prop,
  k = k_value,
  d = d_value,
  refined = TRUE,
  reduced_dims = "PCA"
)

milo <- countCells(
  milo,
  meta.data = data.frame(colData(milo)),
  samples = "sample_id"
)

design_df <- data.frame(colData(milo)) %>%
  select(sample_id, disease) %>%
  distinct()
design_df$disease <- factor(design_df$disease, levels = c("NO_BM", "BM"))
rownames(design_df) <- design_df$sample_id
design_df <- design_df[colnames(nhoodCounts(milo)), , drop = FALSE]

milo <- calcNhoodDistance(milo, d = d_value, reduced.dim = "PCA")

da <- testNhoods(
  milo,
  design = ~ disease,
  design.df = design_df,
  fdr.weighting = "graph-overlap"
)

milo <- buildNhoodGraph(milo)
da_annotated <- annotateNhoods(milo, da, coldata_col = "celltype")
da_annotated$celltype_clean <- ifelse(
  da_annotated$celltype_fraction < 0.70,
  "Mixed",
  da_annotated$celltype
)

write.csv(da_annotated, file.path(tab_dir, "milo_differential_abundance.csv"), row.names = FALSE)
write.csv(
  da_annotated %>% filter(SpatialFDR < 0.05),
  file.path(tab_dir, "milo_significant_neighborhoods.csv"),
  row.names = FALSE
)

p_milo <- ggplot(da_annotated, aes(celltype_clean, logFC)) +
  ggbeeswarm::geom_quasirandom(aes(color = SpatialFDR < 0.05), size = 1.2, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(
    x = NULL,
    y = "log fold change (BM vs non-BM)",
    color = "Spatial FDR < 0.05"
  ) +
  theme_classic()
ggsave(file.path(plot_dir, "Milo_differential_abundance.pdf"), p_milo, width = 7, height = 6)

saveRDS(milo, file.path(obj_dir, "milo_object.rds"))

# ------------------------------------------------------------------
# 6. CellChat
# ------------------------------------------------------------------
CellChatDB.use <- CellChatDB.human

get_normalized_matrix <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  try({ obj[["RNA"]] <- JoinLayers(obj[["RNA"]]) }, silent = TRUE)
  NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  GetAssayData(obj, assay = "RNA", layer = "data")
}

run_cellchat_group <- function(seu_obj, disease_use) {
  x <- subset(seu_obj, subset = disease == disease_use & !is.na(celltype))
  Idents(x) <- "celltype"

  counts <- table(Idents(x))
  keep <- names(counts[counts >= cellchat_min_cells])
  x <- subset(x, idents = keep)

  Idents(x) <- "celltype"
  x <- subset(x, downsample = cellchat_downsample)

  data_input <- get_normalized_matrix(x)
  meta <- data.frame(
    celltype = x$celltype,
    disease = x$disease,
    sample_id = x$sample_id,
    row.names = colnames(x)
  )

  cellchat <- createCellChat(data_input[, rownames(meta), drop = FALSE], meta = meta, group.by = "celltype")
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE, type = "truncatedMean", trim = 0.1)
  cellchat <- filterCommunication(cellchat, min.cells = cellchat_min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  cellchat
}

cellchat_non_bm <- run_cellchat_group(seu, "NO_BM")
cellchat_bm <- run_cellchat_group(seu, "BM")

saveRDS(cellchat_non_bm, file.path(obj_dir, "cellchat_nonBM.rds"))
saveRDS(cellchat_bm, file.path(obj_dir, "cellchat_BM.rds"))

comm_non_bm <- subsetCommunication(cellchat_non_bm, thresh = 0.05)
comm_bm <- subsetCommunication(cellchat_bm, thresh = 0.05)

write.csv(comm_non_bm, file.path(tab_dir, "CellChat_nonBM_significant_interactions.csv"), row.names = FALSE)
write.csv(comm_bm, file.path(tab_dir, "CellChat_BM_significant_interactions.csv"), row.names = FALSE)

selected_interactions <- c(
  "SELE_GLG1",
  "POSTN - (ITGAV+ITGB3)",
  "POSTN - (ITGAV+ITGB5)",
  "ICAM1 - (ITGAM+ITGB2)"
)

selected_comm <- comm_bm %>%
  filter(interaction_name %in% selected_interactions)

write.csv(
  selected_comm,
  file.path(tab_dir, "CellChat_BM_selected_interactions.csv"),
  row.names = FALSE
)

build_prob_matrix <- function(comm_df, interaction_use, cell_types) {
  x <- comm_df %>% filter(interaction_name == interaction_use)
  m <- matrix(0, nrow = length(cell_types), ncol = length(cell_types),
              dimnames = list(cell_types, cell_types))
  if (nrow(x) > 0) {
    for (i in seq_len(nrow(x))) {
      src <- as.character(x$source[i])
      tgt <- as.character(x$target[i])
      if (src %in% cell_types && tgt %in% cell_types) m[src, tgt] <- x$prob[i]
    }
  }
  m
}

cell_types <- levels(cellchat_bm@idents)
group_size <- as.numeric(table(cellchat_bm@idents)[cell_types])

pdf(file.path(plot_dir, "CellChat_BM_selected_interactions.pdf"), width = 7, height = 7)
for (interaction_use in selected_interactions) {
  prob_mat <- build_prob_matrix(comm_bm, interaction_use, cell_types)
  if (sum(prob_mat) == 0) next
  netVisual_circle(
    prob_mat,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE,
    edge.width.max = 8,
    title.name = paste0("BM: ", interaction_use)
  )
}
dev.off()

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
