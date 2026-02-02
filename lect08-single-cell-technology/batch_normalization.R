## Batch normalization helper functions and code
## Used by lect08-single-cell-technology.Rmd

library(Seurat)

#' Load multiple PDAC samples
#' @param samples Vector of sample names
#' @param data_dir Path to data directory
load_pdac_samples <- function(samples, data_dir) {
    seu_list <- lapply(samples, function(s) {
        Read10X(file.path(data_dir, s)) |>
            CreateSeuratObject(project = s, min.cells = 3, min.features = 200)
    })
    names(seu_list) <- samples
    seu_list
}

#' Preprocess each sample (normalize and find variable features)
#' @param seu_list List of Seurat objects
preprocess_samples <- function(seu_list) {
    lapply(seu_list, function(seu) {
        seu |>
            NormalizeData() |>
            FindVariableFeatures(nfeatures = 2000)
    })
}

#' Merge samples without integration
#' @param seu_list List of Seurat objects
merge_without_integration <- function(seu_list) {
    seu_merged <- merge(seu_list[[1]], seu_list[-1])
    seu_merged |>
        FindVariableFeatures(nfeatures = 2000) |>
        ScaleData() |>
        RunPCA(npcs = 30) |>
        RunUMAP(dims = 1:20)
}

#' Integrate samples using Seurat CCA
#' @param seu_list List of preprocessed Seurat objects
integrate_seurat_cca <- function(seu_list) {
    features <- SelectIntegrationFeatures(seu_list)
    anchors <- FindIntegrationAnchors(seu_list, anchor.features = features)
    seu_integrated <- IntegrateData(anchors)

    seu_integrated |>
        ScaleData() |>
        RunPCA(npcs = 30) |>
        RunUMAP(dims = 1:20)
}

#' Integrate samples using Harmony
#' @param seu_merged Merged Seurat object (not integrated)
integrate_harmony <- function(seu_merged) {
    if (!requireNamespace("harmony", quietly = TRUE)) {
        stop("harmony package required")
    }
    seu_merged |>
        harmony::RunHarmony("orig.ident") |>
        RunUMAP(reduction = "harmony", dims = 1:20)
}
