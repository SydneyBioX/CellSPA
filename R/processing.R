#' Preprocessing the the spatial experiment data.
#'
#'
#' @param spe A Spatial Experiment object
#' @param qc_range A list vector with two elements indicating
#' the range of the total transciprts of the cells to keep.
#' @param verbose A logical value indicates whether to print out running messages.
#' @param seed A numeric indicates the seed to be set
#' @importFrom SummarizedExperiment colData
#' @return A SpatialExperiment object
#' @export


processingSPE <- function(spe,
                          qc_range = list(total_transciprts = c(20, 2000),
                                          total_genes = c(20, Inf)),
                          verbose = TRUE,
                          seed = 2023) {

    if (!all(names(qc_range) %in% colnames(colData(spe)))) {
        stop("One of the names of qc_range list is not in colData(spe).")
    }

    for (i in names(qc_range)) {
        spe <- spe[, colData(spe)[, i] >= qc_range[[i]][1] &
                       colData(spe)[, i] <= qc_range[[i]][2]]

    }

    tif_output <- spe@metadata$CellSegOutput
    spe@metadata$CellSegOutput <- tif_output[tif_output$cell_id %in% spe$cell_id, ]


    spe <- scater::logNormCounts(spe)
    set.seed(seed)
    spe <- scater::runPCA(spe)
    spe <- scater::runUMAP(spe,
                           dimred = "PCA",
                           min_dist = 0.3,
                           verbose = verbose)
    return(spe)
}



#' Preprocessing the the reference scRNA-seq data
#'
#'
#' @param sce_ref A SingleCellExperiment object
#' @param celltype A vector indicating cell type information.
#' @param subset_row A vector specifying the subset of features to keep for the
#' reference data. This can be a character vector of row names,
#' an integer vector of row indices or a logical vector.
#'
#' @importFrom scater aggregateAcrossCells
#' @import SummarizedExperiment
#' @import SingleCellExperiment
#'
#' @return A SingleCellExperiment object
#'
#' @export



processingRef <- function(sce_ref, celltype, subset_row = NULL) {

    if (!is.null(subset_row)) {
        if (is.character(subset_row)) {
            subset_row <- intersect(subset_row, rownames(sce_ref))
        }
        sce_ref <- sce_ref[subset_row, ]
    }

    ref_mean <- scater::aggregateAcrossCells(sce_ref,
                                             celltype,
                                             use.assay.type = "logcounts",
                                             statistics = "mean")

    ref_prop_detected <- scater::aggregateAcrossCells(sce_ref,
                                                      celltype,
                                                      use.assay.type = "logcounts",
                                                      statistics = "prop.detected")

    ref_celltype_prop <- table(celltype)/length(celltype)

    sce_ref_agg <- SingleCellExperiment(assay = list(mean = assay(ref_mean),
                                                     prop_detected = assay(ref_prop_detected)),
                                        colData = unlist(table(celltype)/length(celltype)))

    rowData(sce_ref_agg) <- rowData(ref_mean)


    return(sce_ref_agg)
}






#' Generating marker list using the reference data
#'
#' @param sce_ref A SingleCellExperiment object
#' @param exprs_values A string indicating which assay of `sce_ref` contains
#' the expression values. Default is set to "mean".
#' @param q A numeric indicates the quantile cutoff for selecting the marker list.
#' By default is 0.9.
#' @param t A numeric indicates the cutoff of that a marker allows to appear
#' in maximum proportion of cell types. By default is 0.25.
#' @param verbose A logical value indicates whether to print out running messages.
#' @param type A string indicates whether to calculate positive or negative markers.
#'
#' @import SummarizedExperiment
#' @import SingleCellExperiment
#'
#' @return A list object
#'
#' @export




generateMarkerList <- function(sce_ref,
                               exprs_values = "mean",
                               q = 0.9,
                               t = 0.25,
                               verbose = TRUE,
                               type = "positive") {


    type <- match.arg(type, c("positive", "negative"))

    ref_exprs <- assay(sce_ref, exprs_values)

    if (type == "positive") {
        gene_marker_list <- lapply(1:ncol(ref_exprs), function(i) {
            w <- (ref_exprs[, i] - rowMeans(ref_exprs[, -i, drop = FALSE]))
            names(which(w > stats::quantile(w, q)))
        })
    }

    if (type == "negative") {
        gene_marker_list <- lapply(1:ncol(ref_exprs), function(i) {
            w <- (ref_exprs[, i] - rowMeans(ref_exprs[, -i, drop = FALSE]))
            names(which(w < stats::quantile(w, 1 - q)))
        })
    }

    names(gene_marker_list) <- colnames(ref_exprs)

    if (verbose) {
        print("For the top 10% of genes, the overlap freq between cell types")
        print(table(table(unlist(gene_marker_list))))
    }


    remove_genes <- names(which(table(unlist(gene_marker_list)) >=  t * ncol(ref_exprs)))
    gene_marker_list <- lapply(gene_marker_list, function(x) x[!x %in% remove_genes])



    if (verbose) {
        print("Length of the positive gene marker")
        print(unlist(lapply(gene_marker_list, length)))
    }


    return(gene_marker_list)

}
