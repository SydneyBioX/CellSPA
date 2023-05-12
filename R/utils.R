
#' Subset Spatial experiment object by column
#'
#'
#' @param spe A Spatial Experiment object
#' @param col_idx The column index to keep.
#'
#' @return A SpatialExperiment object
#' @export
subset <- function(spe, col_idx) {
    spe <- spe[, col_idx]
    if (!is.null(spe@metadata$CellSegOutput)) {
        keep <- spe@metadata$CellSegOutput$cell_id %in% spe$cell_id
        spe@metadata$CellSegOutput <- spe@metadata$CellSegOutput[keep, ]
    }

    return(spe)
}

.add_metrics <- function(spe, metrics_name, type) {
    spe@metadata$CellSPA$metrics[[type]] <- sort(unique(c(spe@metadata$CellSPA$metrics[[type]],
                                                          metrics_name)))
    return(spe)
}


.add_dataset_metrics <- function(spe) {
    spe@metadata$CellSPA$dataset_level_metrics <- list(num_cells = ncol(spe),
                                                       num_genes = nrow(spe))
    spe <- .add_metrics(spe, c("num_cells", "num_genes"), "dataset_level")
    return(spe)
}


.initialise_CellSPA_list <- function(spe) {
    spe@metadata$CellSPA <- list()
    spe@metadata$CellSPA$metrics <- vector("list", 4)
    names(spe@metadata$CellSPA$metrics) <- c("dataset_level",
                                             "gene_level",
                                             "cell_level",
                                             "bin_level")
    return(spe)
}



