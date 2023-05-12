
#' Calculate the baseline metrics of segmented cell from cell segmentation output
#'
#'
#' @param spe A Spatial Experiment object
#' @param metrics A character vector indicating the metrics to be calculated
#' @param rerun A logical value indicates whether to recalculate existing metrics.
#' @param verbose A logical value indicates whether to print out running messages.
#'
#' @importFrom SingleCellExperiment counts
#' @return A SpatialExperiment object
#' @export

calBaselineCellMetrics <- function(spe,
                                   metrics = c("total_transciprts", "total_genes",
                                               "total_cells", "meanExprsPct_cells"),
                                   rerun = TRUE,
                                   verbose = TRUE) {

    metrics <- match.arg(metrics, c("total_transciprts", "total_genes",
                                    "total_cells", "meanExprsPct_cells"), several.ok = TRUE)
    if (!rerun) {
        metrics <- setdiff(metrics, unlist(spe@metadata$CellSPA$metrics))
    }

    if (verbose) {
        print(paste("Metrics to run: ", paste(metrics, collapse = ", ")))
    }

    spe <- .cal_baseline(spe, metrics = metrics)
    return(spe)
}

#' Calculate the baseline metrics of cell shapes from cell segmentation output
#'
#'
#' @param spe A Spatial Experiment object
#' @param metrics A character vector indicating the metrics to be calculated
#' @param use_BPPARAM a BiocParallelParam instance determining the parallel
#' back-end to be used during evaluation.
#' @param rerun A logical value indicates whether to recalculate existing metrics.
#' @param verbose A logical value indicates whether to print out running messages.
#'
#' @import BiocParallel

#' @return A SpatialExperiment object
#' @export

calBaselineCellShapeMetrics <- function(spe,
                                        metrics = c("cell_area", "elongation",
                                                    "compactness", "eccentricity",
                                                    "sphericity", "solidity",
                                                    "convexity", "circularity"),
                                        use_BPPARAM = BiocParallel::SerialParam(),
                                        rerun = TRUE,
                                        verbose = TRUE) {

    metrics <- match.arg(metrics, c("cell_area", "elongation",
                                    "compactness", "eccentricity",
                                    "sphericity", "solidity",
                                    "convexity", "circularity"), several.ok = TRUE)

    if (!rerun) {
        metrics <- setdiff(metrics, unlist(spe@metadata$CellSPA$metrics))
        if (verbose) {
            print(paste("Metrics to run: ", paste(metrics, collapse = ", ")))
        }

    }

    if (verbose) {
        BiocParallel::bpprogressbar(use_BPPARAM) <- TRUE
    }

    spe <- .cal_baseline_tiff(spe, metrics = metrics, use_BPPARAM = use_BPPARAM,
                              verbose = verbose)

    return(spe)
}




#' Calculate the baseline metrics of cell from cell segmentation output
#'
#'
#' @param spe A Spatial Experiment object
#' @param metrics A character vector indicating the metrics to be calculated
#' @param use_BPPARAM a BiocParallelParam instance determining the parallel
#' back-end to be used during evaluation.
#' @param rerun A logical value indicates whether to recalculate existing metrics.
#' @param verbose A logical value indicates whether to print out running messages.
#'
#' @import BiocParallel

#' @return A SpatialExperiment object
#' @export

calBaselineAllMetrics <- function(spe,
                                  metrics = c("total_transciprts", "total_genes",
                                              "total_cells", "meanExprsPct_cells",
                                              "cell_area", "elongation",
                                              "compactness", "eccentricity",
                                              "sphericity", "solidity",
                                              "convexity", "circularity",
                                              "density"),
                                  use_BPPARAM = BiocParallel::SerialParam(),
                                  rerun = TRUE,
                                  verbose = TRUE) {

    non_cell_shape_metrics <- intersect(metrics,
                                        .CellSPARenvir[["non_cell_shape_baseline_metrics"]])


    if (length(non_cell_shape_metrics) != 0) {
        spe <- calBaselineCellMetrics(spe,
                                      metrics = non_cell_shape_metrics,
                                      rerun = rerun,
                                      verbose = verbose)
    }

    cell_shape_metrics <- intersect(metrics,
                                    .CellSPARenvir[["cell_shape_baseline_metrics"]])

    if (length(non_cell_shape_metrics) != 0 & !is.null(spe@metadata$CellSegOutput)) {
        spe <- calBaselineCellShapeMetrics(spe,
                                           metrics = cell_shape_metrics,
                                           use_BPPARAM = use_BPPARAM,
                                           rerun = rerun,
                                           verbose = verbose)
    }

    if ("density" %in% metrics & !is.null(spe@metadata$CellSegOutput)) {
        spe <- density(spe)
    }

    return(spe)
}



#' @importFrom SummarizedExperiment colData rowData
#' @importFrom Matrix rowMeans rowSums colSums

.cal_baseline <- function(spe, metrics = c("total_transciprts", "total_genes",
                                           "total_cells", "meanExprsPct_cells")) {
    if (!"CellSPA" %in% names(spe@metadata)) {
        spe@metadata$CellSPA <- list()
    }


    if ("total_transciprts" %in% metrics) {

        total_transciprts <- Matrix::colSums(counts(spe))
        colData(spe)$total_transciprts <- total_transciprts
        spe@metadata$CellSPA$metrics <- c(spe@metadata$CellSPA$metrics,
                                          "total_transciprts")
        spe <- .add_metrics(spe, "total_transciprts", "cell_level")
    }

    if ("total_genes" %in% metrics) {
        total_genes <- Matrix::colSums(counts(spe) != 0)
        colData(spe)$total_genes <- total_genes
        spe <- .add_metrics(spe, "total_genes", "cell_level")
    }

    if ("total_cells" %in% metrics) {
        total_cells <- Matrix::rowSums(counts(spe) != 0)
        SummarizedExperiment::rowData(spe)$total_cells <- total_cells
        spe <- .add_metrics(spe, "total_cells", "gene_level")
    }

    if ("meanExprsPct_cells" %in% metrics) {
        meanExprsPct_cells <- Matrix::rowMeans(counts(spe) != 0)
        SummarizedExperiment::rowData(spe)$meanExprsPct_cells <- meanExprsPct_cells
        spe <- .add_metrics(spe, "meanExprsPct_cells", "gene_level")
    }



    return(spe)
}



#' Generating polygon
#'
#'
#' @param spe A Spatial Experiment object
#' @param use_BPPARAM a BiocParallelParam instance determining the parallel
#' back-end to be used during evaluation.
#' @param rerun A logical value indicates whether to recalculate existing metrics.
#' @param verbose A logical value indicates whether to print out running messages.
#'
#' @import BiocParallel

#' @return A SpatialExperiment object
#'
#' @export

generatePolygon <- function(spe,
                            use_BPPARAM = BiocParallel::SerialParam(),
                            rerun = TRUE,
                            verbose = TRUE) {


    tiff_res <- spe@metadata$CellSegOutput
    tiff_res_list <- split(tiff_res[, c(1:2)], tiff_res$cell_id)

    if (verbose) {
        BiocParallel::bpprogressbar(use_BPPARAM) <- TRUE
    }


    if (is.null(spe@metadata$CellSPA$poly_objects) | rerun) {
        spe@metadata$CellSPA$poly_objects <- BiocParallel::bplapply(tiff_res_list, function(x) {
            get_ashape_chull_poly(x)
        }, BPPARAM = use_BPPARAM)
    }


    return(spe)
}




#' @import BiocParallel

.cal_baseline_tiff <- function(spe,
                               metrics = c("cell_area", "elongation",
                                           "compactness", "eccentricity",
                                           "sphericity", "solidity",
                                           "convexity", "circularity"),
                               use_BPPARAM = SerialParam(),
                               verbose = verbose) {
    tiff_res <- spe@metadata$CellSegOutput
    tiff_res_list <- split(tiff_res[, c(1:2)], tiff_res$cell_id)
    # filter <- unlist(lapply(tiff_res_list, nrow)) <= 3
    # tiff_res_list <- tiff_res_list[!filter]

    if ("cell_area" %in% metrics) {
        if (!is.null(tiff_res)) {
            cell_area <- table(tiff_res$cell_id)
            cell_area <- as.numeric(cell_area[as.character(spe$cell_id)])
            colData(spe)$cell_area <- cell_area
            spe <- .add_metrics(spe, "cell_area", "cell_level")
        }
    }


    if (is.null(spe@metadata$CellSPA$poly_objects)) {
        if (verbose) {
            print("Generating ashape and chull...")
        }
        spe@metadata$CellSPA$poly_objects <- BiocParallel::bplapply(tiff_res_list, function(x) {
            get_ashape_chull_poly(x)
        }, BPPARAM = use_BPPARAM)

    }

    poly_objects <- spe@metadata$CellSPA$poly_objects

    for (i in setdiff(metrics, "cell_area")) {
        spe <- .run_cell_shape_metrics(spe,
                                       poly_objects,
                                       method = i,
                                       use_BPPARAM = use_BPPARAM,
                                       verbose = verbose)
    }

    return(spe)
}


density <- function(spe) {
    colData(spe)$density <- spe$total_transciprts/spe$cell_area
    spe <- .add_metrics(spe, "density", "cell_level")
    return(spe)
}
