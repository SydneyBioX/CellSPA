#' A wrapper function to calculate all metrics
#'
#'
#' @param spe A Spatial Experiment object
#' @param spe_celltype A string or a vector indicating the cell type information.
#' If it is a string, it will retrive the information from colData.
#' @param sce_ref A SingleCellExperiment object for the reference data
#' @param ref_celltype A vector with length of the sample indicating
#' the cell type label.
#' @param positive_marker_list A list of positive markers
#' @param negative_marker_list A list of negative markers
#' @param nn_celltype_pair A pair of cell type.
#' @param nn_neg_markers_list A list of negative marker corresponding to each cell type in `nn_celltype_pair`.
#' @param exprs_values 	A string indicating which assay of `spe` contains the expression values.
#' @param use_BPPARAM A BiocParallelParam instance determining the parallel
#' back-end to be used during evaluation.
#' @param bin_width A numeric indicates the width of the bin.
#' @param distance_breaks A numeric vector of two or more unique cut points for distance between the pair of cell type.
#' @param verbose A logical value indicates whether to print out running messages.
#'
#' @import SpatialExperiment
#'
#' @return A SpatialExperiment object
#' @export


CellSPA <- function(spe,
                    spe_celltype = NULL,
                    sce_ref,
                    ref_celltype = sce_ref$celltype,
                    positive_marker_list = NULL,
                    negative_marker_list = NULL,
                    nn_celltype_pair = NULL,
                    nn_neg_markers_list = NULL,
                    exprs_values = "logcounts",
                    use_BPPARAM = BiocParallel::SerialParam(),
                    bin_width = 500,
                    distance_breaks = c(seq(0, 50, 10), 100),
                    verbose = TRUE) {

    spe <- generatePolygon(spe, use_BPPARAM = use_BPPARAM)
    spe <- calBaselineAllMetrics(spe,
                                 verbose = verbose,
                                 use_BPPARAM = use_BPPARAM)

    if (verbose) {
        print("Calculating expression correlation metrics")
    }
    # Expression correlation
    spe <- calExpressionCorrelation(spe,
                                    sce_ref,
                                    ref_celltype = ref_celltype,
                                    method = c("pearson", "cosine"),
                                    spe_exprs_values = exprs_values,
                                    ref_exprs_values = "mean")

    spe <- calExpressionCorrelation(spe,
                                    sce_ref,
                                    ref_celltype = ref_celltype,
                                    method = c("pearson", "cosine"),
                                    spe_exprs_values = exprs_values,
                                    ref_exprs_values = "prop_detected")

    if (is.null(positive_marker_list)) {
        positive_marker_list <- generateMarkerList(sce_ref, type = "positive")
    }

    if (is.null(negative_marker_list)) {
        negative_marker_list <- generateMarkerList(sce_ref, type = "negative", t = 1)
    }

    if (is.null(spe_celltype)) {
        spe_celltype = "mean_celltype_correlation"
    }

    if (verbose) {
        print("Calculating marker purity metrics")
    }

    spe <- calMarkerPurity(spe,
                           celltype = spe_celltype,
                           marker_list = positive_marker_list,
                           marker_list_name = "positive",
                           use_BPPARAM = use_BPPARAM,
                           exprs_values = exprs_values)


    spe <- calMarkerPurity(spe,
                           celltype = spe_celltype,
                           marker_list = negative_marker_list,
                           marker_list_name = "negative",
                           use_BPPARAM = use_BPPARAM,
                           exprs_values = exprs_values)



    spe <- calMarkerPct(spe,
                        celltype = spe_celltype,
                        marker_list = positive_marker_list,
                        marker_list_name = "positive",
                        use_BPPARAM = use_BPPARAM,
                        exprs_values = exprs_values)

    spe <- calMarkerPct(spe,
                        celltype = spe_celltype,
                        marker_list = negative_marker_list,
                        marker_list_name = "negative",
                        use_BPPARAM = use_BPPARAM,
                        exprs_values = exprs_values)


    if (verbose) {
        print("Calculating spatial diversity metrics")
    }

    spe <- calSpatialMetricsDiversity(spe,
                                      celltype = spe_celltype,
                                      bin_width = bin_width)

    if (!is.null(nn_celltype_pair) & !is.null(nn_neg_markers_list)) {
        if (verbose) {
            print("Calculating nn marker metrics")
        }

        spe <- calNegMarkerVsDist(spe,
                                  celltype = spe_celltype,
                                  nn_celltype_pair = nn_celltype_pair,
                                  nn_neg_markers_list = nn_neg_markers_list,
                                  exprs_values = exprs_values,
                                  distance_breaks = distance_breaks)
    }



    return(spe)

}


