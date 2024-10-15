#' Calculate marker purity for each cell given a marker list
#'
#' @param spe A Spatial Experiment object
#' @param celltype A string or a vector indicating the cell type information.
#' If it is a string, it will retrive the information from colData.
#' @param marker_list A list of marker, the name of the list corresponding to
#' the cell type information provided in `celltype`.
#' @param marker_list_name A string specifying the name to be used to store the result.
#' @param exprs_values 	A string indicating which assay of `spe` contains the expression values.
#' @param use_BPPARAM A BiocParallelParam instance determining the parallel
#' back-end to be used during evaluation.
#'
#' @import BiocParallel
#' @importFrom SummarizedExperiment assay
#' @importFrom Metrics recall precision
#' @return A SpatialExperiment object
#'
#' @export


calMarkerPurity <- function(spe,
                            celltype,
                            marker_list,
                            marker_list_name = "",
                            exprs_values = "logcounts",
                            use_BPPARAM = BiocParallel::SerialParam()) {

    if ("character" %in% class(celltype) & length(celltype) == 1) {
        celltype <- colData(spe)[, celltype]
    }

    exprsMat <- assay(spe, exprs_values)
    res_mat <- .calculate_marker_purity(exprsMat,
                                        celltype,
                                        marker_list,
                                        use_BPPARAM = use_BPPARAM)

    if (marker_list_name != "") {
        marker_list_name <- paste(marker_list_name,
                                  colnames(res_mat), sep = "_")
    }

    spe@colData[, marker_list_name] <- res_mat
    spe <- .add_metrics(spe, marker_list_name, "cell_level")

    return(spe)

}



.calculate_marker_purity <- function(exprsMat,
                                     celltype,
                                     marker_list,
                                     use_BPPARAM) {

    marker_f1 <- BiocParallel::bplapply(1:length(marker_list), function(i) {
        gene_idx <- rownames(exprsMat) %in%  marker_list[[i]]
        idx <- celltype == names(marker_list)[i]
        res <- apply(exprsMat[, idx, drop = FALSE],  2, function(x) {
            actual <- as.numeric(gene_idx)
            predicted <- as.numeric(x > stats::quantile(x, 1 - mean(gene_idx)))
            precision <- Metrics::precision(actual, predicted)
            recall <- Metrics::recall(actual, predicted)
            f1 <- 2 * precision * recall/(precision + recall)
            f1[is.na(f1)] <- 0
            rbind(f1, precision, recall)
        })
        return(res)
    }, BPPARAM = use_BPPARAM)

    marker_f1 <- do.call(cbind, marker_f1)
    not_calculate_cells <- colnames(exprsMat)[!colnames(exprsMat) %in% colnames(marker_f1)]

    if (length(not_calculate_cells) > 0) {
        marker_f1_no_zero <- matrix(0, nrow = nrow(marker_f1),
                                    ncol = length(not_calculate_cells))
        colnames(marker_f1_no_zero) <- not_calculate_cells
        marker_f1 <- cbind(marker_f1, marker_f1_no_zero)
    }

    marker_f1 <- marker_f1[, colnames(exprsMat)]
    rownames(marker_f1) <- c("F1", "Precision", "Recall")
    marker_f1 <- data.frame(t(marker_f1))
    return(marker_f1)
}



#' Calculate marker expression purity for each cell given a marker list
#'
#' @param spe A Spatial Experiment object
#' @param celltype A string or a vector indicating the cell type information.
#' If it is a string, it will retrive the information from colData.
#' @param marker_list A list of marker, the name of the list corresponding to
#' the cell type information provided in `celltype`.
#' @param marker_list_name A string specifying the name to be used to store the result.
#' @param exprs_values 	A string indicating which assay of `spe` contains the expression values.
#' @param threshold A numeric indicating the threshold to be considered as expressed.
#' @param use_BPPARAM A BiocParallelParam instance determining the parallel
#' back-end to be used during evaluation.
#'
#' @import BiocParallel
#' @importFrom SummarizedExperiment assay
#' @return A SpatialExperiment object
#'
#' @export
calMarkerPct <- function(spe,
                         celltype,
                         marker_list,
                         marker_list_name = "",
                         exprs_values = "logcounts",
                         threshold = 0,
                         use_BPPARAM = BiocParallel::SerialParam()) {

    if ("character" %in% class(celltype) & length(celltype) == 1) {
        celltype <- colData(spe)[, celltype]
    }

    exprsMat <- assay(spe, exprs_values)
    res_mat <- .calculate_marker_exprPct(exprsMat,
                                         celltype,
                                         marker_list,
                                         threshold = threshold,
                                         use_BPPARAM = use_BPPARAM)

    if (marker_list_name != "") {
        marker_list_name <- paste(marker_list_name,
                                  colnames(res_mat), sep = "_")
    }

    spe@colData[, marker_list_name] <- res_mat
    spe <- .add_metrics(spe, marker_list_name, "cell_level")

    return(spe)

}




.calculate_marker_exprPct <- function(exprsMat,
                                      celltype,
                                      marker_list,
                                      threshold = 0,
                                      use_BPPARAM) {

    marker_pct <- BiocParallel::bplapply(1:length(marker_list), function(i) {
        gene_idx <- rownames(exprsMat) %in%  marker_list[[i]]
        idx <- celltype %in% names(marker_list)[i]
        Matrix::colMeans(exprsMat[gene_idx, idx, drop = FALSE] > threshold)
    }, BPPARAM = use_BPPARAM)
    marker_pct <- do.call(c, marker_pct)

    not_calculate_cels <- colnames(exprsMat)[!colnames(exprsMat) %in%
                                                 colnames(marker_pct)]
    if (length(not_calculate_cels) > 0) {
        marker_pct_no_zero <- rep(0, length(not_calculate_cels))
        names(marker_pct_no_zero) <- not_calculate_cels
        marker_pct <- c(marker_pct, marker_pct_no_zero)
    }

    marker_pct <- marker_pct[colnames(exprsMat)]
    marker_pct <- data.frame(exprsPct = marker_pct)
    return(marker_pct)
}








#'  Calculate association with the reference
#'
#'
#'
#' @param spe A SpatialExperiment object
#' @param sce_ref A SingleCellExperiment object for the reference data
#' @param ref_celltype A vector with length of the sample indicating
#' the cell type label.
#' @param method method to compute similarity or distance ("pearson" or
#' "cosine"). Default is set to "pearson".
#' @param spe_exprs_values A string indicating which assay of `spe` contains
#' the expression values. Default is set to "logcounts".
#' @param ref_exprs_values A string indicating which assay of `sce_ref` contains
#' the expression values. Default is set to "mean".
#'
#' @importFrom SummarizedExperiment colData assay
#'
#' @return A SpatialExperiment object
#' @export

calExpressionCorrelation <- function(spe,
                                     sce_ref,
                                     ref_celltype = NULL,
                                     method = c("pearson"),
                                     spe_exprs_values = "logcounts",
                                     ref_exprs_values = "mean") {



    method <- match.arg(method, c("pearson", "cosine"), several.ok = TRUE)

    method <- ifelse(method == "pearson", "correlation", method)

    spe_exprs <- assay(spe, spe_exprs_values)
    ref_exprs <- assay(sce_ref, ref_exprs_values)

    common_genes <- intersect(rownames(ref_exprs), rownames(spe_exprs))

    spe_exprs <- spe_exprs[common_genes, ]
    ref_exprs <- ref_exprs[common_genes, ]

    for (i in method) {
        res <- .cal_expression_correlation(spe_exprs, ref_exprs,
                                           ref_celltype, i)
        colnames(res) <- paste(ref_exprs_values, colnames(res), sep = "_")
        colData(spe)[, colnames(res)] <- res
        spe <- .add_metrics(spe, colnames(res), "cell_level")
    }



    return(spe)
}



#' @importFrom proxyC simil
#' @importFrom RcppParallel setThreadOptions


.cal_expression_correlation <- function(spe_exprs,
                                        ref_exprs,
                                        ref_celltype,
                                        method) {
    RcppParallel::setThreadOptions(1)
    cor_mat <- proxyC::simil(t(ref_exprs),
                             t(spe_exprs),
                             method = method, use_nan = FALSE)
    celltype <- ref_celltype[apply(cor_mat, 2, which.max)]
    cor_res <- apply(cor_mat, 2, max)
    res <- data.frame(cor = cor_res, celltype = celltype)
    colnames(res) <- c(paste("cor", method, sep = "_"),
                       paste("celltype", method, sep = "_"))
    return(res)

}







#'  Calculate association with the reference for the aggregated profile
#'
#'
#'
#' @param spe A SpatialExperiment object
#' @param celltype A string or a vector indicating the cell type information.
#' @param sce_ref A SingleCellExperiment object for the reference data
#' @param ref_celltype A vector with length of the sample indicating
#' the cell type label.
#' @param method method to compute similarity or distance ("pearson" or
#' "cosine"). Default is set to "pearson".
#' @param spe_exprs_values A string indicating which assay of `spe` contains
#' the expression values. Default is set to "logcounts".
#' @param ref_exprs_values A string indicating which assay of `sce_ref` contains
#' the expression values. Default is set to "mean".
#'
#' @importFrom SummarizedExperiment colData assay
#' @importFrom proxyC simil
#'
#' @return A SpatialExperiment object
#' @export

calAggExpressionCorrelation <- function(spe,
                                        celltype,
                                        sce_ref,
                                        ref_celltype = NULL,
                                        method = c("pearson"),
                                        spe_exprs_values = "logcounts",
                                        ref_exprs_values = "mean") {


    if ("character" %in% class(celltype) & length(celltype) == 1) {
        celltype <- colData(spe)[, celltype]
    }


    if ("character" %in% class(ref_celltype) & length(ref_celltype) == 1) {
        ref_celltype <- colData(sce_ref)[, ref_celltype]
    }

    method <- match.arg(method, c("pearson", "cosine"), several.ok = TRUE)

    method <- ifelse(method == "pearson", "correlation", method)

    spe_exprs <- assay(spe, spe_exprs_values)

    spe_agg_exprs <- sapply(names(table(ref_celltype)), function(x) {
        Matrix::rowMeans(spe_exprs[, celltype == x, drop = FALSE] != 0)
    })


    ref_exprs <- assay(sce_ref, ref_exprs_values)
    common_genes <- intersect(rownames(ref_exprs), rownames(spe_exprs))

    spe_agg_exprs <- spe_agg_exprs[common_genes, ]
    spe_agg_exprs[is.na(spe_agg_exprs)] <- 0
    ref_exprs <- ref_exprs[common_genes, ]


    for (i in method) {
        res <- proxyC::simil(as(t(ref_exprs), "CsparseMatrix"),
                             as(t(spe_agg_exprs), "CsparseMatrix"),
                             method = i, use_nan = FALSE)

        metrics_name <- paste('agg', ref_exprs_values, i, sep = "_")
        spe@metadata$CellSPA$similarity_metrics[[metrics_name]] <- as.matrix(res)
        #spe <- .add_metrics(spe, colnames(res), "celltype_level")

    }



    return(spe)
}




