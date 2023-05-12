

#' Calculate the spatial diversity of base line metrics
#'
#'
#' @param spe A Spatial Experiment object
#' @param celltype A string or a vector indicating the cell type information.
#' If it is a string, it will retrive the information from colData.
#' @param nn_celltype_pair A pair of cell type.
#' @param nn_neg_markers_list A list of negative marker corresponding to each cell type in `nn_celltype_pair`.
#' @param exprs_threshold A numeric indicating the threshold to be considered as expressed.
#' @param distance_breaks A numeric vector of two or more unique cut points for distance between the pair of cell type.
#' @param exprs_values 	A string indicating which assay of `spe` contains the expression values.
#'
#' @import SpatialExperiment
#' @importFrom stats sd cor aggregate
#' @importFrom proxyC dist
#'
#' @return A SpatialExperiment object
#'
#' @export


calNegMarkerVsDist <- function(spe,
                               celltype,
                               nn_celltype_pair,
                               nn_neg_markers_list,
                               exprs_threshold = 0,
                               distance_breaks = c(seq(0, 50, 10), 100),
                               exprs_values = "logcounts") {



    if ("character" %in% class(celltype) & length(celltype) == 1) {
        celltype <- colData(spe)[, celltype]
    }

    nn_neg_markers_list <- lapply(nn_neg_markers_list, function(x) intersect(x, rownames(spe)))


    coord <- SpatialExperiment::spatialCoords(spe)


    exprsMat <- assay(spe, exprs_values)

    neg_res_list <- lapply(names(nn_neg_markers_list), function(l) {
        # Calculate the distance for a pair of cell type
        dist_mat <- proxyC::dist(as(coord[grepl(l, celltype), , drop = FALSE], "CsparseMatrix"),
                                 as(coord[grepl(setdiff(nn_celltype_pair, l), celltype), , drop = FALSE], "CsparseMatrix"))
        dist_to_nn_celltype <- apply(dist_mat, 1, min)
        neg_mat <- exprsMat[nn_neg_markers_list[[l]],
                            names(dist_to_nn_celltype),
                            drop = FALSE]
        distance_bin <- cut(dist_to_nn_celltype,
                            breaks = c(distance_breaks, max(dist_to_nn_celltype)),
                            include.lowest = TRUE)


        neg_res <- apply(neg_mat, 1, function(m) {
            stats::aggregate(m, list(distance_bin), function(exprs) {
                mean(exprs > exprs_threshold)
            }, drop = FALSE)}$x)

        rownames(neg_res) <- names(table(distance_bin))


        return(neg_res)
    })

    names(neg_res_list) <- names(nn_neg_markers_list)

    spe@metadata$CellSPA$`negMarkerExprs_vs_dist` <- neg_res_list

    return(spe)
}


