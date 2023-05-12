

#' Calculate the spatial diversity of base line metrics
#'
#'
#' @param spe A Spatial Experiment object
#' @param celltype A string or a vector indicating the cell type information.
#' If it is a string, it will retrive the information from colData.
#' @param bin_width A numeric indicates the width of the bin.
#'
#' @import SpatialExperiment
#' @importFrom stats sd cor
#'
#' @return A SpatialExperiment object
#' @export


calSpatialMetricsDiversity <- function(spe,
                                       celltype,
                                       bin_width = 500) {


    if ("character" %in% class(celltype) & length(celltype) == 1) {
        celltype <- colData(spe)[, celltype]
    }


    # Create spatial bins

    coord <- SpatialExperiment::spatialCoords(spe)
    coord_x <- coord[, 1]
    coord_y <- coord[, 2]
    num_bin_x <- ceiling(diff(range(coord_x))/bin_width)
    x_bin_groups <- cut(coord_x, breaks = max(num_bin_x, 2))
    num_bin_y <- ceiling(diff(range(coord_y))/bin_width)
    y_bin_groups <- cut(coord_y, breaks = max(num_bin_y, 2))


    bin_groups <- paste(x_bin_groups, y_bin_groups, sep = '|')
    df_bin <-  data.frame(x_bin_groups, y_bin_groups, bin_groups)
    colData(spe)[, colnames(df_bin)] <- data.frame(x_bin_groups, y_bin_groups, bin_groups)

    res <- table(celltype,
                 bin_groups)

    # apply(res, 2, function(x) chisq.test(x)$residuals)
    bin_entropy <- apply(res, 2, function(x) entropy::entropy(x))

    df_entropy_bin <- data.frame(bin = names(bin_entropy),
                                 entropy = bin_entropy)
    df_entropy_bin <- cbind(df_entropy_bin,
                            do.call(rbind, strsplit(as.character(df_entropy_bin$bin), "\\|")))
    colnames(df_entropy_bin) <- c("bin", "entropy",  "x_bin", "y_bin")
    df_entropy_bin$x_bin <- factor(df_entropy_bin$x_bin, levels = levels(x_bin_groups))
    df_entropy_bin$y_bin <- factor(df_entropy_bin$y_bin, levels = levels(y_bin_groups))

    # Combine with cell type proportion
    celltype_prop <- apply(res, 2, function(x) x/sum(x))
    celltype_prop <- t(celltype_prop)
    colnames(celltype_prop) <- paste("cellTypeProp", colnames(celltype_prop), sep = "_")
    df_entropy_bin <- cbind(df_entropy_bin, celltype_prop[as.character(df_entropy_bin$bin), ])


    common_cell_metrics <- intersect(.CellSPARenvir[["cell_level_baseline_metrics"]],
                                     colnames(colData(spe)))


    for (i in common_cell_metrics) {
        if (i %in% c("total_transciprts", "total_genes", "cell_area")) {
            df_bin_stats <- stats::aggregate(log10(colData(spe)[, i]),
                                             list(bin_groups), function(x) sd(x)/mean(x))
            colnames(df_bin_stats) <- c("bin", paste("cv", i, sep = "_"))
            df_entropy_bin <- merge(df_entropy_bin, df_bin_stats, by = "bin")

            df_bin_stats <- stats::aggregate((colData(spe)[, i]),
                                             list(bin_groups), function(x) sd(x))
            colnames(df_bin_stats) <- c("bin", paste("sd", i, sep = "_"))
            df_entropy_bin <- merge(df_entropy_bin, df_bin_stats, by = "bin")

        } else {
            df_bin_stats <- stats::aggregate((colData(spe)[, i]),
                                             list(bin_groups), function(x) sd(x)/mean(x))
            colnames(df_bin_stats) <- c("bin", paste("cv", i, sep = "_"))
            df_entropy_bin <- merge(df_entropy_bin, df_bin_stats, by = "bin")


            df_bin_stats <- stats::aggregate((colData(spe)[, i]),
                                             list(bin_groups), function(x) sd(x))
            colnames(df_bin_stats) <- c("bin", paste("sd", i, sep = "_"))
            df_entropy_bin <- merge(df_entropy_bin, df_bin_stats, by = "bin")
        }


    }

    cor_cv <- cor(df_entropy_bin[, "entropy"],
                  df_entropy_bin[, paste("cv", common_cell_metrics, sep = "_")],
                  use = "complete.obs")[1, ]
    cor_sd <- cor(df_entropy_bin[, "entropy"],
                  df_entropy_bin[, paste("sd", common_cell_metrics, sep = "_")],
                  use = "complete.obs")[1, ]

    cor_res <- c(cor_cv, cor_sd)

    spe@metadata$CellSPA$spatialMetricsDiversity <- list(statistics = cor_res,
                                                         results = df_entropy_bin)
    spe <- .add_metrics(spe,
                        names(spe@metadata$CellSPA$spatialMetricsDiversity$statistics),
                        "bin_level")

    return(spe)


}






