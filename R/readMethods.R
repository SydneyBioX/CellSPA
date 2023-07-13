
#' Read Xenium output
#'
#'
#' @param data_dir a character specifying a single output directory as
#' returned by Space Ranger, corresponding to one 10x Genomics Xenium sample.
#' @param keep_type a character vector indicating the type of features to be
#' included in the object. By default it is set to `"Gene Expression"`.
#' @param feature_named_by a character indicating the way that the features are named.
#' By default it is set to `"Symbol"`.
#' @param tiff_path a character specifying the path of a tiff file output from BIDCell.
#' @param method_name a character specifying the name of the method.
#'
#' @import SpatialExperiment
#' @importFrom utils read.csv
#' @importFrom DropletUtils read10xCounts
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment rowData
#'
#' @return A SpatialExperiment object
#' @examples
#'
#' tenX_output_dir <- system.file("extdata/10x_output_subset", package = "CellSPA")
#' tenX_output_tif <- system.file("extdata/10x_output_subset/10x_from_csv_subset.tif",
#'                                 package = "CellSPA")
#' spe_10x <- readXenium(tenX_output_dir,
#'                       tiff_path = tenX_output_tif)
#'
#' @export


readXenium <- function(data_dir,
                       keep_type = "Gene Expression",
                       feature_named_by = "Symbol",
                       tiff_path = NULL,
                       method_name = "10x") {

    # Read count matrix
    sce <- DropletUtils::read10xCounts(file.path(data_dir, "cell_feature_matrix"))
    cell_meta <- read.csv(file.path(data_dir, "cells.csv.gz"))
    colnames(sce) <- paste("Cell", cell_meta$cell_id, sep = "_")

    spe <- SpatialExperiment::SpatialExperiment(assay = list(counts = counts(sce)),
                                                colData =  DataFrame(cell_meta),
                                                rowData = rowData(sce),
                                                spatialCoordsNames = c("x_centroid",
                                                                       "y_centroid"))

    if (!is.null(keep_type)) {
        keep_genes <- rowData(spe)$Type %in% keep_type
        spe <- spe[keep_genes, ]
    }

    if (feature_named_by %in% colnames(rowData(spe))) {
        rownames(spe) <- rowData(spe)[, feature_named_by]
    }




    spe <- .initialise_CellSPA_list(spe)
    spe <- .add_dataset_metrics(spe)

    if (!is.null(tiff_path)) {
        tif_output <- readTiffOutput(tiff_path)
        tif_output <- tif_output[tif_output$cell_id %in% spe$cell_id, ]
        spe@metadata$CellSegOutput <- tif_output

        if (!all(spe$cell_id %in% spe@metadata$CellSegOutput$cell_id)) {
            message("There are cells that are not in the tiff file,
                    please check your tiff file")
            spe <- spe[, spe$cell_id %in% spe@metadata$CellSegOutput$cell_id]
        }
    }

    spe <- .cal_baseline(spe)
    spe@metadata$CellSPA$method <- method_name

    return(spe)
}


#' Read BIDCell output
#'
#'
#' @param data_dir a character specifying a single output directory as
#' returned by BIDCell, corresponding to one sample.
#' @param meta_idx a numeric vector indicating the indices of the column that
#' corresponding to the meta data in the output files.
#' @param tiff_path a character specifying the path of a tiff file output from BIDCell.
#' @param method_name a character specifying the name of the method.
#'
#' @import SpatialExperiment
#' @importFrom methods as
#' @importFrom Matrix t
#' @return A SpatialExperiment object
#' @examples
#'
#' data_dir <- system.file("extdata/BIDCell_csv_output", package = "CellSPA")
#' tiff_path <- system.file("extdata/BIDCell_output_subset.tif", package = "CellSPA")
#' spe <- readBIDCell(data_dir,
#'                    tiff_path = tiff_path,
#'                    method_name = "BIDCell")
#'
#' @export


readBIDCell <- function(data_dir,
                        meta_idx = NULL,
                        tiff_path = NULL,
                        method_name = NULL) {

    all_files <- list.files(data_dir)

    cell_outputs <- lapply(1:length(all_files), function(i) {
        res <- read.csv(file.path(data_dir, all_files[i]))
        res$slide <- i
        res
    })





    cell_outputs <- do.call(rbind, cell_outputs)

    if (is.null(meta_idx)) {
        meta_idx <- c(1:8, (ncol(cell_outputs)-2):ncol(cell_outputs))
    }


    data <- cell_outputs[, -meta_idx]
    meta <- cell_outputs[, meta_idx]


    duplicated_cell_id <- duplicated(cell_outputs$cell_id)



    if (sum(duplicated_cell_id) > 0) {
        warning(paste("There are", sum(duplicated_cell_id),
                      "cells with duplicated cell id"))
        data <- aggregate.Matrix(data, cell_outputs$cell_id)
        rownames(data) <- paste("Cell", rownames(data), sep = "_")
        meta <- meta[!duplicated_cell_id, ]
        rownames(meta) <- paste("Cell", meta$cell_id, sep = "_")
        meta <- meta[rownames(data), ]
    } else {
        rownames(data) <- paste("Cell", cell_outputs$cell_id, sep = "_")
    }



    data <- as(as.matrix(data), "CsparseMatrix")



    spe <- SpatialExperiment(assay = list(counts = Matrix::t(data)),
                             colData = data.frame(meta),
                             spatialCoordsNames = c("cell_centroid_x",
                                                    "cell_centroid_y"))

    spe <- .initialise_CellSPA_list(spe)
    spe <- .add_dataset_metrics(spe)

    if (!is.null(tiff_path)) {
        tif_output <- readTiffOutput(tiff_path)
        tif_output <- tif_output[tif_output$cell_id %in% spe$cell_id, ]
        spe@metadata$CellSegOutput <- tif_output

        if (!all(spe$cell_id %in% spe@metadata$CellSegOutput$cell_id)) {
            warning("There are cells that are not in the tiff file,
                    please check your tiff file")
            spe <- spe[, spe$cell_id %in% spe@metadata$CellSegOutput$cell_id]
        }
    }

    spe <- .cal_baseline(spe)
    spe@metadata$CellSPA$method <- method_name

    return(spe)
}



#' Read Tiff output from BIDCell
#'
#'
#' @param tiff_path a character specifying the path of a tiff file.
#'
#'
#' @importFrom tiff readTIFF
#' @importFrom reshape2 melt
#' @return A SpatialExperiment object
#'
#' tiff_path <- system.file("extdata/BIDCell_output_subset.tif", package = "CellSPA")
#' tiff_res <- readTiffOutput(tiff_path = tiff_path)
#' @export

readTiffOutput <- function(tiff_path) {

    tiff_res <- readTIFF(tiff_path, as.is = TRUE)
    tiff_res <- reshape2::melt(tiff_res)
    tiff_res <- tiff_res[tiff_res$value != 0, ]

    colnames(tiff_res) <- c("coord_y", "coord_x", "cell_id")
    tiff_res <- tiff_res[, c("coord_x", "coord_y", "cell_id")]
    tiff_res$cell_id <- as.numeric(tiff_res$cell_id)
    tiff_res$coord_x <- as.numeric(tiff_res$coord_x)
    tiff_res$coord_y <- as.numeric(tiff_res$coord_y)
    tiff_res <- tiff_res[order(tiff_res$cell_id), ]
    rownames(tiff_res) <- NULL
    return(tiff_res)
}




#' Read Vizgen output
#'
#'
#' @param data_dir a character specifying a single output directory from Vizgen.
#' The data_dir should include one cell_by_gene.csv and one cell_metadata.csv.
#' @param filter_gene_pattern Pattern of the genes to be filtered out.
#' @param tiff_path a character specifying the path of a tiff file output from BIDCell.
#' @param method_name a character specifying the name of the method.
#'
#' @import SpatialExperiment
#' @importFrom data.table fread
#' @importFrom utils read.csv
#' @importFrom DropletUtils read10xCounts
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment rowData
#'
#' @return A SpatialExperiment object
#' @export



readVizgen <- function(data_dir,
                       filter_gene_pattern = "Blank",
                       tiff_path = NULL,
                       method_name = "Vizgen") {

    exprsMat <- data.table::fread(file.path(data_dir, "cell_by_gene.csv"))
    cell_idx <- exprsMat[, 1]
    exprsMat <- exprsMat[, -1]
    exprsMat <- as(as.matrix(t(exprsMat)), "dgCMatrix")
    colnames(exprsMat) <- paste("Cell", unlist(cell_idx), sep = "_")
    exprsMat <- exprsMat[!grepl(filter_gene_pattern, rownames(exprsMat)), ]

    metadata <- read.csv(file.path(data_dir, "cell_metadata.csv"))
    metadata <- metadata[order(metadata$X), ]
    metadata$cell_id <- metadata$X
    rownames(metadata) <- paste("Cell", metadata$cell_id, sep = "_")
    metadata <- metadata[colnames(exprsMat), ]


    spe <- SpatialExperiment::SpatialExperiment(assay = list(counts = exprsMat),
                                                colData = DataFrame(metadata),
                                                spatialCoordsNames = c("center_x",
                                                                       "center_y"))


    spe <- .initialise_CellSPA_list(spe)
    spe <- .add_dataset_metrics(spe)
    if (!is.null(tiff_path)) {
        tif_output <- readTiffOutput(tiff_path)
        tif_output <- tif_output[tif_output$cell_id %in% spe$cell_id, ]
        spe@metadata$CellSegOutput <- tif_output
        if (!all(spe$cell_id %in% spe@metadata$CellSegOutput$cell_id)) {
            message("There are cells that are not in the tiff file,\n
                    please check your tiff file")
            spe <- spe[, spe$cell_id %in% spe@metadata$CellSegOutput$cell_id]
        }
    }
    spe <- .cal_baseline(spe)
    spe@metadata$CellSPA$method <- method_name
    return(spe)

}

