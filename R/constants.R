.CellSPARenvir <- new.env(parent = emptyenv())
.CellSPARenvir[["cell_level_baseline_metrics"]] <- c("total_transciprts",
                                                     "total_genes",
                                                     "cell_area", "elongation",
                                                     "compactness", "eccentricity",
                                                     "sphericity", "solidity",
                                                     "convexity", "circularity")


.CellSPARenvir[["cell_shape_baseline_metrics"]] <- c("cell_area", "elongation",
                                                     "compactness", "eccentricity",
                                                     "sphericity", "solidity",
                                                     "convexity", "circularity")


.CellSPARenvir[["non_cell_shape_baseline_metrics"]] <- c("total_transciprts",
                                                         "total_genes",
                                                         "total_cells",
                                                         "meanExprsPct_cells")
