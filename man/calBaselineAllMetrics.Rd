% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/baseline_metrics.R
\name{calBaselineAllMetrics}
\alias{calBaselineAllMetrics}
\title{Calculate the baseline metrics of cell from cell segmentation output}
\usage{
calBaselineAllMetrics(
  spe,
  metrics = c("total_transciprts", "total_genes", "total_cells", "meanExprsPct_cells",
    "cell_area", "elongation", "compactness", "eccentricity", "sphericity", "solidity",
    "convexity", "circularity", "density"),
  use_BPPARAM = BiocParallel::SerialParam(),
  rerun = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{spe}{A Spatial Experiment object}

\item{metrics}{A character vector indicating the metrics to be calculated}

\item{use_BPPARAM}{a BiocParallelParam instance determining the parallel
back-end to be used during evaluation.}

\item{rerun}{A logical value indicates whether to recalculate existing metrics.}

\item{verbose}{A logical value indicates whether to print out running messages.}
}
\value{
A SpatialExperiment object
}
\description{
Calculate the baseline metrics of cell from cell segmentation output
}
