% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expression_metrics.R
\name{calMarkerPct}
\alias{calMarkerPct}
\title{Calculate marker expression purity for each cell given a marker list}
\usage{
calMarkerPct(
  spe,
  celltype,
  marker_list,
  marker_list_name = "",
  exprs_values = "logcounts",
  threshold = 0,
  use_BPPARAM = BiocParallel::SerialParam()
)
}
\arguments{
\item{spe}{A Spatial Experiment object}

\item{celltype}{A string or a vector indicating the cell type information.
If it is a string, it will retrive the information from colData.}

\item{marker_list}{A list of marker, the name of the list corresponding to
the cell type information provided in `celltype`.}

\item{marker_list_name}{A string specifying the name to be used to store the result.}

\item{exprs_values}{A string indicating which assay of `spe` contains the expression values.}

\item{threshold}{A numeric indicating the threshold to be considered as expressed.}

\item{use_BPPARAM}{A BiocParallelParam instance determining the parallel
back-end to be used during evaluation.}
}
\value{
A SpatialExperiment object
}
\description{
Calculate marker expression purity for each cell given a marker list
}
