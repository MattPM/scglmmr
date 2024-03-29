% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pseudobulk_helpers.r
\name{Normalize}
\alias{Normalize}
\title{Normalize - normalize summed gene counts across samples. This is a wrapper around edgeR functions `calcNormFactors` and `filterByExpr`.}
\usage{
Normalize(
  pseudobulk.list,
  design,
  group = NULL,
  normalization.method = "RLE",
  minimum.gene.count = 1
)
}
\arguments{
\item{pseudobulk.list}{object created with PseudobulkList only use this function if argument to PseudobulkList avg_or_sum was 'sum' computes normalization for pseudobulk libraries}

\item{design}{see edgeR `filterByExpr` function. This is the design matrix created using model.matrix and e.g. the metadata created  AggregateCellMetadata.}

\item{normalization.method}{see edgeR function calcNormFactors, this is the argument to `method`}

\item{minimum.gene.count}{see edgeR function `filterbyExpr` thie is the argument to `min.count`}
}
\value{
a list of dgeList objects indexed by celltype
}
\description{
Normalize - normalize summed gene counts across samples. This is a wrapper around edgeR functions `calcNormFactors` and `filterByExpr`.
}
\examples{
#' #'\dontrun{
}
}
